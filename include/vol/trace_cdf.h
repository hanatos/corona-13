#pragma once
typedef struct vol_cdf_entry_t
{
  float tmin;
  float tmax;
  float density;  // average density in this segment (assumed to be homogeneous)
  float emission; // average density * L_e in this segment (no sigma_e multiplied yet)
  float p;        // probability of this segment
}
vol_cdf_entry_t;

static inline int vol_cdf_create(
    vol_cdf_entry_t *cdf,
    const float sigma_t,
    const float sigma_e,
    const int cdf_len)
{
  // for each voxel, compute normalised probability p = tau * (sigma_e * emission)/(sigma_t * density) * (1-voxel transmittance) / #vox
  float transmittance = 1.0f;
  float sum = 0.0f;
  for(int i=0;i<cdf_len;i++)
  {
    float dist = cdf[i].tmax - cdf[i].tmin;
    float transmittance_voxel = expf(-cdf[i].density * sigma_t * dist);
    if(cdf[i].density > 0)
      cdf[i].p = transmittance * sigma_e * cdf[i].emission / (sigma_t * cdf[i].density) * (1.0f-transmittance_voxel) / i;
    else cdf[i].p = 0;
    transmittance *= transmittance_voxel;
    sum += cdf[i].p;
  }
  if(!(sum > 0.0f)) return 1;
  for(int i=0;i<cdf_len-1;i++)
  {
    cdf[i].p /= sum;
    cdf[i+1].p += cdf[i].p;
  }
  cdf[cdf_len-1] = 1.0f;
  return 0;
}

// returns pdf of given distance
static inline float vol_cdf_pdf(
    vol_cdf_entry_t *cdf,
    const int cdf_len,
    float dist)
{
  unsigned int min = 0, max = cdf_len;
  unsigned int t = max/2;
  while (t != min)
  {
    if(cdf[t].tmin <= dist) min = t;
    else max = t;
    t = (min + max)/2;
  }
  if(max < cdf_len && cdf[t].tmin <= dist) t = max;
  return 1.0f/(cdf[t].tmax - cdf[t].tmin) * (t ? cdf[t].p - cdf[t-1].p : cdf[0].p);
}

// returns the distance
static inline float vol_cdf_sample(
    vol_cdf_entry_t *cdf,
    const int cdf_len,
    const float rand,
    float *pdf)
{
  // sample
  unsigned int min = 0, max = cdf_len;
  unsigned int t = max/2;
  while (t != min)
  {
    if(cdf[t].p <= rand) min = t;
    else max = t;
    t = (min + max)/2;
  }
  // last step: decide between min and max one more time (min is rounding default),
  // but consider that if max is still out of bounds, it's invalid.
  // (rand == 1.0 and a cdf[0]=1, cdf_len=1 would break otherwise)
  if(max < cdf_len && cdf[t].p <= rand) t = max;

  // now sample uniform distance within voxel:
  float dist = 0.0f;
  if(t == 0) dist = cdf[0].tmin + (cdf[0].tmax - cdf[0].tmin) * rand/cdf[0].p;
  dist = cdf[t].tmin + (cdf[t].tmax - cdf[t].tmin) * (rand - cdf[t-1].p)/(cdf[t].p-cdf[t-1].p);
  assert(dist==dist);
  if(pdf) *pdf = 1.0f/(cdf[t].tmax - cdf[t].tmin) * (t ? cdf[t].p - cdf[t-1].p : cdf[0].p);
  return dist;
}

// returns the solid angle pdf of sampling this direction with NEE point sampling
static inline float vol_lighthierarchy_collect_cdf(
    const vol_tree_t *tree,              // volume tree
    const float *pos_ws,                 // world space position of query ray
    const float *dir_ws,                 // world space direction of query ray
    const float lambda,                  // wavelength
    float time,                          // time for motion blur access
    vol_interpolation_t interpolation,   // access with interpolation or motion blur?
    const int lod,                       // LOD: 0 is finest, 1 is 8x coarser, ..
    vol_cdf_entry_t *cdf,                // output will be written to here
    int *cdf_len,                        // return length of cdf here
    const int cdf_max_len)               // max length of output array
{
  assert(lod); // don't have aggregate nodes on leaf level
  const int use_max = 1;
  float accum_pdf = 0.0f;
  int cdflen = 0;
  *cdf_len = 0;
  float pos[3] = {pos_ws[0], pos_ws[1], pos_ws[2]};
  float dir[3] = {dir_ws[0], dir_ws[1], dir_ws[2]};
  vol_transform_w2o(tree, pos, 1); // transform to object space
  vol_transform_w2o(tree, dir, 0);
  for(int k=0;k<3;k++) // this is terrible:
    if(fabsf(dir[k]) < 1e-3f) dir[k] = copysignf(1e-3f, dir[k]);

  // intersect bounding box, get entry and exit distances tmin, tmax
  float tmin = 0.0f;
  float tmax = FLT_MAX;
  for(int k=0;k<3;k++)
  {
    const float t0 = (tree->content_box[k+3] - pos[k])/dir[k];
    const float t1 = (tree->content_box[k]   - pos[k])/dir[k];
    if(t0 <= t1)
    {
      tmin = t0 > tmin ? t0 : tmin;
      tmax = t1 < tmax ? t1 : tmax;
    }
    else
    {
      tmin = t1 > tmin ? t1 : tmin;
      tmax = t0 < tmax ? t0 : tmax;
    }
  }
  // miss bounding box?
  if(tmin > tmax) return 0.0f;
  const float tout = tmax;
  
  // used for path coding. stores decision for lower/upper block in i-th dimension.
  int x[3] = {0};

  // create xor mask from ray direction sign bits
  const int flip_x = (dir[0] < 0.0f)? 7 : 0;
  const int flip_y = (dir[1] < 0.0f)? 7 : 0;
  const int flip_z = (dir[2] < 0.0f)? 7 : 0;

  const int max_depth = 3*tree->header->depth;
  assert(max_depth < 20);
  // fine grained stack of pdfs and octree traversal things:
  float t0[20][3], tm[20][3], t1[20][3];
  float pdf[20] = {1.0f};
  float L[20][8];
  int sublevel[20];
  int isstatic[20] = {tree->header->isstatic};
  int32_t sp = 0;
  const vol_node_t *n[20];
  n[0] = tree->root_node;
  sublevel[0] = 0;
  const vol_payload_t *payload[20] = {tree->root_payload};
  const vol_payload_flux_t *node[20] = {tree->root_light};
  const vol_payload_aggregate_t *aggr[20] = {vol_node_get_aggregate(tree, tree->root_node)};

  for(int i=0;i<3;i++)
  {
    t0[0][i] = ((dir[i] > 0.0f ? tree->aabb[i]   : tree->aabb[3+i]) - pos[i])/dir[i];
    t1[0][i] = ((dir[i] > 0.0f ? tree->aabb[3+i] : tree->aabb[i])   - pos[i])/dir[i];
  }

  // voxel volume, adjusted for LOD/max level
  const float vs = tree->voxel_size * (1<<(3*lod));
  const float volume = vs*vs*vs; // is applied when returning

  // trace through the octree:
  while(1)
  {
    // get one from stack:
    tm[sp][0] = (t0[sp][0]+t1[sp][0])*.5;
    tm[sp][1] = (t0[sp][1]+t1[sp][1])*.5;
    tm[sp][2] = (t0[sp][2]+t1[sp][2])*.5;
    const float tx0 = t0[sp][0], txm = tm[sp][0];
    const float ty0 = t0[sp][1], tym = tm[sp][1];
    const float tz0 = t0[sp][1], tzm = tm[sp][2];

    // find entry plane:
    int max = tx0 > ty0 ? (tx0 > tz0 ? 0 : 2) : (ty0 > tz0 ? 1 : 2);
    switch(max)
    {
      case 2: // entry plane XY
        x[0] |= (txm < tz0) & 1;
        x[1] |= (tym < tz0) & 1;
        // tmin = MAX(0.0f, tz0);
        break;
      case 1: // entry plane XZ
        x[0] |= (txm < ty0) & 1;
        x[2] |= (tzm < ty0) & 1;
        // tmin = MAX(0.0f, ty0);
        break;
      default: // case 0: // entry plane YZ
        x[1] |= (tym < tx0) & 1;
        x[2] |= (tzm < tx0) & 1;
        // tmin = MAX(0.0, tx0);
        break;
    }

    // if depth is at a coarse level jump, get 8 values here and normalise pdf:
    if(sublevel[sp] == 2)
    {
      const int child =  // un-flip x coord
        (((x[2]^flip_z)&6)<<6)|(((x[1]^flip_y)&6)<<3)|((x[0]^flip_x)&6);
      for(int ii=0;ii<8;ii++)
      {
        const int idx = child | (((ii&4)<<4)|((ii&2)<<2)|(ii&1));
        float result[2] = {0.0f};
        vol_payload_get(payload[sp], idx, time, isstatic[sp], result);
        if(sp+1 < max_depth)
          L[sp][ii] = use_max ? result[0] : result[1];
        else if(result[0] > 0.0f && result[1] > 0.0f)
          L[sp][ii] = result[0] * tree->shader(result[0], result[1], lambda);
      }
      const float total = L[sp][0]+L[sp][1]+L[sp][2]+L[sp][3]+
                          L[sp][4]+L[sp][5]+L[sp][6]+L[sp][7];
      if(total > 0.0) for(int k=0;k<8;k++) L[sp][k] /= total;
      for(int k=0;k<8;k++)
        assert(L[sp][k] == L[sp][k]);
    }

    // step through child nodes on same level of current node:
    while(1)
    {
      const float tx1 = t1[sp][0], txm = tm[sp][0];
      const float ty1 = t1[sp][1], tym = tm[sp][1];
      const float tz1 = t1[sp][2], tzm = tm[sp][2];
      const float tmaxx = (x[0]&1) ? tx1 : txm, tmaxy = (x[1]&1) ? ty1 : tym, tmaxz = (x[2]&1) ? tz1 : tzm;

      max = 0; // dimension of closest exit plane
      if     (tmaxy <= tmaxx && tmaxy <= tmaxz) max = 1;
      else if(tmaxz <= tmaxx && tmaxz <= tmaxy) max = 2;
      assert(max != 0 || (tmaxx <= tmaxy && tmaxx <= tmaxz));
      tmax = MIN(tmaxx, MIN(tmaxy, tmaxz));

      // skip voxels behind the ray
      if(tmax > tmin)
      {
        const int i = (x[0]^flip_x)&1, j = (x[1]^flip_y)&1, k = (x[2]^flip_z)&1;
        if(sublevel[sp] == 0)
        {
          // get light hierarchy node->p[8][.] and push pdf to stack
          const uint16_t *data = use_max ? node[sp]->q[8] : node[sp]->p[8];
          pdf[sp+1] = pdf[sp];
          if(k) pdf[sp+1] *= 1.0f-data[0      ]/(float)0x10000;
          else  pdf[sp+1] *=      data[0      ]/(float)0x10000;
          if(j) pdf[sp+1] *= 1.0f-data[1+k    ]/(float)0x10000;
          else  pdf[sp+1] *=      data[1+k    ]/(float)0x10000;
          if(i) pdf[sp+1] *= 1.0f-data[3+2*k+j]/(float)0x10000;
          else  pdf[sp+1] *=      data[3+2*k+j]/(float)0x10000;
          assert(pdf[sp+1] == pdf[sp+1]);
          if(pdf[sp+1] > 0.0f)
            break; // descend to next finer level
        }
        else if(sublevel[sp] == 1)
        {
          // get light hierarchy node->p[c][.] and push pdf to stack
          const int c0 = // un-flip x coord
            (((x[2]^flip_z)&2)<<1)|((x[1]^flip_y)&2)|(((x[0]^flip_x)&2)>>1);
          assert(c0 >= 0);
          assert(c0 < 8);
          const uint16_t *data = use_max ? node[sp]->q[c0] : node[sp]->p[c0];
          pdf[sp+1] = pdf[sp];
          if(k) pdf[sp+1] *= 1.0f-data[0      ]/(float)0x10000;
          else  pdf[sp+1] *=      data[0      ]/(float)0x10000;
          if(j) pdf[sp+1] *= 1.0f-data[1+k    ]/(float)0x10000;
          else  pdf[sp+1] *=      data[1+k    ]/(float)0x10000;
          if(i) pdf[sp+1] *= 1.0f-data[3+2*k+j]/(float)0x10000;
          else  pdf[sp+1] *=      data[3+2*k+j]/(float)0x10000;
          assert(pdf[sp+1] == pdf[sp+1]);
          if(pdf[sp+1] > 0.0f)
            break; // next finer level
        }
        else if(sublevel[sp] == 2)
        { // collect pdf with tmin and tmax:
          const int c = (k<<2)|(j<<1)|i;
          const int child =
            (((x[2]&7)^flip_z)<<6)|(((x[1]&7)^flip_y)<<3)|((x[0]&7)^flip_x);
          if(sp+3*lod+1 >= max_depth)
          { // leaf level/lod
            assert(accum_pdf == accum_pdf);
            assert(pdf[sp] == pdf[sp]);
            assert(L[sp][c] == L[sp][c]);
            accum_pdf += (tmax*tmax*tmax-tmin*tmin*tmin)/3.0f * pdf[sp] * L[sp][c];
            assert(accum_pdf == accum_pdf);
            cdf[cdflen].tmin = tmin;
            cdf[cdflen].tmax = tmax;
            // if lod==0, we would need to ask the payload here
            cdf[cdflen].density = aggr[sp]->density[child];
            cdf[cdflen].emission = L[sp][c];
            cdflen++;
            assert(cdflen <= cdf_max_len);
          }
          else
          {
            if(!vol_node_child_empty(n[sp], child))
            {
              pdf[sp+1] = pdf[sp] * L[sp][c];
              if(pdf[sp+1] > 0.0f)
                break; // descend further
            }
          }
        }
        else assert(0);
        tmin = tmax; // step over this voxel
      }

      // step out of child node through closest exit plane
      uint32_t oldx = x[max]++;
      if(tmax >= tout) // stepped out of content box
      {
        *cdf_len = cdflen;
        return accum_pdf/volume;
      }

      // pop stack, go up one level in the octree:
      while((x[max] ^ oldx) > 1)
      {
        sp--;
        for(int i=0;i<3;i++) x[i] >>= 1;
        oldx >>= 1;
      }
    }
    sp++;
    sublevel[sp] = sublevel[sp-1]+1;
    if(sublevel[sp] == 3)
    { // grid level jump
      const int child =
        (((x[2]&7)^flip_z)<<6)|(((x[1]&7)^flip_y)<<3)|((x[0]&7)^flip_x);
      payload[sp] = vol_node_get_payload(tree, n[sp-1], child);
      isstatic[sp] = vol_node_child_static(n[sp-1], child);
      assert(payload[sp]);
      n[sp] = vol_node_get_child(tree, n[sp-1], child);
      node[sp] = vol_node_get_flux(tree, n[sp-1], child);
      aggr[sp] = vol_node_get_aggregate(tree, n[sp]);
      sublevel[sp] = 0;
    }
    else
    { // stay on same grid level
      n[sp] = n[sp-1];
      payload[sp] = payload[sp-1];
      node[sp] = node[sp-1];
      isstatic[sp] = isstatic[sp-1];
    }

    // copy (t0[sp], t1[sp]) <- (t0[sp-1], t1[sp-1], tm[sp-1])
    for(int i=0;i<3;i++)
    {
      t0[sp][i] = (x[i] & 1) ? tm[sp-1][i] : t0[sp-1][i];
      t1[sp][i] = (x[i] & 1) ? t1[sp-1][i] : tm[sp-1][i];
      x[i] <<= 1;
    }
  }
  assert(0);
  return 0.0f; // never reached
}
