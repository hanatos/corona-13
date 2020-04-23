static inline mf_t vol_lighthierarchy_pdf_direction(
    const vol_tree_t *tree,              // volume tree
    const float *pos_ws,                 // world space position of query ray
    const float *dir_ws,                 // world space direction of query ray
    const mf_t lambda,                   // wavelength
    float time,                          // time for motion blur access
    vol_interpolation_t interpolation,   // access with interpolation or motion blur?
    const int lod,                       // LOD: 0 is finest, 1 is 8x coarser, ..
    const int use_max,
    const float *cullx_ws,
    const float *culln_ws)
{
  float accum_pdf = 0.0f;
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
  if(tmin > tmax) return mf_set1(0.0f);
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

  for(int i=0;i<3;i++)
  {
    t0[0][i] = ((dir[i] > 0.0f ? tree->aabb[i]   : tree->aabb[3+i]) - pos[i])/dir[i];
    t1[0][i] = ((dir[i] > 0.0f ? tree->aabb[3+i] : tree->aabb[i])   - pos[i])/dir[i];
  }

  // voxel volume, adjusted for LOD/max level
  const float vs = tree->voxel_size * (1<<(3*lod));
  const float volume = vs*vs*vs; // is applied when returning

  float cullv[8];
  float binary[7];
  float cullxs[3], *cullx = 0;
  float cullns[3], *culln = 0;
  if(cullx_ws && culln_ws)
  {
    for(int k=0;k<3;k++) cullxs[k] = cullx_ws[k];
    for(int k=0;k<3;k++) cullns[k] = culln_ws[k];
    vol_transform_o2w(tree, cullxs, 1); // transform object to world space
    vol_transform_o2w(tree, cullns, 0); // det == 1, so transforming normals like that is okay
    for(int k=0;k<3;k++)
      cullxs[k] -= tree->aabb[k];
    cullx = cullxs;
    culln = cullns;
  }

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
      const int mask = ~-(1<<sp);
      _eval_culln(cullx, culln,
          flip_x ? (x[0]/2) ^ mask : (x[0]/2),
          flip_y ? (x[1]/2) ^ mask : (x[1]/2),
          flip_z ? (x[2]/2) ^ mask : (x[2]/2),
          tree->voxel_size, max_depth-sp, cullv);
      for(int ii=0;ii<8;ii++)
      {
        const int idx = child | (((ii&4)<<4)|((ii&2)<<2)|(ii&1));
        float result[2] = {0.0f};
        vol_payload_get(payload[sp], idx, time, isstatic[sp], result);
        if(sp+1 < max_depth)
          L[sp][ii] = cullv[ii]*(use_max ? result[0] : result[1]);
        else if(result[0] > 0.0f && result[1] > 0.0f)
          // XXX FIXME: actually should all be spectral here!
          L[sp][ii] = cullv[ii] * result[0] * mf(tree->shader(result[0], result[1], lambda),0);
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
          // _eval_culln(cullx, culln, x[0]/2, x[1]/2, x[2]/2, tree->voxel_size, max_depth-sp, cullv);
          const int mask = ~-(1<<sp);
          _eval_culln(cullx, culln,
              flip_x ? (x[0]/2) ^ mask : (x[0]/2),
              flip_y ? (x[1]/2) ^ mask : (x[1]/2),
              flip_z ? (x[2]/2) ^ mask : (x[2]/2),
              tree->voxel_size, max_depth-sp, cullv);
          _eval_cullp(cullv, data, binary);
          pdf[sp+1] = pdf[sp];
          if(k) pdf[sp+1] *= 1.0f-binary[0      ];//data[0      ]/(float)0x10000;
          else  pdf[sp+1] *=      binary[0      ];//data[0      ]/(float)0x10000;
          if(j) pdf[sp+1] *= 1.0f-binary[1+k    ];//data[1+k    ]/(float)0x10000;
          else  pdf[sp+1] *=      binary[1+k    ];//data[1+k    ]/(float)0x10000;
          if(i) pdf[sp+1] *= 1.0f-binary[3+2*k+j];//data[3+2*k+j]/(float)0x10000;
          else  pdf[sp+1] *=      binary[3+2*k+j];//data[3+2*k+j]/(float)0x10000;
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
          //_eval_culln(cullx, culln, x[0]/2, x[1]/2, x[2]/2, tree->voxel_size, max_depth-sp, cullv);
          const int mask = ~-(1<<sp);
          _eval_culln(cullx, culln,
              flip_x ? (x[0]/2) ^ mask : (x[0]/2),
              flip_y ? (x[1]/2) ^ mask : (x[1]/2),
              flip_z ? (x[2]/2) ^ mask : (x[2]/2),
              tree->voxel_size, max_depth-sp, cullv);
          _eval_cullp(cullv, data, binary);
          if(k) pdf[sp+1] *= 1.0f-binary[0      ];//data[0      ]/(float)0x10000;
          else  pdf[sp+1] *=      binary[0      ];//data[0      ]/(float)0x10000;
          if(j) pdf[sp+1] *= 1.0f-binary[1+k    ];//data[1+k    ]/(float)0x10000;
          else  pdf[sp+1] *=      binary[1+k    ];//data[1+k    ]/(float)0x10000;
          if(i) pdf[sp+1] *= 1.0f-binary[3+2*k+j];//data[3+2*k+j]/(float)0x10000;
          else  pdf[sp+1] *=      binary[3+2*k+j];//data[3+2*k+j]/(float)0x10000;
          assert(pdf[sp+1] == pdf[sp+1]);
          if(pdf[sp+1] > 0.0f)
            break; // next finer level
        }
        else if(sublevel[sp] == 2)
        { // collect pdf with tmin and tmax:
          const int c = (k<<2)|(j<<1)|i;
          if(sp+3*lod+1 >= max_depth)
          { // leaf level/lod
            assert(accum_pdf == accum_pdf);
            assert(pdf[sp] == pdf[sp]);
            assert(L[sp][c] == L[sp][c]);
            accum_pdf += (tmax*tmax*tmax-tmin*tmin*tmin)/3.0f * pdf[sp] * L[sp][c];
            assert(accum_pdf == accum_pdf);
          }
          else
          {
            const int child =
              (((x[2]&7)^flip_z)<<6)|(((x[1]&7)^flip_y)<<3)|((x[0]&7)^flip_x);
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
        return mf_set1(accum_pdf/volume);

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
  return mf_set1(0.0f); // never reached
}
