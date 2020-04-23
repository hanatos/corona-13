  // setup ray traversal:
  // intersect content box, if no intersection, return
  // it is sufficient to get tmin/tmax from content_box and do all the rest with the voxel aabb
  const float *cbox = tree->content_box;
  float invdir[3];
  for(int k=0;k<3;k++) invdir[k] = 1.0f/dir[k];
  const int step[3] = {invdir[0] > 0 ? 1 : -1, invdir[1] > 0 ? 1 : -1, invdir[2] > 0 ? 1 : -1 };
  // stack with stuff for nested grids on each level:
  const int md = tree->header->depth;
  const int leaf_d = md - 1 > lod ? md-1-lod : 0;
  float aabb[md][6], t[md][3], delta[md][3];
  const vol_node_t *node[md];
  const vol_payload_t *data[md];
  int isstatic[md];
  int mind[md], p[md][3];
  int sp = 0;   // stack pointer
  for(int k=0;k<6;k++) aabb[0][k] = tree->aabb[k];
  // node[0] = leaf_d ? tree->root_node : 0; // can't do this if doing stochastic tracking
  node[0] = tree->root_node;
  data[0] = tree->root_payload;
  isstatic[0] = tree->header->isstatic;

  // size of voxel in world space:
  float vwd[md];
  vwd[md-1] = tree->voxel_size;
  for(int d=md-2;d>=0;d--) vwd[d] = vwd[d+1]*8;

  float pos2[3] = {pos[0], pos[1], pos[2]};
  if(interpolation & s_vol_smooth) vol_interpolate_smooth(&pos2[0], &pos2[1], &pos2[2], &time, vwd[leaf_d], 1.0/VOL_MOTION_SAMPLES);

  // intersect bounding box, get entry and exit distances tmin, tmax
  float tmin = 0.0f;
  float tout = max_dist;
  for(int k=0;k<3;k++)
  {
    const float t0 = (cbox[k+3] - pos[k])*invdir[k];
    const float t1 = (cbox[k]   - pos[k])*invdir[k];
    if(t0 <= t1)
    {
      tmin = t0 > tmin ? t0 : tmin;
      tout = t1 < tout ? t1 : tout;
    }
    else
    {
      tmin = t1 > tmin ? t1 : tmin;
      tout = t0 < tout ? t0 : tout;
    }
  }
  // miss bounding box?
  if(tmin > tout)
  {
    VOL_TRACE_MISS
  }

  // step delta for t[] in each dimension (could only init up to <= leaf_d)
  for(int d=0;d<md;d++) for(int k=0;k<3;k++)
    delta[d][k] = vwd[d]*fabsf(invdir[k]);

  // TODO OPT: this needs to be done every time we step into a new grid!
  // TODO OPT: find out when exactly (only when stepping down, and only if not adjacent grid of same size..)
  // TODO OPT: and optimise!
  // entry point as grid index p[]:
  // init t[] to distance offsets to next voxel in every dimension by intersecting a voxel-sized aabb:
  for(int k=0;k<3;k++)
  {
    p[sp][k] = (int)CLAMP((pos2[k] + tmin*dir[k] - aabb[sp][k])/vwd[sp], 0, 7);
    t[sp][k] = (aabb[sp][k] + (p[sp][k] + ((step[k]==1)&1)) * vwd[sp] - pos2[k])*invdir[k];
  }

  VOL_TRACE_INIT
  while(1)
  {
    // compute entry point distance t to next voxel:
    if(t[sp][2] < t[sp][1] && t[sp][2] < t[sp][0]) mind[sp] = 2;
    else if(t[sp][1] < t[sp][0]) mind[sp] = 1;
    else mind[sp] = 0;
    // assert(t[sp][mind[sp]] >= tmin); // when entering a new nested grid, this sometimes happens.
    // hence the CLAMP() clamping when computing the density below.

    // step one voxel: tmin..t[mind]:
    // this can mean descending into a finer grid level:
    vol_index_t idx = {0}; idx.i = p[sp][0]; idx.j = p[sp][1]; idx.k = p[sp][2];
    const float tmax = t[sp][mind[sp]];
    assert(idx.idx >= 0 && idx.idx < 512);
    // no speed difference: if(sp < md-1)
    if(node[sp] && sp < leaf_d)
    { // inner level
      if(!vol_node_child_empty(node[sp], idx.idx))
      { // only descent if not empty
        sp++;
        isstatic[sp] = vol_node_child_static(node[sp-1], idx.idx);
        if(sp < md - 1 )//leaf_d) // also init for = leaf_d so we can access it in LEAF part
          node[sp] = vol_node_get_child(tree, node[sp-1], idx.idx);
        else node[sp] = 0;
        data[sp] = vol_node_get_payload(tree, node[sp-1], idx.idx);
        // TODO OPT: optimise this!
        for(int k=0;k<3;k++)
        {
          aabb[sp][k+0] = aabb[sp-1][k] +  p[sp-1][k]    * (aabb[sp-1][k+3]-aabb[sp-1][k])/8.0f;
          aabb[sp][k+3] = aabb[sp-1][k] + (p[sp-1][k]+1) * (aabb[sp-1][k+3]-aabb[sp-1][k])/8.0f;
          p[sp][k] = (int)CLAMP((pos2[k] + tmin*dir[k] - aabb[sp][k])/vwd[sp], 0, 7);
          t[sp][k] = (aabb[sp][k] + (p[sp][k] + ((step[k]==1)&1)) * vwd[sp] - pos2[k])*invdir[k];
        }
        continue; // restart grid traversal
      }
      tmin = t[sp][mind[sp]];
      t[sp][mind[sp]] += delta[sp][mind[sp]];
      p[sp][mind[sp]] += step[mind[sp]];
    }
    else // node[sp] == 0 => leaf and only payload data is valid
    {    // leaf level or LOD reached
      float result[2] = {0}; // for trilinear interpolation
      vol_payload_get(data[sp], idx.idx, time, isstatic[sp], result);
      VOL_TRACE_LEAF
      // jump to next leaf voxel
      tmin = tmax;
      t[sp][mind[sp]] += delta[sp][mind[sp]];
      p[sp][mind[sp]] += step[mind[sp]];
    }

    while(p[sp][mind[sp]] >= 8 || p[sp][mind[sp]] < 0)
    {
      sp--;
      if(sp < 0 || tmin >= tout)
      {
        VOL_TRACE_RETURN
      }

      // returned to parent grid, now jump over the big voxel we just processed
      tmin = t[sp][mind[sp]];
      t[sp][mind[sp]] += delta[sp][mind[sp]];
      p[sp][mind[sp]] += step[mind[sp]];
    }
    if(tmin >= tout)
    {
      VOL_TRACE_RETURN
    }
  }
#undef VOL_TRACE_INIT
#undef VOL_TRACE_LEAF
#undef VOL_TRACE_MISS
#undef VOL_TRACE_RETURN
