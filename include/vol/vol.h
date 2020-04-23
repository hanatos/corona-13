#pragma once

#include "vol/types.h"
#include "vol/lighthierarchy_create.h"
#include "vol/shaders.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <assert.h>
#include <inttypes.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

static inline int vol_node_child_empty(const vol_node_t *node, const int i)
{
  // 255 signifies empty, in all but the last nodes in each block (where 255 is actually a valid value).
  if(i == 255) return node->off255_empty;
  if(i == 511) return node->off511_empty;
  else return node->off[i] == 255;
}

static inline int vol_node_child_static(const vol_node_t *node, const int i)
{
  if(i < 256) return node->data_static0;
  else        return node->data_static1;
}

static inline int vol_node_leaf(const vol_node_t *node)
{
  return node->node_leaf;
}

static inline vol_node_t *vol_node_get_child(const vol_tree_t *tree, const vol_node_t *node, const int i)
{
  if(vol_node_child_empty(node, i)) return 0;
  if(vol_node_leaf(node)) return 0;
  if(i > 255) return tree->nodes + node->node_offset1 + node->off[i];
  else        return tree->nodes + node->node_offset0 + node->off[i];
}

static inline vol_payload_t *vol_node_get_payload(const vol_tree_t *tree, const vol_node_t *node, const int i)
{
  if(vol_node_child_empty(node, i)) return 0;
  if(i > 255)
  {
    size_t psize = node->data_static1 ? vol_payload_static_size() : sizeof(vol_payload_t);
    return (vol_payload_t *)(tree->payload + node->data_offset1 + psize * node->off[i]);
  }
  else
  {
    size_t psize = node->data_static0 ? vol_payload_static_size() : sizeof(vol_payload_t);
    return (vol_payload_t *)(tree->payload + node->data_offset0 + psize * node->off[i]);
  }
}

static inline vol_payload_flux_t *vol_node_get_flux(
    const vol_tree_t *tree, const vol_node_t *node, const int i)
{
  if(vol_node_child_empty(node, i)) return 0;
  if(i > 255) return tree->light + node->lh_offset1 + node->off[i];
  else        return tree->light + node->lh_offset0 + node->off[i];
}

static inline void vol_tree_get_node_aabb(const vol_tree_t *tree, const vol_index_t *i, const int level, float *aabb)
{
  memcpy(aabb, tree->aabb, sizeof(float)*6);
  for(int l=1;l<level;l++)
  {
    const float wd = aabb[3] - aabb[0];
    const float o0 = aabb[0];
    aabb[0] = o0 +  i[l].i    * wd/8.0f;
    aabb[3] = o0 + (i[l].i+1) * wd/8.0f;
    const float ht = aabb[4] - aabb[1];
    const float o1 = aabb[1];
    aabb[1] = o1 +  i[l].j    * ht/8.0f;
    aabb[4] = o1 + (i[l].j+1) * ht/8.0f;
    const float dp = aabb[5] - aabb[2];
    const float o2 = aabb[2];
    aabb[2] = o2 +  i[l].k    * dp/8.0f;
    aabb[5] = o2 + (i[l].k+1) * dp/8.0f;
  }
}

static inline vol_payload_type_t vol_create_leaf(
    vol_payload_uncompressed_t *payload,
    float aabb[6],
    int force_static,
    vol_payload_fill_t payload_fill,
    void *data)
{
  memset(payload, 0, sizeof(*payload));
  return payload_fill(data, payload, aabb, force_static);
}

static inline vol_payload_type_t vol_create_subtree(
    vol_tree_t *tree,
    vol_node_t *node,
    vol_payload_uncompressed_t *payload,
    vol_index_t *const child,
    const int depth,
    const int max_depth,
    int fd_nodes,
    vol_payload_fill_t payload_fill,
    void *fill_data,
    const int force_static,
    int fd_light)
{
  // early out empty test:
  float aabb[6];
  vol_tree_get_node_aabb(tree, child, depth, aabb);
  for(int k=0;k<3;k++) if(aabb[k+3] < tree->header->content_box[k+0] ||
                          aabb[k+0] > tree->header->content_box[k+3]) return s_vol_empty;

  // TODO: root node is a leaf special case:
  // TODO: set all child pointers in node to `empty' and init payload as a leaf
  // TODO: not sure we can achieve this.

  memset(node, 0, sizeof(*node));
  vol_payload_type_t fill = s_vol_static;

  vol_node_t *cnodes = 0;
  const int leaf_level = (depth + 1 == max_depth);
  if(leaf_level)
    // mark as leaf in parent: only payload data pointers are valid, no more tree structure.
    node->node_leaf = 1;
  else cnodes = (vol_node_t *)malloc(sizeof(vol_node_t)*512);

  // could alloc once for this tree depth and reuse, but whatever.
  vol_payload_uncompressed_t *cdata = (vol_payload_uncompressed_t *)malloc(sizeof(vol_payload_uncompressed_t)*512);
  memset(cdata, 0, sizeof(*cdata)*512);
  vol_payload_t *compressed = (vol_payload_t *)malloc(sizeof(vol_payload_t)*512);

  vol_index_t vx = {0};
  // 1) parallelise leaf level / go through inner nodes in serial
  vol_payload_type_t empty[512] = {s_vol_empty};
  if(leaf_level)
  {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) default(none) shared(aabb, tree, empty, fill_data, payload_fill, cdata)
#endif
    for(int i=0;i<512;i++)
    { // handle leaf nodes in parallel (sample/compress)
      vol_index_t vx = {0};
      vx.idx = i;
      vol_index_t nchild[max_depth+1];
      memcpy(nchild, child, sizeof(vol_index_t)*depth);
      nchild[depth] = vx;

      // create only leaf payload data
      float caabb[6];
      vol_tree_get_node_aabb(tree, nchild, depth+1, caabb);
      empty[vx.idx] = vol_create_leaf(cdata+vx.idx, caabb, force_static, payload_fill, fill_data);
    }
  }
  else for(vx.idx=0;vx.idx<512;vx.idx++)
  { // else handle inner nodes, serial
    if(depth == 1) fprintf(stderr, "\r[vol progress] %0.2f %%    ", 100.0 * vx.idx/512.0);
    child[depth] = vx;
    empty[vx.idx] = vol_create_subtree(tree, cnodes+vx.idx, cdata+vx.idx, child, depth + 1, max_depth, fd_nodes, payload_fill, fill_data, force_static, fd_light);
  }

  // 2) compress in parallel:
  node->data_static0 = 1;
  node->data_static1 = 1;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) default(none) shared(empty, cdata, compressed, node, fill)
#endif
  for(int i=0;i<512;i++)
  {
    if(empty[i] != s_vol_empty)
    {
      const int s = vol_payload_compress(cdata+i, compressed+i, empty[i] == s_vol_static);
      if(!s)
      {
        if(i < 256) node->data_static0 = 0;
        else        node->data_static1 = 0;
        fill = s_vol_full;
      }
    }
  }

  int cnt = 0;
  // 3) compactify offsets
  for(int i=0;i<512;i++)
  {
    if(i == 256) cnt = 0;
    if(empty[i] == s_vol_empty)
    {
      if(i == 255) node->off255_empty = 1;
      if(i == 511) node->off511_empty = 1;
      node->off[i] = 255;
    }
    else
    {
      // store offset to child:
      node->off[i] = cnt++;
    }
  }

  // create light hierarchy aggregates. we do this the memory expensive way, to
  // be able to parallelise it:
  vol_payload_flux_t light[512] = {{{{0}}}}; // bracket overkill
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) default(shared)
#endif
  for(int i=0;i<512;i++)
  {
    if(!vol_node_child_empty(node, i))
      // XXX need to pass child_static and force_static, really:
      vol_lighthierarchy_create_node(payload, cdata+i, i, light+i,
          leaf_level, tree->shader, vol_node_child_static(node, i));
  }

  // now all children have written their stuff, it's our turn to dump the memory.
  cnt = 0;
  node->lh_offset0 = lseek(fd_light, 0, SEEK_END)/sizeof(vol_payload_flux_t);
  node->data_offset0 = lseek(tree->fd, 0, SEEK_END)-sizeof(vol_header_t);
  node->node_offset0 = lseek(fd_nodes, 0, SEEK_END)/sizeof(vol_node_t);
  for(int i=0;i<512;i++)
  {
    if(i == 256)
    {
      node->lh_offset1 = lseek(fd_light, 0, SEEK_END)/sizeof(vol_payload_flux_t);
      node->data_offset1 = lseek(tree->fd, 0, SEEK_END)-sizeof(vol_header_t);
      node->node_offset1 = lseek(fd_nodes, 0, SEEK_END)/sizeof(vol_node_t);
    }
    if(!vol_node_child_empty(node, i))
    {
      cnt++;
      // write out light hierarchy file
      if(write(fd_light, light+i, sizeof(light[i])) < (ssize_t)sizeof(light[i])) goto fail;

      if(!vol_node_leaf(node) &&
          write(fd_nodes, cnodes+i, sizeof(vol_node_t)) < (ssize_t)sizeof(vol_node_t))
      {
        fprintf(stderr, "[vol_create_subtree] writing child node [%d] to disk failed!\n", i);
        goto fail;
      }
      const size_t psize = vol_node_child_static(node, i) ?
        vol_payload_static_size() : sizeof(vol_payload_t);
      if(write(tree->fd, compressed+i, psize) < (ssize_t)psize)
      {
        fprintf(stderr, "[vol_create_subtree] writing child data [%d] to disk failed!\n", i);
        goto fail;
      }
    }
  }
  free(cnodes);
  free(cdata);
  free(compressed);

  // all empty
  if(cnt == 0) return s_vol_empty;
  return fill; // return static or full payloads

fail:
  free(cnodes);
  free(cdata);
  free(compressed);
  return s_vol_empty;
}


static inline int vol_create_tree(
    const char *filename,
    const float *aabb,
    const float voxel_size,
    const float *loc,
    const float *rot,
    vol_payload_fill_t payload_fill,
    void *fill_data,
    const vol_shader_id_t shaderid,
    int force_static)
{
  vol_header_t header;
  memset(&header, 0, sizeof(header));
  header.magic = VOL_HEADER_MAGIC;
  header.version = VOL_HEADER_VERSION | (VOL_MOTION_SAMPLES << 16);
  memcpy(header.content_box, aabb, sizeof(float)*6);
  header.voxel_size = voxel_size;
  vol_emission_shader_t shader = vol_shader_get(shaderid);
  // init identity transform or copy what we are passed in:
  for(int k=0;k<3;k++) header.rot[k] = header.loc[k] = 0.0f;
  if(loc) for(int k=0;k<3;k++) header.loc[k] = loc[k];
  if(rot) for(int k=0;k<3;k++) header.rot[k] = rot[k];

  // now find depth such that resulting aabb and square voxels cover full content.
  const int max_depth = 8;
  int depth = 2; // depth == 1 means root node is a leaf, we don't currently support this.
  for(;depth<=max_depth;depth++)
  {
    // size of tree bounding box given size of smallest voxel and tree depth
    const float sz = voxel_size * powf(8, depth);
    // if that covers the whole content box, we're good.
    if(sz >= aabb[3]-aabb[0] && sz >= aabb[4]-aabb[1] && sz >= aabb[5]-aabb[2])
      break;
  }
  header.depth = depth;
  for(int k=0;k<3;k++) header.aabb[k]   = header.content_box[k];
  for(int k=0;k<3;k++) header.aabb[3+k] = header.content_box[k] + voxel_size * powf(8, depth);
  const double sizeMB = sizeof(vol_payload_t)*(1-powf(512, depth))/(1.0-512.0)/1024.0/1024.0;
  if(sizeMB < 1024.0)
    fprintf(stderr, "[vol_create_tree] creating tree (depth %d), may need %.02f MB.\n",
        depth, sizeMB);
  else if(sizeMB < 1024.0 * 1024.0)
    fprintf(stderr, "[vol_create_tree] creating tree (depth %d), may need %.02f GB.\n",
        depth, sizeMB/1024.0);
  else // hope this is is..
    fprintf(stderr, "[vol_create_tree] creating tree (depth %d), may need %.02f TB.\n",
        depth, sizeMB/1024.0/1024.0);

  // call tree construction recursively. for that, init a tree_t struct:
  vol_tree_t tree;
  memcpy(tree.aabb, header.aabb, sizeof(float)*6);
  memcpy(tree.content_box, header.content_box, sizeof(float)*6);
  tree.voxel_size = header.voxel_size;
  tree.header = &header;

  tree.fd = open(filename, O_CREAT|O_WRONLY|O_TRUNC, 0644);
  // write to file
  if(tree.fd < 0 || write(tree.fd, &header, sizeof(header)) < 0)
  {
    fprintf(stderr, "[vol_create_tree] failed to open data file `%s'!\n", filename);
    return 1;
  }

  char tempnam[] = "/tmp/nodesXXXXXX";
  int fd_nodes = mkstemp(tempnam);
  if(fd_nodes < 0)
  {
    fprintf(stderr, "[vol_create_tree] failed to open node file `%s'!\n", tempnam);
    close(tree.fd);
    return 1;
  }
  unlink(tempnam);

  char lhfilename[] = "/tmp/lightXXXXXX";
  int fd_light = mkstemp(lhfilename);
  // write to file
  if(fd_light < 0)
  {
    fprintf(stderr, "[vol_create_tree] failed to open light hierarchy data file `%s'!\n", lhfilename);
    return 1;
  }
  unlink(lhfilename);

  // storage for child indices on stack
  vol_index_t child[depth];
  memset(child, 0, sizeof(child));

  vol_node_t root;
  vol_payload_uncompressed_t upayload = {{{0.0f}}};
  vol_payload_flux_t light = {{{0}}};
  header.shaderid = shaderid;
  tree.shader = vol_shader_get(shaderid);
  header.isstatic = s_vol_static == vol_create_subtree(&tree, &root, &upayload,
      child, 1, depth, fd_nodes, payload_fill, fill_data, force_static, fd_light);
  fprintf(stderr, "\r[vol progress] 100.00 %%    \n");

  // write out root node and payload data
  vol_payload_t payload;
  header.isstatic = vol_payload_compress(&upayload, &payload, header.isstatic);
  const size_t root_payload_size = header.isstatic ? vol_payload_static_size() : sizeof(payload);
  ssize_t wr = write(tree.fd, &payload, root_payload_size);
  wr += write(fd_nodes, &root, sizeof(root));
  if(wr != (ssize_t)(root_payload_size + sizeof(root))) goto fail;

  vol_lighthierarchy_create_node(0, &upayload, 0, &light, depth==1, shader, force_static);
  if(write(fd_light, &light, sizeof(light)) < (ssize_t)sizeof(light)) goto fail;

  // remember where in the file the nodes start:
  header.nodes = lseek(tree.fd, 0, SEEK_END);

  // append temp files to tree file
  {
  uint64_t buf[1024];
  ssize_t remaining = lseek(fd_nodes, 0, SEEK_END);
  lseek(fd_nodes, 0, SEEK_SET); // rewind
  while(remaining > 0)
  {
    ssize_t s = MIN((ssize_t)sizeof(buf), remaining);
    s = read(fd_nodes, buf, s);
    wr = write(tree.fd, buf, s);
    if(wr != s) goto fail;
    remaining -= s;
  }
  header.light = lseek(tree.fd, 0, SEEK_END);
  remaining = lseek(fd_light, 0, SEEK_END);
  lseek(fd_light, 0, SEEK_SET); // rewind
  while(remaining > 0)
  {
    ssize_t s = MIN((ssize_t)sizeof(buf), remaining);
    s = read(fd_light, buf, s);
    wr = write(tree.fd, buf, s);
    if(wr != s) goto fail;
    remaining-= s;
  }
  header.end = lseek(tree.fd, 0, SEEK_END);
  }

  // overwrite first 4k page on disk with useful information (offsets to first level of inner nodes + metadata)
  lseek(tree.fd, 0, SEEK_SET);
  if(write(tree.fd, &header, sizeof(header)) != sizeof(header)) goto fail;

  close(tree.fd);
  close(fd_nodes);
  close(fd_light);
  return 0;
fail:
  fprintf(stderr, "[vol_create_tree] failed writing to `%s'!\n", filename);
  close(tree.fd);
  close(fd_nodes);
  close(fd_light);
  return 1;
}

static inline void vol_close(vol_tree_t *tree)
{
  munmap(tree->data, tree->data_size);
  close(tree->fd);
  free(tree);
}

static inline void vol_create_transform(vol_tree_t *tree, const float *loc, const float *rot)
{
  memset(tree->w2o, 0, sizeof(float)*12);
  tree->w2o[0][0] = tree->w2o[1][1] = tree->w2o[2][2] = 1.0f;
  if(rot)
  {
    const float sx = sinf(rot[0]), cx = cosf(rot[0]);
    const float sy = sinf(rot[1]), cy = cosf(rot[1]);
    const float sz = sinf(rot[2]), cz = cosf(rot[2]);
    // TODO: check signs (esp y rot and translation) and order of execution!
    tree->w2o[0][0] = cy*cz;
    tree->w2o[0][1] = -cx*sz - sx*sy*cz;
    tree->w2o[0][2] = sx*sz - cx*sy*cz;
    tree->w2o[1][0] = cy*sz;
    tree->w2o[1][1] = cx*cz - sx*sy*sz;
    tree->w2o[1][2] = -sx*cz - cx*sy*sz;
    tree->w2o[2][0] = sy;
    tree->w2o[2][1] = -sx*cy;
    tree->w2o[2][2] = cx*cy;
  }
  if(loc)
  {
    tree->w2o[0][3] = loc[0];
    tree->w2o[1][3] = loc[1];
    tree->w2o[2][3] = loc[2];
  }
  tree->transform = (rot && ((rot[0] != 0) || (rot[1] != 0) || (rot[2] != 0))) ||
                    (loc && ((loc[0] != 0) || (loc[1] != 0) || (loc[2] != 0)));
  // fprintf(stderr, "[vol] transformation matrix:\n");
  // fprintf(stderr, "  %g %g %g  %g\n", tree->w2o[0][0], tree->w2o[0][1], tree->w2o[0][2], tree->w2o[0][3]);
  // fprintf(stderr, "  %g %g %g  %g\n", tree->w2o[1][0], tree->w2o[1][1], tree->w2o[1][2], tree->w2o[1][3]);
  // fprintf(stderr, "  %g %g %g  %g\n", tree->w2o[2][0], tree->w2o[2][1], tree->w2o[2][2], tree->w2o[2][3]);
}

// open the tree from disk again
static inline vol_tree_t *vol_open(const char *filename)
{
  vol_tree_t *tree = (vol_tree_t *)malloc(sizeof(*tree));
  tree->fd = open(filename, O_RDONLY);
  if(tree->fd < 0) return 0;
  tree->data_size = lseek(tree->fd, 0, SEEK_END);
  lseek(tree->fd, 0, SEEK_SET);
  common_readahead(tree->fd, 0, tree->data_size);
  tree->data = (uint8_t *)mmap(0, tree->data_size, PROT_READ, MAP_SHARED, tree->fd, 0);
  tree->header = (vol_header_t *)tree->data;
  if((tree->header->magic != VOL_HEADER_MAGIC) || (tree->header->version != (VOL_HEADER_VERSION | (VOL_MOTION_SAMPLES << 16))))
  {
    fprintf(stderr, "[vol_open] failed to open vol file! (magic: %X our %X, version: %d our %d, motion samples %d our %d)\n",
        tree->header->magic, VOL_HEADER_MAGIC, tree->header->version&0xffff, VOL_HEADER_VERSION, tree->header->version>>16, VOL_MOTION_SAMPLES);
    vol_close(tree);
    return 0;
  }
  memcpy(tree->aabb, tree->header->aabb, sizeof(float)*6);
  memcpy(tree->content_box, tree->header->content_box, sizeof(float)*6);
  tree->voxel_size = tree->header->voxel_size;
  vol_create_transform(tree, tree->header->loc, tree->header->rot);

  tree->nodes = (vol_node_t *)(tree->data + tree->header->nodes);
  tree->payload = tree->data + sizeof(vol_header_t);
  tree->light = (vol_payload_flux_t *)(tree->data + tree->header->light);
  tree->shader = vol_shader_get((vol_shader_id_t)tree->header->shaderid);

  // root node is the one we wrote out last:
  tree->root_node = (vol_node_t *)(tree->data + tree->header->light) - 1;
  size_t root_payload_size = tree->header->isstatic ? vol_payload_static_size() : sizeof(vol_payload_t);
  tree->root_payload = (vol_payload_t *)((uint8_t *)tree->nodes - root_payload_size);
  tree->root_light = (vol_payload_flux_t *)(tree->data + tree->header->end)-1;

  fprintf(stdout, "[vol_open] opened tree with depth %d and voxel size %f\n", tree->header->depth, tree->header->voxel_size);
  return tree;
}

static inline uint64_t vol_count_voxels(const vol_tree_t *tree)
{
  uint64_t num = 0;
  const vol_node_t *curr = tree->nodes;
  while(curr <= tree->root_node)
  {
    for(int i=0;i<512;i++)
      if(vol_node_get_payload(tree, curr, i))
        num += 512;
    curr++;
  }
  return num;
}

static inline void vol_transform_w2o(const vol_tree_t *tree, float *v, const int pos)
{
  if(tree->transform)
  {
    float tmp[3];
    for(int k=0;k<3;k++) tmp[k] = v[k] - (pos ? tree->w2o[k][3] : 0);
    v[0] = v[1] = v[2] = 0.0f;
    for(int j=0;j<3;j++)
      for(int i=0;i<3;i++)
        v[j] += tree->w2o[j][i] * tmp[i];
  }
}

static inline void vol_transform_o2w(const vol_tree_t *tree, float *v, const int pos)
{
  if(tree->transform)
  {
    float tmp[3] = {v[0], v[1], v[2]};
    v[0] = v[1] = v[2] = 0.0f;
    for(int j=0;j<3;j++)
      for(int i=0;i<3;i++)
        v[j] += tree->w2o[i][j] * tmp[i];
    if(pos) for(int k=0;k<3;k++) v[k] += tree->w2o[k][3];
  }
}

// point sample the voxel grid
static inline int vol_sample(
    const vol_tree_t *tree,  // voxel data
    int i, int j, int k,     // integer voxel coordinates, at finest scale
    int lod,                 // 0-finest, 1-one coarser, 2-.. tree->header->depth-1
    vol_channel_t channel,   // select channels to return, bitmask
    float time,              // time for motion blur
    float *result)           // store the result here
{
  vol_index_t idx[50];
  const int md = tree->header->depth;
  const int leaf_d = md > lod + 1 ? md-lod-1 : 0;
  for(int d=md-1;d>=0;d--)
  {
    idx[d].idx = 0;
    idx[d].i = i&7;
    idx[d].j = j&7;
    idx[d].k = k&7;
    i>>=3; j>>=3; k>>=3;
  }
  if(i || j || k) return 1; // index out of bounds
  const vol_node_t *n = tree->root_node;
  const vol_payload_t *data = tree->root_payload;
  int isstatic = tree->header->isstatic;
  for(int d=0;d<leaf_d;d++)
  {
    data = vol_node_get_payload(tree, n, idx[d].idx);
    assert(!data || (uint8_t*)data < (uint8_t*)tree->nodes);
    isstatic = vol_node_child_static(n, idx[d].idx);
    n = vol_node_get_child(tree, n, idx[d].idx);
    if(!n) break;
    assert(!n || (uint8_t*)n < tree->data + tree->data_size);
  }
  if(!data) return 1;

  vol_payload_get(data, idx[leaf_d].idx, time, isstatic, result);
  return 0;
}

