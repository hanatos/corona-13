#pragma once

#include "prims.h"
#include "geo.h"
#include "vdata.h"
#include <stdint.h>
#include <float.h>

// #define Float double
#define Float float

// need some more dynamic structure to put mesh into for subdivision:
typedef struct sd_vert_t
{
  uint32_t f;  // face index into faceidx and edgeidx arrays
  uint32_t e;  // edge index into faceidx and edgeidx arrays
  uint32_t nf; // number of adjacent faces
  uint32_t ne; // number of adjacent edges, different to nf for holes
}
sd_vert_t;

typedef struct sd_edge_t
{
  uint32_t v0, v1;   // pointing into verts array
  uint32_t f0, f1;   // pointing into faces array
}
sd_edge_t;

typedef struct sd_face_t
{
  uint32_t vi;    // index of first vertex index in vtxidx array
  uint32_t vcnt;  // number of vertices (would be 4 starting at 2nd iteration)
}
sd_face_t;

typedef struct sd_mesh_t
{
  int num_c;  // number of vertex coordinates per vertex (3 + 3 shutter close and uvs?)
  Float *vc;  // vertex coordinates: world space, motion blur pos, uvs.

  uint32_t num_verts, num_verts_allocated;
  uint32_t num_edges, num_edges_allocated;
  uint32_t num_faces, num_faces_allocated;
  sd_vert_t *verts;   // vertex index structs, parallel to vc above
  sd_edge_t *edges;   // edge list
  sd_face_t *faces;   // face list

  uint32_t num_vertidx, num_vertidx_allocated;
  uint32_t num_edgeidx, num_edgeidx_allocated;
  uint32_t num_faceidx, num_faceidx_allocated;
  uint32_t *vertidx;   // vertex indices, sorted by face
  uint32_t *edgeidx;   // edge indices, sorted by vertex
  uint32_t *faceidx;   // face indices, sorted by vertex

  uint32_t *writeout_order; // permutation to increase locality when writing to disk
  uint32_t *writeout_order_rev; // reverse mapping
}
sd_mesh_t;

static inline void sd_mesh_init(
    sd_mesh_t *m,
    uint32_t nv,
    uint32_t ne,
    uint32_t nf,
    uint32_t nvi,
    uint32_t num_c)
{
  memset(m, 0, sizeof(sd_mesh_t));

  m->vc = (Float *)malloc(sizeof(Float)*num_c*nv);
  m->verts = (sd_vert_t *)malloc(sizeof(sd_vert_t)*nv);
  m->edges = (sd_edge_t *)malloc(sizeof(sd_edge_t)*ne);
  m->faces = (sd_face_t *)malloc(sizeof(sd_face_t)*nf);
  memset(m->verts,  0, sizeof(sd_vert_t)*nv);
  memset(m->edges, -1, sizeof(sd_edge_t)*ne);
  memset(m->faces,  0, sizeof(sd_face_t)*nf);
  m->num_verts_allocated = nv;
  m->num_edges_allocated = ne;
  m->num_faces_allocated = nf;

  m->vertidx = (uint32_t *)malloc(sizeof(uint32_t)*nvi);
  memset(m->vertidx, -1, sizeof(uint32_t)*nvi);
  m->num_vertidx_allocated = nvi;

  m->num_c = num_c;
}

static inline void sd_mesh_cleanup(
    sd_mesh_t *m)
{
  free(m->vc);
  free(m->verts);
  free(m->edges);
  free(m->faces);
  free(m->vertidx);
  free(m->edgeidx);
  free(m->faceidx);
  free(m->writeout_order);
  free(m->writeout_order_rev);
  memset(m, 0, sizeof(sd_mesh_t));
}

static inline uint64_t _morton_pad(uint64_t x)
{
#if 0
  // we may have to use 3x 10 bits if compare only supports integer return?
  x = (x | (x << 16)) & 0x030000FF;
  x = (x | (x <<  8)) & 0x0300F00F;
  x = (x | (x <<  4)) & 0x030C30C3;
  x = (x | (x <<  2)) & 0x09249249;
#endif
  x &= 0x1fffff;
  x = (x | x << 32) & 0x1f00000000ffff;
  x = (x | x << 16) & 0x1f0000ff0000ff;
  x = (x | x << 8) & 0x100f00f00f00f00f;
  x = (x | x << 4) & 0x10c30c30c30c30c3;
  x = (x | x << 2) & 0x1249249249249249;
  return x;
}

// compute 64-bit morton code from 3 21-bit integers:
static inline uint64_t _morton(uint32_t xx, uint32_t yy, uint32_t zz)
{
  const uint64_t x = _morton_pad(xx);
  const uint64_t y = _morton_pad(yy);
  const uint64_t z = _morton_pad(zz);
  return x | (y<<1) | (z<<2);
}

static inline int _morton_compare(const void *aa, const void *bb)
{
  const uint64_t *a = (const uint64_t *)aa, *b = (const uint64_t *)bb;
  const int64_t diff = a[1] - b[1];
  // can only return 32-bit int:
  return diff > 0 ? 1 : (diff < 0 ? -1 : 0);
}

// create morton array to sort vertices by spatial locality
static inline void sd_mesh_create_writeout_order(
    sd_mesh_t *mesh)
{
  // find (shutter open) bounding box
  float aabb[6] = {FLT_MAX,FLT_MAX,FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX};
  for(int v=0;v<mesh->num_verts;v++) for(int k=0;k<3;k++)
  {
    aabb[k  ] = MIN(aabb[k], mesh->vc[mesh->num_c*v + k]);
    aabb[3+k] = MAX(aabb[k], mesh->vc[mesh->num_c*v + k]);
  }
  // create array with morton indices of vertex:
  // (could pass mesh to sort instead and only alloc 1xu32)
  uint64_t *sorted = (uint64_t *)malloc(2*mesh->num_verts*sizeof(uint64_t));
  for(int v=0;v<mesh->num_verts;v++)
  {
    sorted[2*v] = v;
    sorted[2*v+1] = _morton(
        (mesh->vc[mesh->num_c*v+0] - aabb[0])*0x1fffff/(aabb[3]-aabb[0]),
        (mesh->vc[mesh->num_c*v+1] - aabb[1])*0x1fffff/(aabb[4]-aabb[1]),
        (mesh->vc[mesh->num_c*v+2] - aabb[2])*0x1fffff/(aabb[5]-aabb[2]));
  }
  qsort(sorted, mesh->num_verts, 2*sizeof(uint64_t), &_morton_compare);
  // this array shall encode where on disk you will actually find
  // the vertex if you look for vertex v: writeout_order[v].
  mesh->writeout_order = (uint32_t *)malloc(mesh->num_verts*sizeof(uint32_t));
  mesh->writeout_order_rev = (uint32_t *)malloc(mesh->num_verts*sizeof(uint32_t));
  for(int v=0;v<mesh->num_verts;v++)
  {
    mesh->writeout_order[sorted[2*v]] = v;
    mesh->writeout_order_rev[v] = sorted[2*v];
  }
  free(sorted);
}

static inline int sd_mesh_write_geo(
    const sd_mesh_t *mesh,
    const char *geoname)
{
  FILE *o = fopen(geoname, "wb");
  if(!o) return 1;

  prims_header_t header;
  header.magic = GEO_MAGIC;
  header.version = GEO_VERSION;
  header.num_prims = mesh->num_faces;
  uint64_t cnt = sizeof(prims_header_t) + sizeof(primid_t)*mesh->num_faces;
  header.vtxidx_offset = cnt;
  cnt += sizeof(prims_vtxidx_t)*mesh->num_vertidx;
  header.vertex_offset = (cnt + 0xf) & ~0xf; // round up to next multiple of 16 bytes for sse alignment
  assert(header.vertex_offset >= cnt);
  assert((header.vertex_offset & 0xf) == 0);
  fwrite(&header, sizeof(prims_header_t), 1, o);
  for(int f=0;f<mesh->num_faces;f++)
  {
    primid_t p = {0};
    p.vcnt = mesh->faces[f].vcnt;
    p.vi = mesh->faces[f].vi;
    fwrite(&p, sizeof(primid_t), 1, o);
  }
  for(int i=0;i<mesh->num_vertidx;i++)
  {
    prims_vtxidx_t vi = {0};
    vi.v = mesh->vertidx[i];
    vi.uv = geo_encode_uv(mesh->vc[mesh->num_c*vi.v + mesh->num_c-2], mesh->vc[mesh->num_c*vi.v + mesh->num_c-1]);
    if(mesh->writeout_order) vi.v = mesh->writeout_order[mesh->vertidx[i]];
    fwrite(&vi, sizeof(prims_vtxidx_t), 1, o);
  }
  if(header.vertex_offset > cnt) // pad to sse alignment
    fwrite(mesh->verts, header.vertex_offset - cnt, 1, o);
  for(int vv=0;vv<mesh->num_verts;vv++)
  {
    int v = vv;
    if(mesh->writeout_order) v = mesh->writeout_order_rev[vv];
    prims_vtx_t vtx = {0};
    for(int i=0;i<3;i++) vtx.v[i] = mesh->vc[mesh->num_c*v+i];
    fwrite(&vtx, sizeof(vtx), 1, o);
    if(mesh->num_c == 8)
    { // motion blur
      for(int i=0;i<3;i++) vtx.v[i] = mesh->vc[mesh->num_c*v+3+i];
      fwrite(&vtx, sizeof(vtx), 1, o);
    }
  }
  fclose(o);

  // create vertex normals on geo by loading it back in
  prims_t prims;
  prims_init(&prims);
  prims_allocate(&prims, 1);
  char basename[1024];
  strncpy(basename, geoname, 1024);
  basename[strlen(geoname)-4] = '\0';
  if(prims_load_with_flags(&prims, basename, "none", 0, 'w', 0))
  {
    fprintf(stderr, "[sd_mesh_write_geo] failed to read back `%s' to compute normals!\n", geoname);
    return 2;
  }
  geo_recompute_normals(prims.shape);
  prims_cleanup(&prims);
  return 0;
}

static inline int sd_mesh_add_vert(
    sd_mesh_t *mesh,
    float *x0,
    float *x1,
    float *uv)
{
  int v = mesh->num_verts++;
  if(uv) for(int i=0;i<2;i++) mesh->vc[mesh->num_c*v+mesh->num_c-2+i] = uv[i];
  for(int i=0;i<3;i++) mesh->vc[mesh->num_c*v+i] = x0[i];
  if(mesh->num_c == 8) for(int i=0;i<3;i++) mesh->vc[mesh->num_c*v+3+i] = x1[i];
  mesh->verts[v].f = -1;
  mesh->verts[v].e = -1;
  mesh->verts[v].nf = 0;
  mesh->verts[v].ne = 0;
  assert(mesh->num_verts <= mesh->num_verts_allocated);
  return v;
}

static inline void sd_mesh_add_face(
    sd_mesh_t *mesh,
    uint32_t vcnt,     // only support insertion of tris or quads for subd
    uint32_t v0,
    uint32_t v1,       // vi point into vertex struct,
    uint32_t v2,       // this function will create the vertidx entry.
    uint32_t v3)
{
  int f = mesh->num_faces++;
  int vi = mesh->num_vertidx;
  mesh->num_vertidx += vcnt;

  mesh->faces[f].vi = vi;
  mesh->faces[f].vcnt = vcnt;
  mesh->vertidx[vi+0] = v0;
  mesh->vertidx[vi+1] = v1;
  mesh->vertidx[vi+2] = v2;
  if(vcnt > 3) mesh->vertidx[vi+3] = v3;

  assert(v0 < mesh->num_verts);
  assert(v1 < mesh->num_verts);
  if(vcnt > 2) assert(v2 < mesh->num_verts);
  if(vcnt > 3) assert(v3 < mesh->num_verts);

  // for now only count number of adjacent faces, we'll come back to this after
  // all faces have been created, namely in sd_mesh_update_indices().
  mesh->verts[v0].nf++;
  mesh->verts[v1].nf++;
  mesh->verts[v2].nf++;
  if(vcnt > 3) mesh->verts[v3].nf++;

  assert(mesh->num_faces <= mesh->num_faces_allocated);
  assert(mesh->num_vertidx <= mesh->num_vertidx_allocated);
}

// given two vertex indices, search the first vertex for a matching edge reference
static inline int sd_mesh_find_edge(
    const sd_mesh_t *m,
    const int vi0,
    const int vi1)
{
  const int v0 = m->vertidx[vi0];
  const int v1 = m->vertidx[vi1];
  const int e  = m->verts[v0].e;    // offset in edgeidx array for vertex 0
  const int n  = m->verts[v0].nf+1; // number of entries in edgeidx array is bounded by this (for holes)
  int ei = e;
  for(;ei<e+n&&m->edgeidx[ei]!=-1;ei++)
  { // find all edges of current vertex and compare
    sd_edge_t *edge = m->edges + m->edgeidx[ei];
    if((edge->v0 == v0 && edge->v1 == v1) ||
        (edge->v0 == v1 && edge->v1 == v0))
      return ei;
  }
  // this may happen for two reasons: this edge not found yet (v0/v1 differ)
  // or this vertex is at a boundary, where you actually have one edge more than faces.
  assert(m->edgeidx[ei] == -1);
  return ei;
}

static inline void sd_mesh_update_indices(
    sd_mesh_t *mesh)
{
  // go through all vert_t structs and init pointers to edge and face index arrays
  uint32_t fcnt = 0, ecnt = 0;
  for(int v=0;v<mesh->num_verts;v++)
  {
    mesh->verts[v].f = fcnt;
    mesh->verts[v].e = ecnt;
    mesh->verts[v].ne = 0; // init later, after deduping
    fcnt += mesh->verts[v].nf;
    ecnt += mesh->verts[v].nf + 1; // security margin, we may be a border vertex with for instance only one face and two edges
  }
  if(mesh->edgeidx == 0 || mesh->num_edgeidx != ecnt)
  {
    free(mesh->edgeidx);
    mesh->edgeidx = (uint32_t *)malloc(sizeof(uint32_t)*ecnt);
  }
  if(mesh->faceidx == 0 || mesh->num_faceidx != fcnt)
  {
    free(mesh->faceidx);
    mesh->faceidx = (uint32_t *)malloc(sizeof(uint32_t)*fcnt);
  }
  memset(mesh->edgeidx, -1, sizeof(uint32_t)*ecnt);
  memset(mesh->faceidx, -1, sizeof(uint32_t)*fcnt);
  mesh->num_edgeidx_allocated = mesh->num_edgeidx = ecnt;
  mesh->num_faceidx_allocated = mesh->num_faceidx = fcnt;
  // go through all faces, access all associated vert_t structs and write face index to faceidx array
  for(int f=0;f<mesh->num_faces;f++)
  {
    for(int v=0;v<mesh->faces[f].vcnt;v++)
    {
      const sd_vert_t *vx = mesh->verts + mesh->vertidx[mesh->faces[f].vi + v];
      for(int ff=0;ff<vx->nf;ff++)
      {
        if(mesh->faceidx[vx->f+ff] == -1)
        {
          mesh->faceidx[vx->f+ff] = f;
          break;
        }
        else if(ff == vx->nf-1) fprintf(stderr, "[subd mesh update]: no more free face index slot!\n");
      }
    }
  }

  mesh->num_edges = 0;
  for(int f=0;f<mesh->num_faces;f++)
  { // for all faces, try to insert all edges
    const int vcnt = mesh->faces[f].vcnt;
    for(int k=0;k<vcnt;k++)
    { // go through all vertices and adjacent edge
      const int vi0 = mesh->faces[f].vi+k;
      const int vi1 = mesh->faces[f].vi+((k+1)%vcnt);
      int ei0 = sd_mesh_find_edge(mesh, vi0, vi1);
      int ei1 = sd_mesh_find_edge(mesh, vi1, vi0);
      if(mesh->edgeidx[ei0] != -1)
      { // found one at v0, update entry
        sd_edge_t *edge = mesh->edges + mesh->edgeidx[ei0];
        assert(edge->f0 != -1);
        assert(edge->f1 == -1);
        edge->f1 = f;
        if(mesh->edgeidx[ei1] == -1)
        {
          mesh->verts[mesh->vertidx[vi1]].ne++;
          mesh->edgeidx[ei1] = mesh->edgeidx[ei0];
        }
        else assert(mesh->edgeidx[ei1] == mesh->edgeidx[ei0]);
      }
      else if(mesh->edgeidx[ei1] != -1)
      { // found one at v1, update entry
        sd_edge_t *edge = mesh->edges + mesh->edgeidx[ei1];
        assert(edge->f0 != -1);
        assert(edge->f1 == -1);
        edge->f1 = f;
        if(mesh->edgeidx[ei0] == -1)
        {
          mesh->verts[mesh->vertidx[vi0]].ne++;
          mesh->edgeidx[ei0] = mesh->edgeidx[ei1];
        }
        else assert(mesh->edgeidx[ei1] == mesh->edgeidx[ei0]);
      }
      else
      { // allocate new slot
        mesh->verts[mesh->vertidx[vi0]].ne++;
        mesh->verts[mesh->vertidx[vi1]].ne++;
        mesh->edgeidx[ei0] = mesh->edgeidx[ei1] = mesh->num_edges;
        sd_edge_t *edge = mesh->edges + mesh->num_edges++;
        assert(mesh->num_edges <= mesh->num_edges_allocated);
        edge->v0 = mesh->vertidx[vi0];
        edge->v1 = mesh->vertidx[vi1];
        edge->f0 = f;
        edge->f1 = -1;
      }
      assert(mesh->edgeidx[ei0] != -1);
      assert(mesh->edgeidx[ei1] != -1);
      assert(mesh->edgeidx[ei0] == mesh->edgeidx[ei1]);
    }
  }
}

// write extra vertex data along with shape.geo file that will be lost
// otherwise because it's not needed for rendering.  in particular, this
// currently writes per-vertex uv-space ellipses for anisotropic texture
// lookups.
static inline int sd_mesh_write_vdata(
    const sd_mesh_t *m,
    const char *filename)
{
  const int c = m->num_c;         // number of coordinates per vertex
  const int uvo = c == 8 ? 6 : 3; // offset to uv coordinates
  vdata_t vdata;
  // XXX TODO: check whether a covariance matrix (3 floats) would be enough!
  vdata_alloc(&vdata, m->num_verts, 5);
  for(int vv=0;vv<m->num_verts;vv++)
  {
    int v = vv;
    if(m->writeout_order) v = m->writeout_order_rev[vv];
    // go through all edges of this vert, collect uv deltas
    // then get approximate ewa axes:
    // get max edge => da, get max projection to perpendicular => db
    float da[2] = {0.0f}, db[2] = {0.0f};
    const int ne = m->verts[v].ne;
    float duv[ne][2];
    for(int e=0;e<ne;e++)
    {
      const int v0 = m->edges[m->edgeidx[m->verts[v].e + e]].v0;
      const int v1 = m->edges[m->edgeidx[m->verts[v].e + e]].v1;
      const int ve = v0 == v ? v1 : v0; // other vertex on edge
      for(int i=0;i<2;i++) duv[e][i] = m->vc[c*ve + uvo + i] - m->vc[c*v + uvo + i];
      if(duv[e][0]*duv[e][0] + duv[e][1]*duv[e][1] > da[0]*da[0] + da[1]*da[1])
      {
        da[0] = duv[e][0];
        da[1] = duv[e][1];
      }
    }
    float pa[2] = {-da[1], da[0]};
    const float ilen_da = 1.0f/sqrtf(da[0]*da[0] + da[1]*da[1]);
    for(int e=0;e<ne;e++)
    {
      if(fabsf(duv[e][0]*pa[0] + duv[e][1]*pa[1]) > fabsf(db[0]*pa[0] + db[1]*pa[1]))
      {
        // orthogonalise
        const float dota = (da[0]*duv[e][0] + da[1]*duv[e][1])*ilen_da;
        db[0] = duv[e][0] - dota*da[0];
        db[1] = duv[e][1] - dota*da[1];
      }
    }
    // TODO: at least store all 6 in half?
    float *vd = vdata_get(&vdata, vv);
    ((uint32_t *)vd)[0] = geo_encode_uv(m->vc[c*v + uvo], m->vc[c*v + uvo+1]);
    vd[1] = da[0];
    vd[2] = da[1];
    vd[3] = db[0];
    vd[4] = db[1];
  }
  int res = vdata_write(&vdata, filename);
  vdata_cleanup(&vdata);
  return res;
}

static inline void sd_mesh_subdiv(
    const sd_mesh_t *in,
    sd_mesh_t *out)
{
  // will chop every face into vcnt new faces
  // will introduce one new vertex per face and per edge
  int num_out_faces = 0; // == 4*in->num_faces after first iteration
  for(int f=0;f<in->num_faces;f++)
    num_out_faces += in->faces[f].vcnt;
  // allocate storage for output:
  sd_mesh_init(out,
      in->num_verts + in->num_edges + in->num_faces,
      2*in->num_edges + num_out_faces,
      num_out_faces,
      num_out_faces*4, // will all be quads
      in->num_c);

  // face vertices at center of input faces:
  const int fvoff = in->num_verts;
#pragma omp parallel for schedule(static) default(shared)
  for(int f=0;f<in->num_faces;f++)
  {
    // average all vertices
    for(int i=0;i<out->num_c;i++) out->vc[out->num_c*(fvoff+f)+i] = 0.0;
    for(int k=0;k<in->faces[f].vcnt;k++)
      for(int i=0;i<out->num_c;i++)
        out->vc[out->num_c*(fvoff+f)+i] += in->vc[in->num_c*in->vertidx[in->faces[f].vi+k]+i]/in->faces[f].vcnt;
  }
  // edge vertices at center of two adjacent verts and two faces:
  const int evoff = in->num_verts + in->num_faces;
#pragma omp parallel for schedule(static) default(shared)
  for(int e=0;e<in->num_edges;e++)
  {
    if(in->edges[e].f1 == -1)
      for(int i=0;i<out->num_c;i++)
        out->vc[out->num_c*(evoff+e)+i] = (in->vc[in->num_c*in->edges[e].v0+i] + in->vc[in->num_c*in->edges[e].v1+i])/2.0;
    else
      for(int i=0;i<out->num_c;i++)
        out->vc[out->num_c*(evoff+e)+i] = (in->vc[in->num_c*in->edges[e].v0+i] + in->vc[in->num_c*in->edges[e].v1+i] +
            out->vc[out->num_c*(fvoff+in->edges[e].f0)+i] + out->vc[out->num_c*(fvoff+in->edges[e].f1)+i])/4.0;
  }
  out->num_verts = in->num_verts + in->num_faces + in->num_edges;

  // insert faces
  // could run in parallel with some added index magic in add_face()
  // (need to know output face number and vertidx offset)
  // #pragma omp parallel for schedule(static) default(shared)
  for(int f=0;f<in->num_faces;f++)
  { //  emit vcnt new faces and edges
    const int vcnt = in->faces[f].vcnt;
    uint32_t ev[vcnt];
    for(int k=0;k<vcnt;k++)
      ev[k] = sd_mesh_find_edge(in, in->faces[f].vi+k, in->faces[f].vi+((k+1)%vcnt));
    for(int k=0;k<vcnt;k++)
      sd_mesh_add_face(out, 4,
          in->vertidx[in->faces[f].vi + k],
          evoff + in->edgeidx[ev[k]],
          fvoff + f,
          evoff + in->edgeidx[ev[(k+vcnt-1)%vcnt]]);
  }

  // update references and edges
  sd_mesh_update_indices(out);

  // move original points: first init to 0
  memset(out->vc, 0, sizeof(float)*out->num_c*in->num_verts);
#pragma omp parallel for schedule(static) default(shared)
  for(int v=0;v<in->num_verts;v++)
  {
    const int nf = in->verts[v].nf;
    const int ne = in->verts[v].ne;
    if(ne <= 2) // special case for valence 2 (corners)
      for(int i=0;i<out->num_c;i++)
        out->vc[out->num_c*v+i] = in->vc[in->num_c*v+i];
    else if(ne > nf) // special case for holes (edges)
    {
      int cnt = 0;
      for(int k=0;k<ne;k++) // loop over all adjacent edges and sum up edge midpoints only for boundary edges
        if(in->edges[in->edgeidx[in->verts[v].e + k]].f1 == -1)
        {
          cnt++;
          for(int i=0;i<out->num_c;i++)
          {
            out->vc[out->num_c*v+i] += (in->vc[in->num_c*in->edges[in->edgeidx[in->verts[v].e + k]].v0+i] +
                in->vc[in->num_c*in->edges[in->edgeidx[in->verts[v].e + k]].v1+i]);
          }
        }
      if(cnt) for(int i=0;i<out->num_c;i++) // normalise
        out->vc[out->num_c*v+i] /= 2.0*cnt;
      else for(int i=0;i<out->num_c;i++)
        out->vc[out->num_c*v+i] = in->vc[in->num_c*v+i];
    }
    else // regular case
    {
      for(int k=0;k<nf;k++) // add average face vertex
        for(int i=0;i<out->num_c;i++)
          out->vc[out->num_c*v+i] += out->vc[out->num_c*(fvoff + in->faceidx[in->verts[v].f + k])+i] / nf;
      for(int k=0;k<ne;k++) // loop over all adjacent edges and sum up edge midpoints twice
        for(int i=0;i<out->num_c;i++)
          out->vc[out->num_c*v+i] += (in->vc[in->num_c*in->edges[in->edgeidx[in->verts[v].e + k]].v0+i] +
              in->vc[in->num_c*in->edges[in->edgeidx[in->verts[v].e + k]].v1+i]) / ne;
      for(int i=0;i<out->num_c;i++) // sum up old vertex part
        out->vc[out->num_c*v+i] += (ne-3.0) * in->vc[in->num_c*v+i];
      for(int i=0;i<out->num_c;i++) // normalise edge and face vertex part
        out->vc[out->num_c*v+i] /= ne;
    }
  }
}

static inline int sd_mesh_init_from_geo(
    sd_mesh_t *in,
    const char *filename)
{
  prims_t prims;
  prims_init(&prims);
  prims_allocate(&prims, 1);
  prims_load_with_flags(&prims, filename, "none", 0, 'r', 0);
  if(prims.num_loaded_shapes == 0)
  {
    prims_cleanup(&prims);
    return 1;
  }

  prims_shape_t *s = prims.shape;

  const int mb = s->primid[0].mb+1;
  uint32_t num_orig_verts = (s->data_size - ((uint8_t *)s->vtx - (uint8_t *)s->data))/(sizeof(prims_vtx_t) * mb);
  uint32_t num_orig_faces = s->num_prims;
  uint32_t num_orig_edges = 0;
  for(int f=0;f<num_orig_faces;f++)
    if(s->primid[f].vcnt > 2 && s->primid[f].vcnt < 5)
      num_orig_edges += s->primid[f].vcnt;
  uint32_t num_orig_vertidx = 0; // stupid:
  for(int f=0;f<num_orig_faces;f++)
    if(s->primid[f].vi + s->primid[f].vcnt > num_orig_vertidx)
      num_orig_vertidx = s->primid[f].vi + s->primid[f].vcnt;

  // allocate storage arrays
  // number of edges is a blatant lie (upper bound)
  sd_mesh_init(in, num_orig_verts, num_orig_edges, num_orig_faces, num_orig_vertidx, 3*mb + 2);
  // fill in position buffers and uvs of all vertices:
  for(int v=0;v<num_orig_verts;v+=mb)
    sd_mesh_add_vert(in, s->vtx[v].v, s->vtx[v+1].v, 0);
  for(int vi=0;vi<num_orig_vertidx;vi++)
  {
    float uv[2];
    geo_decode_uv(s->vtxidx[vi].uv, uv, uv+1);
    in->vc[in->num_c*s->vtxidx[vi].v + in->num_c-2] = uv[0];
    in->vc[in->num_c*s->vtxidx[vi].v + in->num_c-1] = uv[1];
  }

  // for all faces: fill in vtxidx and face array
  // needs an atomic++ (variable vcnt)
  for(int f=0;f<num_orig_faces;f++)
    if(s->primid[f].vcnt > 2 && s->primid[f].vcnt < 5)
      sd_mesh_add_face(in, s->primid[f].vcnt,
          s->vtxidx[s->primid[f].vi+0].v,
          s->vtxidx[s->primid[f].vi+1].v,
          s->vtxidx[s->primid[f].vi+2].v,
          s->vtxidx[s->primid[f].vi+3].v);

  // done with input
  prims_cleanup(&prims);

  // init faceidx array to point from vertices to faces:
  sd_mesh_update_indices(in);
  return 0;
}

