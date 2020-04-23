#include "corona_common.h"
#include "texture.h"
#include "prims.h"
#include "tools/geo/vdata.h"

typedef struct data_t
{
  int num_slots;
  tex_slot_t slots[8];
  float mul[8];
  uint32_t shapeid;

  vdata_t vdata;
}
data_t;

float prepare(path_t *p, int v, void *data)
{
  data_t *d = (data_t *)data;
  // TODO: get vertex indices for all corners of the primitive!
  // TODO: for every vertex and every channel, get data:
  const uint32_t vi = rt.prims->shape[p->v[v].hit.prim.shapeid].vtxidx[p->v[v].hit.prim.vi].v;
  // const uint32_t vcnt = p->v[v].hit.prim.vcnt;
  // for vcnt ..
  // TODO: more clever layout with quantised/compressed entries?
  // TODO: at least don't store 3 full channels for everything!
  const float *vd = vdata_get(&d->vdata, vi);
  assert(d->num_slots * (vi+1) * sizeof(float) <= d->vdata.data_size);
  assert(p->v[v].hit.prim.shapeid == d->shapeid);
  // for(int k=0;k<d->num_slots;k++)
  // TODO: respect stride ddd ggg or similar?
  // float tmp[3] = {0.0};
  // _geo_decode_uv(((uint32_t*)vd)[0], tmp, tmp+1);
  // tmp[0] = 4e3f*sqrtf(vd[1]*vd[1]+vd[2]*vd[2]);
  // tmp[1] = 4e3f*sqrtf(vd[3]*vd[3]+vd[4]*vd[4]);

  // yay, renderman style shaders ftw:
  const float scale = 5000.f;
  float tmp[3] = {-scale*vd[1], -scale*vd[2], 1.0};
  float ws[3] = {0.0f};
  normalise(tmp);
  for(int k=0;k<3;k++) ws[k] = p->v[v].hit.n[k] + p->v[v].hit.a[k] * tmp[0] + p->v[v].hit.b[k] * tmp[1];
  float dir[3] = {0, 0, 1.0};
  // quaternion_transform(&(rt.cam->orient), dir);
  tmp[0] = tmp[1] = tmp[2] = powf(-dotproduct(ws, dir), 3.0);
  // tmp[0] = (tmp[0] + 1.0)/2.0;
  // tmp[1] = (tmp[1] + 1.0)/2.0;
  // tmp[2] = 0.0;//vd[3];
#if 0
  // texture stretch:
  tmp[0] = pow(MAX(0.0, 1.0f -
    sqrtf(vd[3]*vd[3]+vd[4]*vd[4])/
    sqrtf(vd[1]*vd[1]+vd[2]*vd[2])), 2.);
  tmp[1] = 0.1f;
  tmp[2] = 0.0;//atan2f(vd[1], vd[2]);
#endif
    tex_set_slot_coeff(p, v, d->slots[1], d->mul[1], tmp);
  // TODO: interpolate with uv coords on primitive
  return 1.0;
}

int init(FILE *s, void **data)
{
  data_t *d = (data_t *)malloc(sizeof(data_t));
  memset(d, 0, sizeof(*d));
  *data = d;
  int err = 0;
  char filename[1024];
  err = fscanf(s, "%s", filename);
  // load channel configuration
  d->num_slots = 8;
  for(int k=0;k<8;k++)
  {
    d->mul[k] = 1.0f;
    char c;
    err = fscanf(s, " %c", &c);
    d->slots[k] = tex_parse_slot(c);
    if(err != 1 || d->slots[k] == s_slot_count)
    {
      d->num_slots = k;
      break;
    }
  }
  // map vertex data
  if(vdata_map(&d->vdata, filename))
  {
    char fname[1024];
    snprintf(fname, 1024, "%s/%s", rt.searchpath, filename);
    vdata_map(&d->vdata, fname);
  }
  if(d->vdata.fd == -1)
  {
    fprintf(stderr, "[vdata] could not load vertex data `%s'!\n", filename);
    return 1;
  }
  d->vdata.num_slots = d->num_slots;
  err = fscanf(s, "%*[^\n]\n"); // munch potential comments
  return err == -1;
}

void cleanup(void *data)
{
  data_t *d = (data_t*)data;
  vdata_cleanup(&d->vdata);
  free(d);
}

int shape_init(uint32_t shapeid, shader_so_t *self)
{
  // TODO: sanity check vertex count shapeid vs data size and slots!
  data_t *d = (data_t *)self->data;
  d->shapeid = shapeid;
  // TODO: vdata.num_verts = ?
  return 0;
}
