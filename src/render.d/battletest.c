#include "corona_common.h"
#include "render.h"
#include "pathspace.h"
#include "view.h"
#include "shader.h"
#include "spectrum.h"
#include "lights.h"

typedef struct render_t
{
  int keep;
}
render_t;

typedef struct render_tls_t
{
  int keep;
}
render_tls_t;

render_t *render_init()
{
  // TODO: grab some cmdline options:
  // - get scene info or shader string(roughness, ior, ..?) ? or always use last shader in list?
  // - mode is bsdf, pdf, or sample (or do all in one pass per channel?)
  // - testing inside out or outside in
  // - test adjoint pdf?
  return 0;
}

void render_cleanup(render_t *r) { }

render_tls_t *render_tls_init()
{
  return 0;
}

void render_tls_cleanup(render_tls_t *r) { }

void render_clear() { }


void render_print_info(FILE *fd)
{
  fprintf(fd, "render   : bsdf battle testing\n");
}

void render_sample_path(uint64_t k)
{
  path_t path;

  // TODO: use anim_frame for wi inclination angle?
  uint64_t num_px = view_width()*view_height();
  uint64_t frame = k/num_px;
  uint64_t px = k - num_px*frame;
  const float y = px/view_width(), x = px - y*view_width();

  float col[3] = {1,1,1};
  path.index = k;
  path_init(&path, path.index, 0);
  path_set_pixel(&path, x, y); // request pixel
  path_set_aperture(&path, 0, 0); // fix aperture point (pinhole)
  // if(path_extend(&path)) return;
  path.lambda = 550.0f;
  path.time = 0.0f;
  path.sensor.camid = 0;

  path.length = 2; // now sampling v[2]
  hit_t hit = {{0}};
  hit.n[2] = 1.0f;
  hit.gn[2] = 1.0f;
  hit.shader = 10; // XXX get from cmdline or something!
  hit.prim.extra = 0;
  hit.prim.shapeid = 0;
  hit.prim.vi = 0;
  hit.prim.mb = 0;
  hit.prim.vcnt = 4;
  // if(!REFLECT) hit->n[2] = hit->gn[2] = -1.0f;
  get_onb(hit.n, hit.a, hit.b);

  path.v[0].hit.prim = INVALID_PRIMID;
  path_volume_vacuum(&path.e[0].vol);
  path.v[0].interior = path.e[0].vol;
  path.e[1].vol = path.e[0].vol;

  vertex_shading_t *sh = &path.v[1].shading;
  sh->rs = 0.06f;
  sh->rd = 1.0f;
  sh->em = 0.0f;
  sh->rg = 1.0f;
  sh->roughness = 0.6; // XXX get from cmdline? will be overwritten by shader_prepare..

  // incoming direction
  path.e[1].omega[0] = 0.5;
  path.e[1].omega[1] = 0.0;
  const int reflect = 0;
  if(reflect)
  {
    path.e[1].omega[2] = -1.0;
    hit.n[2] = 1.0; // would be flipped to make omega point towards it
    hit.gn[2] = 1.0; // we just always hit it from the right side
  }
  else
  {
    path.e[1].omega[2] = 1.0;
    hit.n[2] = hit.gn[2] = -1.0;
  }
  normalise(path.e[1].omega);

  // shader_prepare(&path, 1);
  shader_so_t *s = rt.shader->shader + path.v[1].hit.shader;
  if(s->prepare) s->prepare(&path, 1, s->data);
  path.v[1].hit = hit; // overwrite tangent frame/normals
  path.e[0].vol.shader =
  path.e[1].vol.shader = -1;
  path.e[0].vol.ior =
  path.e[1].vol.ior = 1.0f; // XXX get from cmdline?
  path.v[1].interior.ior = 1.4;
  path.v[1].interior.shader = 0;
  if(reflect)
  { // reflection:
    path.e[2].vol.shader = path.e[1].vol.shader;
    path.e[2].vol.ior = path.e[1].vol.ior;
  }
  else
  { // transmission:
    path.e[2].vol.shader = path.v[1].interior.shader;
    path.e[2].vol.ior = path.v[1].interior.ior;
  }

  path.v[1].mode = s_absorb;

  // compute bsdf:
  path.e[2].omega[0] = 2.0* path.sensor.pixel_i/(float)view_width() - 1.;
  path.e[2].omega[1] = 2.0* path.sensor.pixel_j/(float)view_height() - 1.;
  const float len2 = path.e[2].omega[0]*path.e[2].omega[0]+path.e[2].omega[1]*path.e[2].omega[1];
  if(len2 < 1.0)
  {
    path.e[2].omega[2] = sqrtf(MAX(0.0, 1.0-len2));
    normalise(path.e[2].omega);
    // cos is taken care of by density dictated by i,j loop.
    const float bsdf = shader_brdf(&path, 1);
    assert(bsdf >= 0.0);
    col[0] = bsdf * 4.0f; // we're sampling the square, not the disc
    col[1] = 0.0f;
    col[2] = 0.0f;
    view_splat_col(&path, col);
  }

  // sample:
  col[0] = 0.0f;
  col[1] = shader_sample(&path);
  col[2] = 0.0f;
  if(col[1] > 0 && path.e[2].omega[2] > 0.0)
  {
    path.sensor.pixel_i = view_width()  * (path.e[2].omega[0] + 1.0)/2.0;
    path.sensor.pixel_j = view_height() * (path.e[2].omega[1] + 1.0)/2.0;
    view_splat_col(&path, col);
  }

  // TODO: separate into outgoing channels (bsdf, pdf, histogram of sampling)
  // TODO: pixel is outgoing direction for bsdf and pdf eval
  // TODO: pixel is random for 
}

void render_splat(const path_t *p, const float value)
{
  view_splat(p, value);
}
