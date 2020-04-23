#include "corona_common.h"
#include "render.h"
#include "pathspace.h"
#include "points.h"
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
  fprintf(fd, "render   : geo visualisation\n");
}

void render_sample_path(uint64_t k)
{
  path_t path;

  int tid = common_get_threadid();

  double begin = common_time_wallclock();
  float col[3];
  path.index = k;
  uint64_t num_px = view_width()*view_height();
  uint64_t frame = k/num_px;
  uint64_t px = k - num_px*frame;
  const float y = px/view_width(), x = px - ((int)y)*view_width();
  path_init(&path, path.index, 0);
#if 0
  // test transmittance methods:
  path.lambda = 550.0f;
  path.time = 0.0f;
  path.v[0].hit.x[0] = -1 + (x/view_width())*2.0f;
  path.v[0].hit.x[1] = -2;
  path.v[0].hit.x[2] = (y/view_height())*2.0f;
  path.e[1].omega[0] = 0.0;
  path.e[1].omega[1] = 1.0f;
  path.e[1].omega[2] = 0.0;
  path.e[1].dist = 100.0;
  path.length = 1;
  path.sensor.pixel_i = x;
  path.sensor.pixel_j = y;
  path.sensor.camid = 0;
  shader_so_t *s = rt.shader->shader + 0;
  float transmittance = s->volume_transmittance(&path, 1, s->data);
  double end = common_time_wallclock();
  col[0] = transmittance;
  col[1] = transmittance;
  col[2] = 200000.0f*(end-begin);
  view_splat_col(&path, col);
  return;
#endif
  path_set_pixel(&path, x+points_rand(rt.points, tid),
                        y+points_rand(rt.points, tid)); // request pixel
  // path_set_aperture(&path, 0, 0);
  //     points_rand(rt.points, tid),
  //     points_rand(rt.points, tid)); // fix aperture
  if(path_extend(&path)) return;
  double end = common_time_wallclock();


#if 0
  col[2] = 0.0;
  col[0] = (1.0f+path.e[1].omega[0])/2.0f;
  col[1] = (1.0f+path.e[1].omega[1])/2.0f;
    view_splat_col(&path, col);
    return;
#endif

#if 1
  // time
  for(int k=0;k<3;k++)
    col[k] = 50000.0f * (end - begin);
  view_splat_col(&path, col);
  return;
#else
#if 0
  if(rt.lights->vol)
  {
    const float pdf = rt.lights->vol->volume_pdf_fnee_direction(&path, 1, rt.lights->vol->data);
    for(int k=0;k<3;k++) col[k] = pdf;
    view_splat_col(&path, col);
  }
  else
#endif

  if(!primid_invalid(path.v[1].hit.prim))
  {
    for(int k=0;k<3;k++)
      // normals
      col[k] = path.v[1].hit.n[k]/2.0f+.5f;
      // col[k] = hit->gn[k]/2.0f+.5f;
      // col[k] = -dotproduct(path.e[1].omega, hit->gn);
    view_splat_col(&path, col);
  }
  else if(path.e[1].dist < FLT_MAX)
  {
    for(int k=0;k<3;k++)
      col[k] = path.e[1].dist / 10;
    view_splat_col(&path, col);
  }
  else if(0)
  {
    const float lambda[3] = {400.0, 550.0, 650.0};
    for(int k=0;k<3;k++)
    {
      path.lambda = lambda[k];
      float power = .1*shader_sky_eval(&path, path.length-1);
      spectrum_p_to_camera(lambda[k], power, col);
      view_splat_col(&path, col);
    }
  }
#endif
}

void render_splat(const path_t *p, const mf_t value)
{
  view_splat(p, value);
}
