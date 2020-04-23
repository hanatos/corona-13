#include "render.h"
#include "vmlt_hslt.h"

void render_accum(const path_t *p, const float value)
{
  // reject unwanted paths here
  float e[3];
  if(!isfinite(value)) return;

  // try to mutate using hslt:
  halfvec_stats_t stats[rt.num_threads];
  vmlt_hslt_t data; // XXX ouch. redo that but in better (depends on assumption that this is the only thing in this struct. currently true).
  data.stats = stats;

  int interesting = 1;
  // int interesting = 0;
  // for(int k=1;k<p->length-1;k++) if(p->v[k].mode & s_glossy) interesting = 1;

  if(hslt_suitability(p, &data) <= 0.0 || !interesting)
  {
    spectrum_p_to_camera(p->lambda, value * ERPT_MUTATIONS, e);
    filter_accum(p->sensor.pixel_i, p->sensor.pixel_j, e, rt.render->fb, 3, 0, 3);
    return;
  }
  const int mutations = ERPT_MUTATIONS;
  const int chains = MAX(1, (value/mutations)/0.2);
  const float weight = 1.0/chains;
  path_t path1, path2;
  path_t *curr = &path1, *tent = &path2;
  for(int c=0;c<chains;c++)
  {
    if(p->v[0].mode & s_sensor)
      *curr = *p; // copy, inefficient :(
    else
      path_reverse(curr, p);
    float accum = 0.0f;
    for(int m=0;m<mutations;m++)
    { // stolen metropolis acceptance logic from pointsampler_vmlt.c:
      float a = hslt_mutate(curr, tent, &data);
      if(!(a > 0.0)) a = 0;
      a = fminf(1.0f, a);
      const float w_tent = value * a;
      const float w_curr = value * (1.0f - a);
      if(curr->length > 0) accum += w_curr;
      if(points_rand(rt.points, common_get_threadid()) < a)
      { // accept
        if(accum > 0.0 && curr->length > 0)
        {
          spectrum_p_to_camera(curr->lambda, accum*weight, e);
          filter_accum(curr->sensor.pixel_i, curr->sensor.pixel_j, e, rt.render->fb, 3, 0, 3);
        }
        accum = w_tent;
        // swap paths:
        path_t *tmp = curr;
        curr = tent;
        tent = tmp;
      }
      else
      { // reject
        if(a > 0.0f)
        {
          spectrum_p_to_camera(tent->lambda, w_tent*weight, e);
          filter_accum(tent->sensor.pixel_i, tent->sensor.pixel_j, e, rt.render->fb, 3, 0, 3);
        }
      }
    }
    // accumulate the rest:
    if(accum > 0.0 && curr->length > 0)
    {
      spectrum_p_to_camera(curr->lambda, accum*weight, e);
      filter_accum(curr->sensor.pixel_i, curr->sensor.pixel_j, e, rt.render->fb, 3, 0, 3);
    }
  }
}

// render path corresponding to given sample index
void render_sample(uint64_t sample)
{
  rt_tls.render->tent_path->index = sample;
  pointsampler_mutate(rt_tls.render->curr_path, rt_tls.render->tent_path);
  if(pointsampler_accept(rt_tls.render->curr_path, rt_tls.render->tent_path))
  {
    path_t *tmp = rt_tls.render->curr_path;
    rt_tls.render->curr_path = rt_tls.render->tent_path;
    rt_tls.render->tent_path = tmp;
  }
  if(((sample+1) % (rt.width*rt.height/ERPT_MUTATIONS)) == 0) rt.render->overlays ++;
}

// update pixel index k
void render_update(uint64_t k)
{
  const float inv_o = rt.cam->iso/100.0f * 1.0f/rt.render->overlays;
  // convert camera/framebuffer to display profile:
  float xyz[3], rgb[3];
  for(int i=0;i<3;i++) rgb[i] = rt.render->fb[3*k+i]*inv_o;
  colour_camera_to_xyz(rgb, xyz);
  colour_xyz_to_output(xyz, rt.render->pixel + 4*k + 1);
}
