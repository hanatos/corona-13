#include "render.h"
#include "pathspace/halfvec.h"
#include "pathspace/multichain.h"

void render_accum(const path_t *p, const float value)
{
  // reject unwanted paths here
  float col[3];
  if(value <= 0 || !isfinite(value)) return;
  spectrum_p_to_camera(p->lambda, value, col);
  // fill visualisation into e here
  filter_accum(p->sensor.pixel_i, p->sensor.pixel_j, col, rt.render->fb, 3, 0, 3);

  // now the hard work:
  // duplicate path and sample another four of them through neighbouring pixels.
  halfvec_stats_t d = {0};
  path_t path = *p;

  path_t curr = (path_t *)p;
  path_t tent = &path;

  // we need to init half vectors to be able to perturb them!
  // init current path's half vectors, and constraint derivatives.
  curr->v[0].diffgeo.type = s_pinned_position;
  for(int i=1;i<curr->length;i++)
    curr->v[i].diffgeo.type = s_free;
  if(manifold_compute_tangents(curr, 0, curr->length-1) == 0.0) return;

  // first very wastefully copy the whole path
  // TODO: really not necessary, could only do vertex by vertex (multichain inits them anyhow)
  // TODO: maybe do this, but then construct multichain on different, empty path
  // and use this to converge the halfvector space part in the same index range?
  *tent = *curr;

  // determine vertex indices a=0, b=?, c=?:
  int b = tent->length - 1; // default to all lens perturbation
  float pdf_b_tent = 1.0f;
  float b_cdf[PATHSPACE_MAX_VERTS];
  if(tent->length > 2)
  {
    hslt_get_b_cdf(tent, 1, b_cdf);
    b = sample_cdf(b_cdf, tent->length, points_rand(rt.points, tid));
    pdf_b_tent = b ? b_cdf[b] - b_cdf[b-1] : b_cdf[0];
  }
  tent->mmlt.pt_verts = b; // XXX DEBUG

  float c_cdf[PATHSPACE_MAX_VERTS];
  float pdf_c_tent = 1.0f;
  int c = tent->length - 1;
  if(b < tent->length - 1)
  {
    hslt_get_c_cdf(tent, b, c_cdf);
    c = b + 1 + sample_cdf(c_cdf+b+1, tent->length-b-1, points_rand(rt.points, tid));
    pdf_c_tent = c ? c_cdf[c] - c_cdf[c-1] : c_cdf[0];
  }
  tent->mmlt.lt_verts = c; // XXX DEBUG

  assert(!(tent->v[b].mode & s_specular));
  assert(!(tent->v[c].mode & s_specular));

  const float offset[4][2] = {{-1,0},{1,0},{0,-1},{0,1}};
  float f[4] = {-1.0, -1.0, -1.0, -1.0};
  for(int k=0;k<4;k++)
  {
    tent->sensor.pixel_i = curr->sensor.pixel_i + offset[k][0];
    tent->sensor.pixel_j = curr->sensor.pixel_j + offset[k][1];
    // mutate point on aperture, after mutating outgoing direction (pixel)
    tent->sensor.aperture_x = sample_mutate_rand(tent->sensor.aperture_x, points_rand(rt.points, tid), 0.01f);
    tent->sensor.aperture_y = sample_mutate_rand(tent->sensor.aperture_y, points_rand(rt.points, tid), 0.01f);

    double T = hslt_perturb(curr, tent, b, c, d->stats+tid);
    if(T <= 0.0f) return 0.0f;
    const int to = tent->length - curr->length;

    if(tent->length > 2)
    {
      hslt_get_b_cdf(tent, 1, b_cdf);
      const float pdf_b_curr = (b+to) ? b_cdf[b+to] - b_cdf[b+to-1] : b_cdf[0];
      hslt_get_c_cdf(tent, b+to, c_cdf);
      const float pdf_c_curr = c+to ? c_cdf[c+to] - c_cdf[c+to-1] : c_cdf[0];
      T *= pdf_b_curr * pdf_c_curr / (pdf_b_tent * pdf_c_tent);
    }
  }



  // ==================================


  // compute ray differential in half vector space, sample half vectors to add up to one pixel step.
  const int s = 0, e = path.length-1;
  path.v[s].diffgeo.type = s_pinned_position;
  for(int i=s+1;i<=e;i++)
    path.v[i].diffgeo.type = s_free;
  // compute abc matrices:
  manifold_compute_tangents(&path, s, e);
  if(raydifferentials_compute_rd_h(&path, path.cache.R, s, e)) return;
  path_t backup = path;

  // reverse-engineer mis weight:
  const double mis_weight = value / (path_measurement_contribution_dx(&path) / path_pdf(&path));

  float sum_roughness = 0.f;
  for(int k=s+1;k<e;k++)
    sum_roughness += path.v[k].shading.roughness;

  const float offset[4][2] = {{-1,0},{1,0},{0,-1},{0,1}};
  float f[4] = {-1.0, -1.0, -1.0, -1.0};
  for(int k=0;k<4;k++)
  {
    path = backup;
    float rd_i[3], rd_j[3];
    // call into camera to get precise per-pixel ray differential
    if(raydifferentials_v1(rt.cam, &path, offset[k][0], offset[k][1], rd_i, rd_j))
    {
      fprintf(stderr, "rd fail\n");
      // goto shift_failed;
      continue;
    }
    // express rd_i and rd_j not in world space but in tangent space of v1 (actually only 2d now):
    const float rduv_i[2] = {dotproduct(rd_i, path.v[1].diffgeo.dpdu), dotproduct(rd_i, path.v[1].diffgeo.dpdv)};
    const float rduv_j[2] = {dotproduct(rd_j, path.v[1].diffgeo.dpdu), dotproduct(rd_j, path.v[1].diffgeo.dpdv)};
    float h[2*(e-s+1)];
    h[0] = h[1] = h[2*(e-s)] = h[2*(e-s)+1] = 0.0f;
    // mutate halfvectors
    for(int i=s+1;i<e;i++)
    {
      if(path.v[i].mode & s_volume)
        assert(0); // not supported
      else
      {
        // sample gaussians based on roughness and ray diffs
        if(path.v[i].mode & s_specular)
        {
          // force half vector to some precise (001) to avoid drift
          h[2*(i-s)+0] = 0.0f;
          h[2*(i-s)+1] = 0.0f;
        }
        else
        {
          // convert tangent space of x[1] to half vector space of h[i]
          float hu[2], hv[2];
          mat2_mulv(path.cache.R + 4*i, rduv_i, hu);
          mat2_mulv(path.cache.R + 4*i, rduv_j, hv);
          const float *currh = backup.v[i].diffgeo.h;
          float *tenth = h + 2*(i-s);
          const float s = path.v[i].shading.roughness/sum_roughness;
          tenth[0] = currh[0] + s * hu[0] + s * hv[0];
          tenth[1] = currh[1] + s * hu[1] + s * hv[1];
        }
      }
    }

    // convert half vectors to world space
    const float pdf = halfvec_to_worldspace(&d, &path, h, s, e);
    if(pdf == 0.0f)
    {
      // if(path.length == 3)
      //fprintf(stderr, "to worldspace fail\n");
      goto shift_failed;
    }
    // evaluate measurement contributions for all four paths.
    // compute transfer matrices:
    if(halfvec_compute_raydifferentials(&d, &path, 0, path.length-1))
    {
      // if(path.length == 3)
      // fprintf(stderr, "rd2 fail\n");
      // float col[3] = {10000, 0, 0};
    // filter_accum(path.sensor.pixel_i, path.sensor.pixel_j, col, rt.render->fb, 3, 0, 3);
      goto shift_failed;
    }
    const double JT = fabs(backup.cache.dh_dx / path.cache.dh_dx); // jacobian of the shift mapping.
    f[k] = mis_weight * path_measurement_contribution_dx(&path) * JT / path_pdf(&backup);
    // fprintf(stderr, "step %d px offs %d %d\n", k,
    //     (int)(path.sensor.pixel_i - backup.sensor.pixel_i),
    //     (int)(path.sensor.pixel_j - backup.sensor.pixel_j));
    // needed?
    // const float vol_rpdf = halfvec_reverse_check(&d->stats[tid], curr, tent, 0, curr->length-1);
#if 0
    float col[3];
    spectrum_p_to_camera(path.lambda, f[k]/5.0, col);
    filter_accum(path.sensor.pixel_i, path.sensor.pixel_j, col, rt.render->fb, 3, 0, 3);
#endif
    // convert to gradient value:
    const float dir[2] = {path.sensor.pixel_i-p->sensor.pixel_i, path.sensor.pixel_j-p->sensor.pixel_j};
    const float dlen = hypotf(dir[0], dir[1]);
    f[k] = (f[k] - value)/dlen * (dir[0]*offset[k][0] + dir[1]*offset[k][1])/dlen;
    // fprintf(stderr, "f[%d] = %g = %g %g\n", k, f[k], dlen,  (dir[0]*offset[k][0] + dir[1]*offset[k][1])/dlen);

    // splat gradient buffers:
    float col[3] = {0};
    float acc = 0.0f;
    if(0)
    { // non-symmetric case, offset path was occluded or similar. gradient is full contribution of base path:
shift_failed:;
      float len = offset[k][0] + offset[k][1];
      acc = len * value;
      // float col[3] = {10000, 0, 0};
    // filter_accum(path.sensor.pixel_i, path.sensor.pixel_j, col, rt.render->fb, 3, 0, 3);
    }
    else
    {
      acc = .5f * f[k];
    }
    spectrum_p_to_camera(p->lambda, acc, col);
    filter_accum(
        p->sensor.pixel_i + .5f*offset[k][0],
        p->sensor.pixel_j + .5f*offset[k][1],
        col,
        (offset[k][0] == 0) ?
        rt.render->grad_y :
        rt.render->grad_x,
        3, 0, 3);
  }
#if 0
  float gx[3], gy[3];
  spectrum_p_to_camera(p->lambda, .5f*(f[1] - f[0]), gx);
  spectrum_p_to_camera(p->lambda, .5f*(f[3] - f[2]), gy);
  // spectrum_p_to_camera(p->lambda, 1.0f, gx);
  // spectrum_p_to_camera(p->lambda, 1.0f, gy);
  if(f[0] > 0.0 && f[1] > 0.0)
    filter_accum(p->sensor.pixel_i, p->sensor.pixel_j, gx, rt.render->grad_x, 3, 0, 3);
  if(f[3] > 0.0 && f[2] > 0.0)
    filter_accum(p->sensor.pixel_i, p->sensor.pixel_j, gy, rt.render->grad_y, 3, 0, 3);
#endif
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
  if(((sample+1) % (rt.width*rt.height)) == 0) rt.render->overlays ++;
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
