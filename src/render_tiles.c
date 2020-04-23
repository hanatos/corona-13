#include "render.h"

#define SPP 16

static inline void render_get_pixel(
    int i,
    int *x,
    int *y)
{
  *x = 0;
  *y = 0;
  // enough for a tile wd of 32:
  for(int k=0;k<5;k++)
  {
    *x |= (i & 1) << k; i >>= 1;
    *y |= (i & 1) << k; i >>= 1;
  }
}

void render_accum(path_t *p, const float value)
{
  float e[3];
  if(!isfinite(value)) return;
  spectrum_p_to_cam(p->lambda, value, e);
  filter_accum(p->sensor.pixel_i, p->sensor.pixel_j, e, rt.render->fb, 3, 0, 3);
}

// render path corresponding to given sample index
void render_sample(uint64_t sample)
{
  const int tilewd = 32; // pixel width of a tile

  // pixel to sample this time
  int px, py;
  const int tid = rt_tls.tid;
  int pixel_index = __sync_fetch_and_add(rt.render->pixel_index + tid, 1);
  int tile = rt.render->tile_index[tid];
  const int max_px = SPP*tilewd*tilewd;

  if(pixel_index >= max_px)
  {
    // get new tile. this is reserved for our thread, the others don't change our tile. grab one from the sync list:
    tile = __sync_fetch_and_add(&rt.render->tile, 1);
    if(tile >= rt.render->num_tiles)
    { // all tiles given away, now steal tasks!
      // we know there must be a task for us, so try hard to find it:
      while(1)
      { // loop over the gap described in the else branch
        for(int k=0;k<rt.num_threads;k++)
        {
          if(rt.render->pixel_index[k] <= max_px)
          {
            tile = rt.render->tile_index[k];
            pixel_index = __sync_fetch_and_add(rt.render->pixel_index + k, 1);
            if(pixel_index < max_px)
            { // safety check for threads
              break;
            }
          }
        }
        break; // worst case we reduce the number of threads and let one finish do all the work.
      }
    }
    else
    { // we got our new tile, start from pixel 0
      // fprintf(stderr, "[worker %03d] grabbing tile %d\n", tid, tile);
      pixel_index = 0;
      rt.render->tile_index[tid] = tile;
      rt.render->pixel_index[tid] = 1; // do that last to avoid others grabbing our pixel.
      // this results in a gap between blocking the tile number above when fetching it
      // and releasing it for task stealing in the line here.
    }
  }

  render_get_pixel(pixel_index/SPP, &px, &py);
  px += rt.render->offx[tile];
  py += rt.render->offy[tile];

  // TODO: need leo's halton enum:
  rt_tls.render->tent_path->index = sample;
  pointsampler_mutate_with_pixel(rt_tls.render->curr_path, rt_tls.render->tent_path, px, py);
  if(pointsampler_accept(rt_tls.render->curr_path, rt_tls.render->tent_path))
  {
    path_t *tmp = rt_tls.render->curr_path;
    rt_tls.render->curr_path = rt_tls.render->tent_path;
    rt_tls.render->tent_path = tmp;
  }
}

void render_prepare_frame()
{
  if(rt.render->tile >= rt.render->num_tiles)
  {
    rt.render->tile = 0;
    rt.render->overlays += SPP;
  }
}

// update pixel index k
void render_update(uint64_t k)
{
  // TODO: do that per pixel/tile correctly!
  const float inv_o = rt.cam->iso/100.0f * 1.0f/(rt.render->overlays+SPP);
  // convert camera/framebuffer to display profile:
  float rgb[3];
  for(int i=0;i<3;i++) rgb[i] = rt.render->fb[3*k+i]*inv_o;
  colorspace_cam_to_rgb(rgb, rt.render->pixel + 4*k + 1);
}

#undef SPP
