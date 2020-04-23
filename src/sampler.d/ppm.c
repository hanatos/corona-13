/*
    This file is part of corona-13.

    corona-13 is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    corona-13 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with corona-13. If not, see <http://www.gnu.org/licenses/>.
*/
#include "view.h"
#include "sampler.h"
#include "pointsampler.h"
#include "accel.h"
#include "knn.h"

// progressive photon map, just to test the knn queries.

// how many samples in one lookup
#define SAMPLER_KNN 8

typedef struct sampler_t
{
  path_t *eyepath;
  knn_point_t *photons;
  uint32_t photon_cnt;
  heap_t **knn;            // one per thread for the lookup
}
sampler_t;

sampler_t *sampler_init()
{
  const int num_knn = SAMPLER_KNN; // how many samples in one lookup of density estimation
  sampler_t *s = (sampler_t *)malloc(sizeof(sampler_t));
  memset(s, 0, sizeof(*s));
  s->knn = (heap_t **)malloc(sizeof(heap_t*)*rt.num_threads);
  for(int k=0;k<rt.num_threads;k++)
    s->knn[k] = heap_init(num_knn);
  return s;
}

void sampler_prepare_frame(sampler_t *s)
{
  // TODO: will have to think about memory here.
  const int num_samples = view_width()*view_height()/4;
  if(!s->photons)
  {
    s->photons = (knn_point_t *)malloc(sizeof(knn_point_t)*num_samples);
    s->photon_cnt = 0;
    s->eyepath = (path_t *)malloc(sizeof(path_t)*num_samples);
    if(!s->eyepath)
    {
      fprintf(stderr, "[ppm] failed to allocate %.02f MB for a %d sample photon map!\n", sizeof(path_t)*num_samples/(1024.0*1024.0), num_samples);
      exit(1);
    }
  }
  // sppm:
  s->photon_cnt = 0;
  if(s->photon_cnt == 0)
  {
  // cast rays from the eye
    // TODO: implement in parallel via task scheduler for pthreads!
  for(int k=0;k<num_samples;k++)
  {
    path_init(s->eyepath + k, k, 0);
    if(path_extend(s->eyepath + k)) continue;
#if 0
    do
    {
      path_extend(s->eyepath + k);
      if(s->eyepath[k].v[s->eyepath[k].length-1].throughput <= 0.0f) break;
      if(s->eyepath[k].v[s->eyepath[k].length-1].flags & s_environment) break;
      if(s->eyepath[k].v[s->eyepath[k].length-1].mode & s_emit) break;
    }
    while(s->eyepath[k].v[s->eyepath[k].length-1].mode & s_specular);
#endif
    if(s->eyepath[k].v[s->eyepath[k].length-1].flags & s_environment) continue;
    if(s->eyepath[k].v[s->eyepath[k].length-1].mode & s_emit) continue;
    // atomic ++
    const int pindex = __sync_add_and_fetch(&s->photon_cnt, 1);
    knn_point_t *p = s->photons + pindex;
    p->payload = k;
    for(int i=0;i<3;i++)
      p->pos[i] = s->eyepath[k].v[s->eyepath[k].length-1].hit.x[i];
    p->pos[3] = 0.0f;//s->eyepath[k].lambda;
    p->pos[4] = 0.0f;//s->eyepath[k].time;
  }
  // build kd tree
  knn_build_rec(s->photons, 0, s->photon_cnt);
  // fprintf(stderr, "[ppm] have %d photons\n", s->photon_cnt);
  }
}

void sampler_cleanup(sampler_t *s)
{
  for(int k=0;k<rt.num_threads;k++)
    heap_cleanup(s->knn[k]);
  free(s->knn);
  free(s->photons);
  free(s);
}


void sampler_clear(sampler_t *s)
{
  s->photon_cnt = 0;
}

void sampler_create_path(path_t *path)
{
  // fraction of bounding box that we deem acceptable bias.
  // this will only affect the maximum bias and thus impact performance a lot.
  // need to make sure it's larger than a pixel footprint though, we're not casting more
  // importons than that (else your picture turns black)
  const float *aabb = accel_aabb(rt.accel);
  const float max_radius = MAX(aabb[3] - aabb[0],
      MAX(aabb[4] - aabb[1], aabb[5] - aabb[2]))*1e-2f;

  path->v[0].mode = s_emit; // start at the light
  const int tid = common_get_threadid();
  while(1)
  {
    if(path_extend(path)) return;
    const int v = path->length-1;

    // find a couple of importon bins for density estimation:
    float query[5];
    for(int k=0;k<3;k++) query[k] = path->v[v].hit.x[k];
    query[3] = 0.0f;//path->lambda;
    query[4] = 0.0f;//path->time;
    knn_find(rt.sampler->photons, rt.sampler->photon_cnt, query, rt.sampler->knn[tid], max_radius);
    // iterate over them, merge what doesn't belong together, and accumulate to framebuffer.
    sampler_t *s = rt.sampler;
    float maxrad2 = 0;
    if(s->knn[tid]->end)
      for(int k=0;k<3;k++) maxrad2 += (s->photons[s->knn[tid]->keys[0]].pos[k] - query[k])*(s->photons[s->knn[tid]->keys[0]].pos[k] - query[k]);
    for(int k=0;k<s->knn[tid]->end;k++)
    {
      // const float dist = s->knn[tid]->vals[k]; // squared cartesian distance (x,y,z,lambda,time)
      const uint64_t p = s->knn[tid]->keys[k];
      const uint32_t path_index = s->photons[p].payload;
      // TODO: do a cartesian product density estimation (disk radius for space, something similar for lambda and time
      // const float density = 1.0f/(M_PI*maxrad2*s->knn[tid]->end*s->photon_cnt);
      const float density = 1.0f/(M_PI*maxrad2*s->photon_cnt);
      const float throughput = path_merge(path, s->eyepath + path_index);
      // TODO: set path's lambda and time to cached one? or not?
      // path->lambda = s->eyepath[path_index].lambda; // more decorrelated but effectively white
      pointsampler_splat(path, density * throughput);
      // restore path, clip connected portion:
      path->length = v+1;
    }
  }
}

void sampler_print_info(FILE *fd)
{
  fprintf(fd, "sampler  : progressive photon map\n");
}
