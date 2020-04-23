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

#include "corona_common.h"
#include "pathspace.h"
#include "pathspace/tech.h"
#include "pathspace/nee.h"
#include "points.h"
#include "pointsampler.h"
#include "sampler.h"
#include "sampler_common.h"
#include "shader.h"
#include "spectrum.h"
#include "view.h"
#include "camera.h"
#include "ext/halton/halton.h"

#include <stdio.h>
#include <float.h>

// #define DEFERRED_EXPLORATION

#define POINTSAMPLER_MAX_RANDS_PER_VERT 10
#define POINTSAMPLER_NUM_DIMS 2*(PATHSPACE_MAX_VERTS*POINTSAMPLER_MAX_RANDS_PER_VERT)
#define POINTSAMPLER_NUM_MAX_MUTATIONS 313 // lattice has 7, 23, 97 or 313 points
#define POINTSAMPLER_NUM_LATTICE_DIM 21

// defines for specific confs
#ifndef randomtiles
#define randomtiles 1
#ifndef RANDSET
#define RANDSET 0
#endif
#endif
#ifndef skip_nee
#define skip_nee 0
#endif
//DOF_cutoff value to set for cutoff on DOF and circle of confusion
//BREAKUP_POINT
//MINV, MAXV minimum and maximum vertex to choose from for breakup point
//USE_HALTON sample pss paths' random numbers with halton points
#define showptype 0
#ifndef debug_DOF
#define debug_DOF 0
#endif
#define debug_caustic 0


typedef struct pointsampler_thr_t
{
  int mutated_samples;
  float *rand_buf;     // buffers to be swapped.
  float *tent_rand;    // tentative state
  float *curr_rand;    // current (old) sample

  // mutated paths data
  float *pixels;       // pixel positions buffer
  double *pdf;         // pdfs buffer
  double *value;       // throughput buffer
  double *mis_weight;  // mis weight of perturbed
  float *lambda;

  // stats
  int failed;
  int total_samples;

  // per-frame cranley patterson shift
  float cp_shift[100];

  path_t nee_cache;    // cache bifurcations in splat until accept
  uint8_t nee_valid[PATHSPACE_MAX_VERTS+1];
}
pointsampler_thr_t;

typedef struct pointsampler_t
{
  pointsampler_thr_t *t;
  halton_t h;
}
pointsampler_t;

// lut for generator vectors for n=7,23,97 and dimension s=1..21
// as constructed by dirk nuyens' fastrank1pt routine%
// (c) 2004, <dirk.nuyens@cs.kuleuven.ac.be>
static inline const int *get_generator(const int n, const int dim)
{
  assert(dim > 1 && dim <= 21);
  assert(n > 1 && n <= 313);

  int num = 0;
  switch(n)
  {
    case   7: num = 0; break;
    case  23: num = 1; break;
    case  97: num = 2; break;
    case 313: num = 3; break;
    default:
      assert(0);
      return 0;
  }

  static const int gen[4][21][21] = {
    { // n=7
    {1}, // dummy, something co-prime with 23 will result in equispaced numbers (i guess i'm saying anything would)
    {1,3},
    {1,3,2},
    {1,3,2,1},
    {1,3,2,1,3},
    {1,3,2,1,3,2},
    {1,3,2,3,1,2,3},
    {1,3,2,1,2,3,1,3},
    {1,3,2,1,2,3,2,1,3},
    {1,3,2,3,1,2,3,1,2,3},
    {1,3,2,2,1,3,2,3,1,2,1},
    {1,3,2,1,3,2,1,2,3,2,3,1},
    {1,3,2,3,1,2,3,2,1,2,1,3,3},
    {1,3,2,1,3,2,1,3,2,3,2,1,3,2},
    {1,3,2,3,2,1,3,2,1,3,2,1,3,2,1},
    {1,3,2,1,3,2,1,3,2,1,3,2,3,1,2,1},
    {1,3,2,3,2,1,3,1,2,1,3,2,3,1,2,3,1},
    {1,3,2,1,3,2,1,3,2,2,1,3,2,1,3,1,3,2},
    {1,3,2,2,3,1,1,2,3,1,3,2,1,3,2,1,3,2,3},
    {1,3,2,3,1,2,1,3,2,1,2,3,1,3,2,3,1,2,3,1},
    {1,3,2,1,3,2,3,1,2,3,1,2,1,3,2,1,3,2,2,1,3},
    },
    { // n=23
    {1},
    {1,10},
    {1,7,5},
    {1,7,5,3},
    {1,10,4,7,5},
    {1,10,4,7,5,6},
    {1,7,5,3,11,4,10},
    {1,10,4,7,5,6,8,3},
    {1,10,4,7,6,9,2,5,3},
    {1,10,4,7,6,5,8,3,11,2},
    {1,7,5,3,4,11,10,2,8,9,6},
    {1,10,4,7,6,5,8,3,11,9,2,7},
    {1,7,5,3,4,11,10,2,8,6,9,5,4},
    {1,7,5,4,3,6,9,11,2,10,8,10,8,6},
    {1,10,4,6,9,7,2,5,3,8,11,4,5,3,7},
    {1,10,4,6,7,9,2,5,3,8,11,3,2,8,11,9},
    {1,10,4,6,9,8,11,7,5,2,3,9,6,1,10,4,8},
    {1,10,4,6,7,9,2,5,3,8,11,2,9,10,8,6,11,5},
    {1,7,5,4,3,11,10,2,8,9,6,11,5,2,3,8,9,4,10},
    {1,7,5,4,3,11,10,2,8,6,9,2,3,8,11,9,10,7,6,1},
    {1,10,4,6,7,9,2,5,3,11,8,11,8,9,2,10,6,5,1,4,7},
    },
    { // n=97
    {1},
    {1,36},
    {1,36,45},
    {1,36,21,45},
    {1,35,41,23,10},
    {1,35,41,23,10,19},
    {1,36,21,45,28,5,39},
    {1,36,21,28,17,31,9,47},
    {1,35,41,10,18,13,24,4,29},
    {1,36,21,28,31,17,9,47,23,43},
    {1,35,41,10,18,13,24,4,29,47,7},
    {1,35,41,10,18,13,24,4,29,7,47,16},
    {1,36,21,28,31,17,9,6,16,39,26,23,47},
    {1,36,21,28,31,17,9,6,16,39,26,23,47,43},
    {1,36,21,28,31,17,9,6,16,39,26,23,47,43,10},
    {1,35,41,10,18,13,34,37,16,11,38,6,14,33,40,9},
    {1,36,21,28,31,17,37,26,6,10,8,22,19,24,15,27,29},
    {1,36,21,31,28,17,37,26,6,10,8,22,19,24,15,27,29,32},
    {1,36,21,31,28,17,37,26,6,10,8,22,19,24,15,27,29,32,47},
    {1,36,21,31,28,17,37,26,6,10,8,22,19,24,15,27,29,32,47,20},
    {1,36,21,31,28,17,37,26,6,10,8,22,19,24,15,27,29,32,47,20,7},
    },
    { // n = 313
    {1},
    {1,121},
    {1,121,82},
    {1,121,82,23},
    {1,121,82,23,65},
    {1,119,55,80,90,69},
    {1,121,82,23,65,102,135},
    {1,121,82,23,65,102,135,18},
    {1,121,82,23,65,102,76,57,35},
    {1,121,82,23,65,102,76,57,35,108},
    {1,121,82,107,47,19,127,131,143,138,113},
    {1,121,82,107,47,19,127,131,143,138,113,23},
    {1,119,55,100,41,70,89,61,115,146,12,47,80},
    {1,119,55,100,41,70,89,61,115,146,12,47,80,138},
    {1,119,55,100,41,70,89,61,115,146,12,47,80,138,18},
    {1,119,55,100,41,70,89,61,115,146,12,47,80,138,18,22},
    {1,121,82,107,47,19,127,131,143,138,113,53,23,109,13,155,43},
    {1,119,91,64,136,38,145,33,126,47,12,20,110,15,87,115,25,43},
    {1,121,56,81,133,97,17,76,91,53,113,84,149,63,115,143,105,118,43},
    {1,119,91,64,136,38,145,33,126,141,6,73,116,14,89,100,148,69,60,21},
    {1,119,91,64,136,38,145,33,126,141,6,73,116,14,89,100,148,69,60,8,21},
    },
  };

  return gen[num][dim-1];
}

// needs the buffer with random numbers, the number of rands you want to change
// and their ids using pathspace.h notation, plus a dist vector to fill
// will return the most close sample id on the lattice
static inline int get_closest(
    const float* tile_rand, // the s-dimensional kelemen point in tile space
    const int r1_n,       // number n of samples in lattice    
    const int r1_s,         // number of dimensions in rank-1 lattice
    float* r1_dist)         // output shortest distance vector
{
#if RANDSET
  memset(r1_dist,0,sizeof(float)*r1_s);
  return 0;
#endif
  float tmp[r1_s], dist = FLT_MAX;
  int id = -1;
  const int* gen = get_generator(r1_n, r1_s);

  // go through the points sampled
  for(int i=0; i<r1_n; i++)
  {
    float tmpdist = 0;

    // compute the distance
    for(int s=0; s<r1_s; s++)
    {
      const float ptmp = i*gen[s]/(float)r1_n;
      assert(ptmp >= 0);
      float p = ptmp - (int)ptmp;
      tmp[s] = p - tile_rand[s];

      // wrap around if needed for modulo 1 dist
      if(tmp[s] > .5f) tmp[s] -= 1;
      else if(tmp[s] < -.5f) tmp[s] += 1;

      tmpdist += tmp[s]*tmp[s];
      if(tmpdist > dist) break;
    }

    // swap buffers
    if(tmpdist < dist)
    {
      id = i;
      dist = tmpdist;
      for(int s=0; s<r1_s; s++) r1_dist[s] = tmp[s];
    }
  }

  assert(id >= 0);
  return id;
}

static inline void shift_point(
    const int i,          // index of rank-1 lattice point
    float* tile_rand,     // overwrite tile coordinate point
    const int r1_n,       // number n of samples in lattice
    const int r1_s,       // number s of dimensions in lattice
    const float* r1_dist) // cranley patterson rotation vector
{
  const int* gen = get_generator(r1_n, r1_s);
  for(int s=0; s<r1_s; s++)
  {
#if RANDSET
    tile_rand[s] = points_rand(rt.points, common_get_threadid()); continue;
#endif    
      
    const float ptmp = i*gen[s]/(float)r1_n;
    assert(ptmp >= 0);
    const float p = ptmp-(int)ptmp;

    tile_rand[s] = p-r1_dist[s];

    // modulo operation
    if     (tile_rand[s] >= 1) tile_rand[s] -= 1;
    else if(tile_rand[s] <  0) tile_rand[s] += 1;
    assert(tile_rand[s] <= 1.f); // may round to 1.0f unfortunately
    assert(tile_rand[s] >= 0.f);
  }  
}

// transforms unit torus random numbers to tile coordinate system
static inline void torus_to_tile(
    const float *pss_rand,    // input kelemen random buffer
    const int   *r1_to_pss,   // rank-1 to kelemen mapping
    const int    r1_s,        // number of dimensions in rank-1 lattice
    const int   *num_tiles,   // tile table
    float       *tile_offset, // output tile offsets
    float       *tile_rand)   // output tile coordinates
{
  for(int s=0;s<r1_s;s++)
  {
    const int pss = r1_to_pss[s];
    // shift tile offset by global cp shift
    const float cp_shift = rt.pointsampler->t[common_get_threadid()].cp_shift[s];
    const int tile = (pss_rand[pss]+cp_shift / num_tiles[s]) * num_tiles[s];
    tile_offset[s] = (tile - cp_shift) / num_tiles[s];
    tile_rand[s] = (pss_rand[pss] - tile_offset[s]) * num_tiles[s];
    // do we have to modulo here as well?
    if     (tile_rand[s] >= 1) tile_rand[s] -= 1;
    else if(tile_rand[s] <  0) tile_rand[s] += 1;
    assert(tile_rand[s] <= 1.f); // float res rounds (!= 0) + 1 = 1
    assert(tile_rand[s] >= 0.f);
  }
}

// transforms tile to unit torus
static inline void tile_to_torus(
    float* pss_rand,           // output kelemen numbers
    const int *r1_to_pss,      // map rank-1 to kelemen dimensions
    const int  r1_s,           // number of dimensions in lattice
    const int *num_tiles,      // tile numbers per r1 dimension
    const float *tile_offset,  // tile offsets for each r1 dimension
    const float *tile_rand)    // input random numbers in tile space
{
  for(int s=0; s<r1_s; s++)
  {
    const int pss = r1_to_pss[s];
    pss_rand[pss] = tile_offset[s] + tile_rand[s]/num_tiles[s];
    // shifted tiles by global cp shift may exceed these bounds and require modulo
    if(pss_rand[pss] >= 1.0f) pss_rand[pss] -= 1.0f;
    if(pss_rand[pss] <  0.0f) pss_rand[pss] += 1.0f;
    assert(pss_rand[pss] <= 1.f); // should be < but float resolution is lower near 1 than near 0
    assert(pss_rand[pss] >= 0.f);
  }
}

static void print_reprj_info()
{
  // random printing, it works fine :D
  if(common_get_threadid() == 0 && (rt.pointsampler->t[0].total_samples % 20000) == 0)
  {
    int failed = 0, total = 0;
    for(int i=0; i<rt.num_threads; i++)
    {
      failed += rt.pointsampler->t[i].failed;
      total += rt.pointsampler->t[i].total_samples;
    }
    printf("failed: %d total: %d average: %.02f\n", failed, total, 1.f - failed/(float)total);
  }
}

void pointsampler_print_info(FILE *f)
{
  fprintf(f, "mutations: local quasi-Monte Carlo\n");
}

pointsampler_t *pointsampler_init(uint64_t frame)
{
  pointsampler_t *s = (pointsampler_t *)malloc(sizeof(pointsampler_t));
  s->t = (pointsampler_thr_t *)malloc(sizeof(pointsampler_thr_t)*rt.num_threads);
  for(int k=0;k<rt.num_threads;k++)
  {
    s->t[k].rand_buf = (float *)malloc(2*sizeof(float)*POINTSAMPLER_NUM_DIMS);
    memset(s->t[k].rand_buf, 0, 2*sizeof(float)*POINTSAMPLER_NUM_DIMS);
    s->t[k].tent_rand = s->t[k].rand_buf;
    s->t[k].curr_rand = s->t[k].rand_buf + POINTSAMPLER_NUM_DIMS;
    s->t[k].curr_rand[s_dim_lambda] = -0.666f;
    // mutated paths data
    s->t[k].pixels = (float*)malloc(2*sizeof(float)*POINTSAMPLER_NUM_MAX_MUTATIONS);
    memset(s->t[k].pixels, -1, 2*sizeof(float)*POINTSAMPLER_NUM_MAX_MUTATIONS);
    s->t[k].pdf = malloc(sizeof(double)*POINTSAMPLER_NUM_MAX_MUTATIONS);
    memset(s->t[k].pdf, 0, sizeof(double)*POINTSAMPLER_NUM_MAX_MUTATIONS);
    s->t[k].value = malloc(sizeof(double)*POINTSAMPLER_NUM_MAX_MUTATIONS);
    memset(s->t[k].value, 0, sizeof(double)*POINTSAMPLER_NUM_MAX_MUTATIONS);
    s->t[k].mis_weight = malloc(sizeof(double)*POINTSAMPLER_NUM_MAX_MUTATIONS);
    memset(s->t[k].mis_weight, 0, sizeof(double)*POINTSAMPLER_NUM_MAX_MUTATIONS);
    s->t[k].lambda = malloc(sizeof(float)*POINTSAMPLER_NUM_MAX_MUTATIONS);
    memset(s->t[k].lambda, 0, sizeof(float)*POINTSAMPLER_NUM_MAX_MUTATIONS);    

    // stats
    s->t[k].failed = 0;
    s->t[k].total_samples = 0;
  }
  // init halton points
  halton_init_random(&s->h, frame);
  return s;
}

void pointsampler_finalize(pointsampler_t *s) {}

void pointsampler_clear()
{
  for(int k=0;k<rt.num_threads;k++)
  {
    memset(rt.pointsampler->t[k].rand_buf, 0, 2*sizeof(float)*POINTSAMPLER_NUM_DIMS);
    memset(rt.pointsampler->t[k].pixels, -1, 2*sizeof(float)*POINTSAMPLER_NUM_MAX_MUTATIONS);
    memset(rt.pointsampler->t[k].pdf, 0, sizeof(double)*POINTSAMPLER_NUM_MAX_MUTATIONS);
    memset(rt.pointsampler->t[k].value, 0, sizeof(double)*POINTSAMPLER_NUM_MAX_MUTATIONS);    
    memset(rt.pointsampler->t[k].mis_weight, 0, sizeof(float)*POINTSAMPLER_NUM_MAX_MUTATIONS); 
    memset(rt.pointsampler->t[k].lambda, 0, sizeof(float)*POINTSAMPLER_NUM_MAX_MUTATIONS);     
  }
}

void pointsampler_cleanup(pointsampler_t *s)
{
  for(int k=0;k<rt.num_threads;k++)
  {
    free(s->t[k].rand_buf);
    free(s->t[k].pixels);
    free(s->t[k].pdf);
    free(s->t[k].value);
    free(s->t[k].mis_weight);
    free(s->t[k].lambda);
  } 
  free(s->t);
  free(s);
}

float pointsampler(path_t *p, const int i)
{
  int v = p->length;
  const int tid = common_get_threadid();
  pointsampler_thr_t *t = rt.pointsampler->t + tid;
  const int cnt = POINTSAMPLER_MAX_RANDS_PER_VERT;
  const int end = p->v[v].rand_beg;
  if(end + i >= t->mutated_samples)
  {
    const int old = t->mutated_samples;
    t->mutated_samples = end + 3*cnt;
    if(t->mutated_samples > PATHSPACE_MAX_VERTS*cnt) t->mutated_samples = PATHSPACE_MAX_VERTS*cnt;
    for(int k=old;k<t->mutated_samples;k++)
    {
      // random points
      #ifndef USE_HALTON
        t->tent_rand[k] = points_rand(rt.points, tid);
      // halton points
      #else
        const int dim = end + k;
        if(dim >= halton_get_num_dimensions())
        {
          // degenerate to pure random mersenne twister
          t->tent_rand[k] = points_rand(rt.points, tid);
        }
        else
          t->tent_rand[k] = halton_sample(&rt.pointsampler->h, dim, p->index);
      #endif
    }
  }
  assert(end + i < POINTSAMPLER_NUM_DIMS);
  return t->tent_rand[end + i];
}

typedef enum path_type_t
{
  s_type_normal = 0,
  s_type_pure_specular = 1,
  s_type_caustic = 2,
  s_type_nee = 3,
  s_type_dof = 4,
}
path_type_t;

// breakup and samples on lattice

// maxv is the max number of vertices we can handle with our generated lattice
static void getBreakup(const path_t *p, const int minv, const int maxv, float *cdf, path_type_t* type)
{
  // by roughness only
  float rough[PATHSPACE_MAX_VERTS] = {1.0};
  // last index is used for the likeliness to connect from v[end-1] to v[end]
  // in case maxv == p->length-1 we can retrace the whole path
  int end = p->length-1 > maxv ? maxv+1 : p->length-1;
  for(int k=1; k<=end; k++)
  {
    if(p->v[k].flags & s_environment) rough[k] = 1.0;
    else if(p->v[k].mode & s_volume)
    {
      const float g = p->v[k].interior.mean_cos;
      rough[k] = sqrtf(1.0f - g*g);
    }
    else
    {
      int indexmatched = ((p->v[k].mode & s_transmit) && fabsf(path_eta_ratio(p, k) - 1.0f) < 1e-3f);
      float step = p->v[k].shading.roughness;
      rough[k] = (indexmatched ? 0 : step);
    }
  }
  
  for(int k=0; k<end; k++)
  {
    cdf[k] = MIN(rough[k],rough[k+1]);
  }

  if(p->v[p->length-1].tech == s_tech_nee) *type = s_type_nee;
  // get if caustic
  {
    int lastb = 0;
    float specval = 0.05;
    for(int i=1; i<end; i++)
    {
      //get when we have the last good bpoint
      // cdf is not normalised yet
      if(cdf[i] > specval) lastb = i;
    }
    // needs to have AT LEAST 1 specular vertex and be in the breakpoint range
    if(lastb < end-1) *type = s_type_caustic;
    if(!lastb)
    { //in this case we want to trace the whole path and cheat a bit
      if(end == p->length-1) end++;
      cdf[end-1] = 50.f; //if all specular, we want to recreate the whole path
    }
  }

  #ifndef SKIP_DOF
  {
    /*get camera
    get f/n = diameter
    focus -> plane
    get interpolation to find the diameter on the point
    cutoff to decide whether is sharp or not
    */

    const camera_t* cam = view_get_camera(p);
    const float view_f_stop[] = {0.5, 0.7, 1.0, 1.4, 2, 2.8, 4, 5.6, 8, 11, 16, 22, 32, 45, 64, 90, 128};
    float diameter = cam->focal_length/view_f_stop[cam->aperture_value];
    float delta = fabsf(cam->focus - p->e[1].dist);
    float cutoff = 0.1;
    #ifdef DOF_cutoff
    cutoff = DOF_cutoff;
    #endif
    if(diameter*delta/cam->focus > cutoff) //diameter is in mm, delta/focus is just a ratio (but they are in dm)
    {
      if(end == p->length) cdf[end-1] = 1.f; //we were retracing the whole path here
      *type = s_type_dof;
      cdf[0] = 5.f;
      #if 1 //we really want to have almost always 0
      cdf[0] = 80.f;
      #endif
    }
  }
  #endif

  // put to zero prob of getting vertices < minv, used for special cases
  for(int k=0;k<minv;k++) cdf[k] = 0.0f;

  for(int k=1;k<end;k++)
    cdf[k] += cdf[k-1];
  for(int k=0;k<end;k++)
    cdf[k] /= cdf[end-1];
}

// currently unused
// static const float *getLatticeSamples(const float v)
// {
//   static const float r1[] = {0.5,0.79,0.99,1};
//   static const float r7[] = {0.15,0.63,0.88,1};
//   static const float r23[] = {0.1,0.38,0.98,1};
//   static const float r97[] = {0.05,0.1,0.2,1};

//  const float logval = logf(v+1.0f)/logf(2.0f)/2.0f;

// //  printf("%f\n",logval);

//  if(logval < 1.5) return r1;
//  if(logval < 3.0) return r7;
//  if(logval < 4.0) return r23;
//  return r97;
// }

// get mapping from reduced rank-1 to full primary sample space.
// also init tile number array.
static int get_r1_to_pss(
    const path_t *path,
    int *breakup,
    int *r1_to_pss,
    int *num_tiles,
    int *img_id,
    const path_type_t ptype)
{
  const int r1_max_s = POINTSAMPLER_NUM_LATTICE_DIM;
  // run-time configuration of perturbation:
  const int mut_lambda = 1;
  const int mut_time = 0;
  const int mut_stereo = 0;

  // set tilesize based on path
  const int ts_aperture = (ptype == s_type_caustic ? 4 : 2);
  const int ts_image = path->e[1].vol.shader >= 0 ? 12 : (ptype == s_type_caustic ? 48 : 24);
  const int ts_free_path = (ptype == s_type_caustic ? 155 : 50);

  // always perturbed, even if only v0
  r1_to_pss[0] = s_dim_aperture_x;
  num_tiles[0] = ts_aperture;
  r1_to_pss[1] = s_dim_aperture_y;
  num_tiles[1] = ts_aperture;
  int s = 2; // number of dimensions in lattice
  if(mut_lambda) r1_to_pss[s] = s_dim_lambda;
  if(mut_lambda) num_tiles[s++] = 1;
  if(mut_time)   r1_to_pss[s] = s_dim_time;
  if(mut_time)   num_tiles[s++] = 1;
  if(mut_stereo) r1_to_pss[s] = s_dim_camid;
  if(mut_stereo) num_tiles[s++] = 1;

  if(*breakup == 0) return s;

  // perturb v1
  r1_to_pss[s] = s_dim_image_x;
  img_id[0] = s;
  num_tiles[s++] = ts_image;
  r1_to_pss[s] = s_dim_image_y;
  img_id[1] = s;
  num_tiles[s++] = (ts_image * view_height())/view_width();
  if(path->e[1].vol.shader >= 0)
  {
    r1_to_pss[s] = path->v[1].rand_beg + s_dim_free_path;
    num_tiles[s++] = ts_free_path;
  }

  // perturb v2++ via outgoing direction
  for(int i=2;i<=*breakup; i++)
  {
    const int olds = s;
    int ts_omega = (ptype == s_type_caustic ? 55 : 155);
    // if traversing a volume
    if(path->e[i].vol.shader >= 0)
    {
      r1_to_pss[s] = path->v[i].rand_beg + s_dim_free_path;
      num_tiles[s++] = 1000;//ts_free_path;
      ts_omega = 1000; //do we really need 1k?
    }
    // if is specular just forget
    if(!(path->v[i].mode & s_specular))
    {
      r1_to_pss[s] = path->v[i].rand_beg + s_dim_omega_x;
      num_tiles[s++] = ts_omega;
      r1_to_pss[s] = path->v[i].rand_beg + s_dim_omega_y;
      num_tiles[s++] = ts_omega;
    }
    if(s>r1_max_s) {*breakup = i-1; return olds;}
  }
  assert(s <= r1_max_s);
  return s;
}

static void swap_tiles(
  int *r1_to_pss,
  int *num_tiles,
  int *img_id,
  const int size
)
{
  // randomly choose two dimensions in the hypercube
  // we're fine with every combination
  int np[2] = {
    points_rand(rt.points,common_get_threadid())*(size-1),
    points_rand(rt.points,common_get_threadid())*(size-1)
  };
  
  // swap dimensions
  for(int k=0; k<2; k++)
  {
    int tmp = r1_to_pss[img_id[k]];
    r1_to_pss[img_id[k]] = r1_to_pss[np[k]];
    r1_to_pss[np[k]] = tmp;

    tmp = num_tiles[img_id[k]];
    num_tiles[img_id[k]] = num_tiles[np[k]];
    num_tiles[np[k]] = tmp;

    img_id[2+k] = np[k];
  }
}

// currently unused
// static void unswap_tiles(
//   float *rand_buff,
//   const int *r1_to_pss,
//   const int *img_id
// )
// {
//   // swap dimensions
//   for(int k=0; k<2; k++)
//   {
//     int mapto = r1_to_pss[img_id[k]];
//     int mapfrom = r1_to_pss[img_id[k+2]];

//     float tmp = rand_buff[mapto];
//     rand_buff[mapto] = rand_buff[mapfrom];
//     rand_buff[mapfrom] = tmp;
//   }  
// }

static void explore(path_t *path, float value)
{
  const int debug_cond = 0;

  float bcdf[PATHSPACE_MAX_VERTS];
  path_type_t ptype = s_type_normal;

  const int tid = common_get_threadid();
  pointsampler_thr_t *t = rt.pointsampler->t + tid;

  //const int r1_n = getLatticeSamples(ptype); // number of lattice samples: 7,23,97
  // max number of dimensions
  const int r1_max_s = POINTSAMPLER_NUM_LATTICE_DIM;
  int camid[POINTSAMPLER_NUM_MAX_MUTATIONS] = {0};
  // map local rank-1 lattce point index to kelemen space
  int r1_to_pss[r1_max_s+3];
  int num_tiles[r1_max_s+3];
  int maxv = path->length -1;
  int img_id[4] = {0}; //indices to swap image space dimensions
  // final number s of dimensions in our rank-1 lattice
  int r1_s = get_r1_to_pss(
      path, &maxv, r1_to_pss, num_tiles, img_id, ptype);

  int minv = 0;
#ifdef MINV
  minv = MINV;
#endif
#ifdef MAXV
  maxv = MAXV;
#endif

  getBreakup(path, minv, maxv, bcdf, &ptype);

  if(path->v[path->length-1].tech == s_tech_nee && ptype != s_type_nee) return;

#if showptype == 1

  float col[3] = {0};
  float fval = skip_nee ? 10.f : 0.03f;
  if(ptype == s_type_normal) for(int i=0; i<3; i++) col[i] = fval;
  else if (ptype == s_type_caustic) col[0] = 10000*fval;
  else if (ptype == s_type_nee) col[1] = fval;
  else if (ptype == s_type_dof) col[2] = fval;

  view_splat_col(path, col); return;
#endif

#ifndef BREAKUP_POINT
  const float brand = points_rand(rt.points,tid);
  const int breakup = sample_cdf(bcdf, maxv+1, brand);

// debug flags
#if debug_DOF
  if(ptype == s_type_dof) view_splat(path,value);return;
#elif debug_caustic
  if(ptype == s_type_caustic) view_splat(path,value);return;
#endif  

  // maybe we have a nicer way of updating dimensionality with respect to breakup?
  // with uniform random this should not be necessary anymore
  maxv = breakup;
  if(ptype != s_type_normal) r1_s = get_r1_to_pss(path, &maxv, r1_to_pss, num_tiles, img_id, ptype);
  // lattice points number selection
  // const float nrand = points_rand(rt.points,tid);
  // const float* ncdf = getLatticeSamples(value);
  // const int nbin = sample_cdf(ncdf, path->length, nrand);
  // const int r1_n = 23;//(nbin == 0 ? 1 : (nbin == 1 ? 7 : (nbin == 2 ? 23 : 97)));
  // optimisation for lattice points number selection
  // if(r1_n == 1)
  // {
  //   view_splat(path, value/ncdf[0]);
  //   return;
  // }
  const int r1_n = (ptype == s_type_nee ? 7 : (ptype == s_type_caustic ? 97 : 23));

#else // compute ref images
  int breakup = BREAKUP_POINT; // path->length-1;
  if(breakup >= path->length) breakup = path->length-1;
  if(path->v[path->length-1].tech == s_tech_nee && breakup == path->length-1)
    breakup = path->length-2; // min length is 3, so we can do this.
  const int r1_n = 23;//313;//23;
#endif

  if((path->v[breakup].mode & s_specular) ||
    ((breakup < path->length-1) && (path->v[breakup+1].mode & s_specular)))
  {
    view_splat(path, value);
    return;
  }

  // swap curr <-> tent buffers
  float* tmp = t->curr_rand;
  t->curr_rand = t->tent_rand;
  t->tent_rand = tmp;

  path_t perpath;

  int failed = 0;
  double pdfsum = 0;

  float r1_dist[r1_max_s];
  float tile_offset[r1_max_s];
  float tile_rand[r1_max_s];

  swap_tiles(r1_to_pss, num_tiles, img_id, r1_s);
  torus_to_tile(t->curr_rand, r1_to_pss, r1_s, num_tiles, tile_offset, tile_rand);
  const int seed = get_closest(tile_rand, r1_n, r1_s, r1_dist);

// splat seed path
#ifndef DEBUG_SEED_PATH
  t->pixels[2*seed+0] = path->sensor.pixel_i;
  t->pixels[2*seed+1] = path->sensor.pixel_j;
  // t->pdf[seed] = path_pdf(path);
  t->value[seed] = path_measurement_contribution_dx(path, 0, path->length-1);
  t->mis_weight[seed] = sampler_mis_weight(path);
  t->lambda[seed] = path->lambda;
  camid[seed] = path->sensor.camid;
  assert(t->value[seed] < DBL_MAX);
  t->pdf[seed] = 1;
  for(int i=0; i<breakup+1; i++)
    t->pdf[seed] *= path_pdf_extend(path,i);
#endif

  float bprob = 1;//(breakup ? bcdf[breakup] - bcdf[breakup-1] : bcdf[breakup]);
  float nprob = 1;//ncdf[nbin] - ncdf[nbin-1];

  // assert(bprob<=1.f);
  // assert(nprob<=1.f);
  // assert(bprob>=0.f);
  // assert(nprob>=0.f);

// #define common
#ifdef common
  double commonpdf = 1;
  for(int i=breakup+3; i<path->length; i++)
    commonpdf *= path_tech_vertex_pdf_as_sampled(path,i);
  assert(commonpdf);

#ifndef DEBUG_SEED_PATH
  pdfsum = commonpdf*path_tech_vertex_pdf_as_sampled(path,breakup+1)*(breakup+2<path->length?path_tech_vertex_pdf_as_sampled(path,breakup+2):1);
#endif
#else
#ifndef DEBUG_SEED_PATH
  pdfsum += path_pdf(path)/t->pdf[seed]*bprob*nprob;
#endif
#endif
  if(t->pdf[seed] == 0) pdfsum = 0;
  if(!(pdfsum < DBL_MAX)) return;
  assert(pdfsum < DBL_MAX);
int printend = 0;
  memcpy(t->tent_rand, t->curr_rand, sizeof(float)*PATHSPACE_MAX_VERTS*POINTSAMPLER_MAX_RANDS_PER_VERT);

  // get current random numbers in tile space:
  torus_to_tile(t->curr_rand, r1_to_pss, r1_s, num_tiles, tile_offset, tile_rand);

  for(int i=0; i<r1_n; i++)
  {
#ifndef DEBUG_SEED_PATH
//     // don't recreate the seed path:
    if(i == seed) continue;
#endif

    path_init(&perpath, path->index, path->sensor.camid);
    perpath.tangent_frame_scrambling = path->tangent_frame_scrambling;

    // shift according to the rank-1 lattice including c/p:
    shift_point(i, tile_rand, r1_n, r1_s, r1_dist);
    // transform back to kelemen space
    tile_to_torus(t->tent_rand, r1_to_pss, r1_s, num_tiles, tile_offset, tile_rand);

    if(i == seed)
      for(int j=0; j<PATHSPACE_MAX_VERTS*POINTSAMPLER_MAX_RANDS_PER_VERT; j++)
        assert(fabs(t->tent_rand[j] - t->curr_rand[j]) < 1e-6);

    while(breakup && perpath.length <= breakup)
    {
      if(perpath.length)
      { // need to overwrite that here:
        perpath.v[perpath.length-1].rand_cnt = path->v[perpath.length-1].rand_cnt;
        perpath.v[perpath.length-1].rand_beg = path->v[perpath.length-1].rand_beg;
      }
      if(path_extend(&perpath))
      {
        if(debug_cond) fprintf(stderr, "failed extend\n");
        goto fail;
      }
      const int v = perpath.length-1;
      if(breakup < path->length-1)
      {
        if(perpath.v[v].mode & s_emit) 
        {
          if(debug_cond) fprintf(stderr, "failed, hit light\n");
          goto fail;
        }
        // todo check if we really want to kill those or they always generate the same set
        // if we start with less vertices but the breakup point is the same should be fine
        if(perpath.v[v].flags & s_environment)
        {
          if(debug_cond) fprintf(stderr, "failed envmap\n");
          goto fail;
        }
      }
      else
      {
        if(perpath.v[v].flags != path->v[v].flags)
        {
          if(debug_cond) fprintf(stderr, "failed last envmap %d %d flag %d %d\n", v, path->length, perpath.v[v].flags, path->v[v].flags);
          goto fail;
        }
      }

      if(perpath.v[v-1].mode != path->v[v-1].mode)
      {
        if(debug_cond) fprintf(stderr, "failed mode %d %d mode %d %d\n", v, path->length, perpath.v[v-1].mode, path->v[v-1].mode);
        goto fail;
      }
    }

    if(breakup == 0)
    {
      // connecting directly to aperture again.
      // we won't be calling path_extend() or anything that inits the path, so
      // we do this now manually:
      perpath.sensor = path->sensor;
      perpath.tangent_frame_scrambling = path->tangent_frame_scrambling;
      perpath.lambda = spectrum_sample_lambda(pointsampler(&perpath, s_dim_lambda), 0);
      perpath.time = view_sample_time(&perpath, pointsampler(&perpath, s_dim_time));
      perpath.sensor.aperture_x = pointsampler(&perpath, s_dim_aperture_x);
      perpath.sensor.aperture_y = pointsampler(&perpath, s_dim_aperture_y);
      perpath.sensor.camid = view_sample_camid(pointsampler(&perpath, s_dim_camid));
      perpath.e[0] = path->e[0];
      perpath.v[0] = path->v[0];
    }

    // connect to rest seed path if we still have a connecting segment:
    if(breakup < path->length-1)
    {
      // check that we can actually sample this event before reprojecting
      if((perpath.v[breakup].material_modes & path->v[breakup].mode) == 0)
      {
        if(debug_cond) fprintf(stderr, "failed breakup modes to reconnect event\n");
        goto fail;
      }

      // first, copy next vertex from seed path
      perpath.e[breakup+1] = path->e[breakup+1];
      perpath.v[breakup+1] = path->v[breakup+1];
      perpath.v[breakup].mode = path->v[breakup].mode;

      // try to connect
      perpath.length = breakup + 2;
      if(path->v[breakup+1].flags & s_environment)
      {
        // keep direction but be sure to move outside aabb
        for(int k=0;k<3;k++)
          perpath.v[breakup+1].hit.x[k] = perpath.v[breakup].hit.x[k]
            + 2.0*(path->v[breakup+1].hit.x[k] - path->v[breakup].hit.x[k]);
      }

      if(path_project(&perpath, breakup+1, s_propagate_mutate) ||
         perpath.v[breakup+1].flags != path->v[breakup+1].flags)
      {
        if(debug_cond) fprintf(stderr, "failed project\n");
        goto fail;
      }
      if(!(path->v[breakup+1].flags & s_environment))
      { // all non envmap vertices need to be sufficiently close:
        const float *oldx = path->v[breakup+1].hit.x;
        const float eps = MAX(MAX(1.0f, fabsf(oldx[0])), MAX(fabsf(oldx[1]), fabsf(oldx[2])));
        if(fabsf(perpath.v[breakup+1].hit.x[0] - oldx[0]) > 1e-3f*eps ||
           fabsf(perpath.v[breakup+1].hit.x[1] - oldx[1]) > 1e-3f*eps ||
           fabsf(perpath.v[breakup+1].hit.x[2] - oldx[2]) > 1e-3f*eps) 
           {
             if(debug_cond) fprintf(stderr, "failed envmap, not close enough\n");
             goto fail;
           }
      }

      if(perpath.v[breakup].mode != path->v[breakup].mode)
      {
        if(debug_cond) fprintf(stderr, "failed breakup mode\n");
        goto fail;
      }      

      // fill the rest
      for(int v=breakup+2; v<path->length; v++)
      {
        perpath.v[v] = path->v[v];
        perpath.e[v] = path->e[v];
      }
      // update for new wavelength:
      for(int v=breakup+2; v<path->length; v++)
        shader_prepare(&perpath, v);


      if(perpath.v[breakup].tech != 1)
        fprintf(stderr, "%d\n", perpath.v[breakup].tech);

    }
    perpath.length = path->length;

    t->pixels[i*2+0] = perpath.sensor.pixel_i;
    t->pixels[i*2+1] = perpath.sensor.pixel_j;
    t->lambda[i] = perpath.lambda;
    camid[i] = perpath.sensor.camid;
    t->value[i] = path_measurement_contribution_dx(&perpath,0,perpath.length-1);
    if(perpath.v[breakup].mode != path->v[breakup].mode) 
    {
      if(debug_cond) fprintf(stderr, "failed breakup mode\n");
      goto fail;
    }
    t->mis_weight[i] = sampler_mis_weight(&perpath);
    if(!(t->value[i] > 0.0))
    {
      if(debug_cond) fprintf(stderr, "failed measurement\n");
      goto fail;
    }
    assert(t->value[i] < DBL_MAX);
    for(int v=0; v<path->length; v++) perpath.v[v].tech = path->v[v].tech;

    t->pdf[i] = 1.0;
    for(int j=0; j<breakup+1; j++)
      t->pdf[i] *= path_pdf_extend(&perpath, j);
    if(!(t->pdf[i] > 0.0)) {
      if(debug_cond) fprintf(stderr, "failed pdf jacobian\n");
      goto fail;
    }

    // check if path respects its group properties
    path_type_t ntype;
    // update cdf for the new path
    float ncdf[32];
    getBreakup(&perpath, minv, maxv, ncdf, &ntype);

    if((ptype == s_type_dof || ntype == s_type_dof) && ntype != ptype) {
      if(debug_cond) fprintf(stderr, "failed path type\n");
      goto fail;
    }

#if showptype == 2

  if(ntype != ptype) goto fail;

  float col[3] = {0};// {0.03, 0.03, 0.03};
  float fval = skip_nee ? .3f : 0.03f;
  if(ptype == s_type_normal) for(int i=0; i<3; i++) col[i] = fval;
  else if (ptype == s_type_caustic) for(int i=0; i<3; i++) col[i] = fval;//col[0] = fval;
  // else if (ptype == s_type_nee) col[1] = fval;
  // else if (ptype == s_type_dof) col[2] = fval;

  // const float[3] col = {fcol[0], fcol[1], fcol[2]};
  view_splat_col(&perpath, col); continue;
#endif

  // const float* nncdf = getLatticeSamples(t->mis_weight[i]*t->value[i]/path_tech_pdf_as_sampled(&perpath));
  float bbprob = 1;
  // compute prob for this path to have taken the same decision as seed
  // int nbk = sample_cdf(bcdf, path->length, brand);
  // bbprob = (nbk ? bcdf[nbk] - bcdf[nbk-1] : bcdf[nbk]);
  // bbprob /= (breakup ? bcdf[breakup] - bcdf[breakup-1] : bcdf[breakup]);
  // int nnbin = sample_cdf(nncdf, path->length, nrand);
  float nnprob = 1;//(nbin? nncdf[nbin] - nncdf[nbin-1] : nncdf[nbin])/(nnbin? ncdf[nnbin] - ncdf[nnbin-1] : ncdf[nnbin]);
  // float bbprob = (nbk ? bcdf[nbk] - bcdf[nbk-1] : bcdf[nbk])/bprob;
  // float bbprob = (nbk ? bcdf[nbk] - bcdf[nbk-1] : bcdf[nbk])/(breakup ? bcdf[breakup] - bcdf[breakup-1] : bcdf[breakup]);
  // float bbprob = (breakup ? bcdf[breakup] - bcdf[breakup-1] : bcdf[breakup])/bprob;
  //if(nbk != breakup) printf("%d %d %d\n",i,breakup,nbk);
  // if(nnbin != nbin) printf("%d ",nnbin);
  // if(nbk != breakup) {if(!printend) printf("%d ",breakup); printend = 1; printf("%d ", nbk);}
  
  // assert(bbprob<=1.f);
  // assert(nnprob<=1.f);
  // assert(bbprob>=0.f);
  // assert(nnprob>=0.f);

  #ifdef common
    pdfsum += commonpdf*path_tech_vertex_pdf_as_sampled(&perpath,breakup+1)*(breakup+2<path->length?path_tech_vertex_pdf_as_sampled(&perpath,breakup+2):1);
  #else
    pdfsum += path_tech_pdf_as_sampled(&perpath)/t->pdf[i]*bbprob*nnprob;
  #endif
    assert(pdfsum < DBL_MAX);

    if(0)
    {
      fail:
      failed++;

      if(i==seed)
      {
        printf("failed to recreate the seed path\n");
        path_print(&perpath, stdout);
        assert(0);
      }

      t->failed++;
      t->pdf[i] = 0;
      t->value[i] = 0;
      t->pixels[i*2+0] = perpath.sensor.pixel_i;
      t->pixels[i*2+1] = perpath.sensor.pixel_j;
      t->lambda[i] = perpath.lambda;      
    }
  }
  t->total_samples += r1_n;

  if(debug_cond) printf("%d\n",failed);
  if(debug_cond) print_reprj_info();

  assert(pdfsum < DBL_MAX);
if(printend) printf("\n");
  // do MIS and splat
  for(int i=0; i<r1_n; i++)
  {

    if(t->pdf[i] == 0) continue;
    float val = t->mis_weight[i] * t->value[i]/(t->pdf[i]*pdfsum);
    if(!(val < DBL_MAX)) continue; // happens in pdf=1e-36 cases
    assert(val>=0);

    perpath.sensor.pixel_i = t->pixels[2*i];
    perpath.sensor.pixel_j = t->pixels[2*i+1];
    perpath.sensor.camid = camid[i];
    perpath.lambda = t->lambda[i];
    perpath.length = path->length;

    #if !showptype
    view_splat(&perpath, val);
    #endif
  }

  // swap back rand buffers
  tmp = t->curr_rand;
  t->curr_rand = t->tent_rand;
  t->tent_rand = tmp;
}

void pointsampler_splat(path_t *path, float value)
{
  if(!(value > 0.0 && value < FLT_MAX)) return;

#if skip_nee
  if(path->length < 3 || path->v[path->length-1].tech == s_tech_nee)
#else
  if(path->length < 3)
#endif  
  {
    #if !debug_DOF && !debug_caustic && !showptype
    view_splat(path, value);
    #endif
    return;
  }

#ifndef DEFERRED_EXPLORATION
  explore(path, value);
// #else

//   const int breakup = 1;
//   if((path->v[breakup].mode & s_specular) ||
//     ((breakup < path->length-1) && (path->v[breakup+1].mode & s_specular)))
//   {
//     view_splat(path, value);
//     return;
//   }

//   pointsampler_t *s = rt.pointsampler;
//   const int tid = common_get_threadid();
//   pointsampler_thr_t *t = s->t + tid;
  
//   if(path->v[path->length-1].tech == s_tech_nee)
//   { // remember nee vertex, it will be popped away
//     t->nee_cache.e[path->length-1] = path->e[path->length-1];
//     t->nee_cache.v[path->length-1] = path->v[path->length-1];
//     t->nee_valid[path->length-1] = 1;
//   }
//   else t->nee_valid[PATHSPACE_MAX_VERTS] = 1; // indicates valid scattering path
#endif
}

int pointsampler_accept(path_t *curr, path_t *tent)
{
#ifndef DEFERRED_EXPLORATION
  return 0;
#else
//   // we'll need a few temporary things. as we explore and splat a group of
//   // paths, they need each one the storage as explore, but some can be shared:
//   const int r1_max_s = POINTSAMPLER_NUM_LATTICE_DIM;
// #ifdef VEACH_DOOR
//   const int r1_n = 313;  // number of points in the lattice
// #else
//   const int r1_n = 23;   // number of points in the lattice
// #endif
//   float pixel_i[r1_n];   // shared for all perturbations:
//   float pixel_j[r1_n];   // pixel position and
//   float lambda[r1_n];    // wavelength.
//   double value[r1_n][PATHSPACE_MAX_VERTS+1];
//   double pdfsum[PATHSPACE_MAX_VERTS+1] = {0.0};
//   double jacobian[PATHSPACE_MAX_VERTS+1] = {0.0};
//   double measurement[PATHSPACE_MAX_VERTS+1] = {0.0};
//   double perpdf[PATHSPACE_MAX_VERTS+1] = {0.0};
//   uint8_t valid[r1_n][PATHSPACE_MAX_VERTS+1];

//   // map local rank-1 lattce point index to kelemen space
//   int r1_to_pss[r1_max_s+3];
//   int num_tiles[r1_max_s+3];
//   const int tid = common_get_threadid();
//   pointsampler_thr_t *t = rt.pointsampler->t + tid;
//   path_t *path = tent;

// #if 1
//   float bcdf[PATHSPACE_MAX_VERTS];
//   path_type_t ptype = s_type_normal;
//   int maxv = path->length -1;
//   // final number s of dimensions in our rank-1 lattice
//   int r1_s = get_r1_to_pss(
//       path, &maxv, r1_to_pss, num_tiles, ptype);
//       // assert(maxv == path->length -1);
//   int minv = 0;
// #ifdef VEACH_DOOR
//   minv = 1;
// #endif
//   getBreakup(path, minv, maxv, bcdf, &ptype);
  
//   const int breakup = sample_cdf(bcdf, maxv+1, points_rand(rt.points,tid));
//   maxv = breakup;
//   if(ptype != s_type_normal) r1_s = get_r1_to_pss(path, &maxv, r1_to_pss, num_tiles, ptype);
//   // const int r1_n = 313;//(ptype == s_type_nee ? 7 : (ptype == s_type_caustic ? 97 : 23));
// #else
//   // determine breakup point
// #ifndef BREAKUP_POINT
// #define BREAKUP_POINT 0
// #endif
//   int breakup = BREAKUP_POINT; // path->length-1;
//   if(breakup >= path->length) breakup = path->length - 1;
//   if(path->v[path->length-1].tech == s_tech_nee && breakup == path->length-1)
//     breakup = path->length-2; // min length is 3, so we can do this.
//   const int r1_s = get_r1_to_pss(
//       path, &breakup, r1_to_pss, num_tiles, 0);
// #endif

//   // now explore all paths we got.
//   // check whether there are any at all
//   int got_paths = t->nee_valid[PATHSPACE_MAX_VERTS]; // that's the scattering one
//   for(int k=2;k<path->length;k++)
//     got_paths |= t->nee_valid[k];
//   if(!got_paths) return 0; // gah.

//   if((path->v[breakup].mode & s_specular) ||
//     ((breakup < path->length-1) && (path->v[breakup+1].mode & s_specular)))
//     return 0; // handled in pointsampler_splat() already.

//   // swap curr <-> tent buffers
//   float* tmp = t->curr_rand;
//   t->curr_rand = t->tent_rand;
//   t->tent_rand = tmp;

//   path_t perpath;

//   float r1_dist[r1_max_s];
//   float tile_offset[r1_max_s];
//   float tile_rand[r1_max_s];

//   torus_to_tile(t->curr_rand, r1_to_pss, r1_s, num_tiles, tile_offset, tile_rand);
//   const int seed = get_closest(tile_rand, r1_n, r1_s, r1_dist);

//   memcpy(t->tent_rand, t->curr_rand, sizeof(float)*PATHSPACE_MAX_VERTS*POINTSAMPLER_MAX_RANDS_PER_VERT);
//   for(int i=0; i<r1_n; i++)
//   {
//     for(int k=0;k<path->length;k++) valid[i][k] = t->nee_valid[k];
//     valid[i][PATHSPACE_MAX_VERTS] = t->nee_valid[PATHSPACE_MAX_VERTS];

//     // clear contributions:
//     for(int j=0; j<=path->length; j++)
//       value[i][j] = jacobian[j] = measurement[j] = perpdf[j] = 0.0;

// #ifndef DEBUG_SEED_PATH
//     // don't recreate the seed path:
//     if(i == seed)
//     {
//       perpath.sensor = path->sensor;
//       perpath.tangent_frame_scrambling = path->tangent_frame_scrambling;
//       perpath.lambda = path->lambda;
//       perpath.time = path->time;
//       perpath.length = path->length;
//       for(int k=0;k<path->length;k++)
//       {
//         perpath.e[k] = path->e[k];
//         perpath.v[k] = path->v[k];
//       }
//     }
//     else
// #endif
//     { // reconstruct path from random numbers:
//       path_init(&perpath, path->index, path->sensor.camid);
//       perpath.tangent_frame_scrambling = path->tangent_frame_scrambling;

//       // get current random numbers in tile space:
//       torus_to_tile(t->curr_rand, r1_to_pss, r1_s, num_tiles, tile_offset, tile_rand);
//       // shift according to the rank-1 lattice including c/p:
//       shift_point(i, tile_rand, r1_n, r1_s, r1_dist);
//       // transform back to kelemen space
//       tile_to_torus(t->tent_rand, r1_to_pss, r1_s, num_tiles, tile_offset, tile_rand);

//       int need_connect = 1;
//       if(breakup == 0)
//       {
//         // connecting directly to aperture again.
//         // we won't be calling path_extend() or anything that inits the path, so
//         // we do this now manually:
//         perpath.sensor = path->sensor;
//         perpath.tangent_frame_scrambling = path->tangent_frame_scrambling;
//         perpath.lambda = spectrum_sample_lambda(pointsampler(&perpath, s_dim_lambda), 0);
//         perpath.time = view_sample_time(&perpath, pointsampler(&perpath, s_dim_time));
//         perpath.sensor.aperture_x = pointsampler(&perpath, s_dim_aperture_x);
//         perpath.sensor.aperture_y = pointsampler(&perpath, s_dim_aperture_y);
//         perpath.sensor.camid = view_sample_camid(pointsampler(&perpath, s_dim_camid));
//         perpath.e[0] = path->e[0];
//         perpath.v[0] = path->v[0];
//         perpath.length = 1;
//       }
//       else while(perpath.length <= breakup)
//       {
//         if(perpath.length)
//         { // need to overwrite that here:
//           perpath.v[perpath.length-1].rand_cnt = path->v[perpath.length-1].rand_cnt;
//           perpath.v[perpath.length-1].rand_beg = path->v[perpath.length-1].rand_beg;
//         }
//         if(path_extend(&perpath) ||
//            ((perpath.v[perpath.length-1].mode & s_emit) != (path->v[perpath.length-1].mode & s_emit)) ||
//            (perpath.v[perpath.length-2].mode != path->v[perpath.length-2].mode) ||
//            (perpath.v[perpath.length-1].flags != path->v[perpath.length-1].flags) ||
//            ((perpath.e[perpath.length-1].dist < FLT_MAX) != (path->e[perpath.length-1].dist < FLT_MAX)))
//         {
//           const int v = perpath.length;
//           for(int k=v;k<path->length;k++) valid[i][k] = 0; // kill bifurcations after this
//           valid[i][PATHSPACE_MAX_VERTS] = 0; // kill scattering contribution

//           // keep path length where it is, v[v-1] can still have a valid nee replacement.
//           // just don't call project to connect next segment:
//           need_connect = 0;
          
//           int got_paths = 0;
//           for(int k=2;k<v;k++) got_paths |= valid[i][k];
//           if(!got_paths) goto fail;
//           break;
//         }
//       }

//       // connect to rest seed path if we still have a connecting segment,
//       // but only if path extension didn't fail above:
//       if(need_connect && breakup < path->length-1)
//       {
//         perpath.v[breakup].mode = path->v[breakup].mode;

//         // try to connect
//         perpath.length = breakup + 2;
//         perpath.v[breakup+1] = path->v[breakup+1];
//         if(path->v[breakup+1].flags & s_environment)
//         {
//           // keep direction but be sure to move outside aabb
//           for(int k=0;k<3;k++)
//             perpath.v[breakup+1].hit.x[k] = perpath.v[breakup].hit.x[k]
//               + 2.0*(path->v[breakup+1].hit.x[k] - path->v[breakup].hit.x[k]);
//         }

//         perpath.e[breakup+1] = path->e[breakup+1];
//         if(path_project(&perpath, breakup+1, s_propagate_mutate) ||
//             perpath.v[breakup+1].flags != path->v[breakup+1].flags)
//         { // restore vertex but mark as dead
//           perpath.v[breakup+1] = path->v[breakup+1];
//           perpath.v[breakup+1].mode = s_absorb;
//           // perpath length stays at breakup + 2
//           for(int k=breakup+2;k<path->length;k++) valid[i][k] = 0; // kill bifurcations after this
//           valid[i][PATHSPACE_MAX_VERTS] = 0; // kill scattering contribution
//           int got_paths = 0;
//           for(int k=2;k<breakup+2;k++) got_paths |= valid[i][k];
//           if(!got_paths) goto fail;
//         }
//         if(!(path->v[breakup+1].flags & s_environment))
//         { // all non envmap vertices need to be sufficiently close:
//           const float *oldx = path->v[breakup+1].hit.x;
//           const float eps = MAX(MAX(1.0f, fabsf(oldx[0])), MAX(fabsf(oldx[1]), fabsf(oldx[2])));
//           if(fabsf(perpath.v[breakup+1].hit.x[0] - oldx[0]) > 1e-3f*eps ||
//              fabsf(perpath.v[breakup+1].hit.x[1] - oldx[1]) > 1e-3f*eps ||
//              fabsf(perpath.v[breakup+1].hit.x[2] - oldx[2]) > 1e-3f*eps)
//           {
//             for(int k=breakup+2;k<path->length;k++) valid[i][k] = 0; // kill bifurcations after this
//             valid[i][PATHSPACE_MAX_VERTS] = 0; // kill scattering contribution
//             int got_paths = 0;
//             for(int k=2;k<breakup+2;k++) got_paths |= valid[i][k];
//             if(!got_paths) goto fail;
//           }
//         }
//         else
//         {
//           // fill the rest
//           for(int v=breakup+2; v<path->length; v++)
//           {
//             perpath.v[v] = path->v[v];
//             perpath.e[v] = path->e[v];
//           }
//           // update for new wavelength:
//           for(int v=breakup+2; v<path->length; v++)
//             shader_prepare(&perpath, v);
//           perpath.length = path->length;
//         }
//       }
//     } // path reconstructed
//     assert(perpath.length <= path->length);

//     // these are shared among all three sub paths
//     pixel_i[i] = perpath.sensor.pixel_i;
//     pixel_j[i] = perpath.sensor.pixel_j;
//     lambda [i] = perpath.lambda;

//     for(int v=0; v<MIN(path->length,breakup+1); v++)
//       perpath.v[v].tech = path->v[v].tech;

//     // compute the jacobians from truncated primary sample space to path space
//     jacobian[0] = 1.0;
//     measurement[0] = view_cam_eval(&perpath);
//     perpdf[0] = path_tech_vertex_pdf_as_sampled(&perpath, 0);
//     for(int j=1; j<=perpath.length; j++)
//     {
//       // precompute jacobians:
//       if(j <= breakup && j < perpath.length)
//       {
//         jacobian[j] = path_pdf_extend(&perpath, j);
//         jacobian[j] *= jacobian[j-1];
//       }
//       else jacobian[j] = jacobian[j-1];

//       // precompute measurement contribution and sampling pdf
//       if(j < perpath.length)
//       {
//         if(measurement[j-1] > 0.0)
//           measurement[j] = measurement[j-1] * path_G(&perpath, j) * shader_vol_transmittance(&perpath, j);
//         else measurement[j] = 0.0;
//         if(j-1) measurement[j] *= shader_brdf(&perpath, j-1);
//         perpdf[j] = perpdf[j-1] * path_tech_vertex_pdf_as_sampled(&perpath, j);

//         if(j-1 == breakup || j-1 == breakup+1)
//         { // brdf eval gives us the real mode flags, check for sanity:
//           if(perpath.v[j-1].mode != path->v[j-1].mode)
//           {
//             for(int k=j+1;k<perpath.length;k++) valid[i][k] = 0;
//             valid[i][PATHSPACE_MAX_VERTS] = 0;
//           }
//         }
//       }
//       else if(j == perpath.length) 
//       {
//         measurement[j] = measurement[j-1] * lights_eval_vertex(&perpath, j-1);
//         perpdf[j] = perpdf[j-1];
//       }
//     }

//     // first record scattering contribution from last vertex without nee, we
//     // will destroy the path later
//     if(valid[i][PATHSPACE_MAX_VERTS] && (jacobian[perpath.length] > 0.0))
//     {
//       value[i][perpath.length] = measurement[perpath.length] *
//         sampler_mis_weight(&perpath) / jacobian[perpath.length];

//       // assert all went well or kill the contribution of this path length:
//       if(!(value[i][perpath.length] > 0.0)) value[i][perpath.length] = 0.0;
//       else // scattering pdf (copied from tentative which would be the scattering one here)
//         pdfsum[perpath.length] += perpdf[perpath.length] / jacobian[perpath.length];
//       if(!(pdfsum[perpath.length] < DBL_MAX) || !(pdfsum[perpath.length] >= 0.0)) pdfsum[perpath.length] = 0.0;
//     }

//     // collect nee branches of the path group, do this by copying over the nee
//     // vertices to the end of the path (which will destroy it)
//     for(int v=perpath.length-1;v>=2;v--)
//     {
//       if(!valid[i][v]) continue;
//       if(!(jacobian[v] > 0.0)) continue;
//       perpath.v[v] = t->nee_cache.v[v];
//       perpath.e[v] = t->nee_cache.e[v];
//       if(v <= breakup + 1)
//       { // these edges are wrong (e.g. for v=2 and breakup=1), fix it:
//         // make sure pdf and transmittances are evaluated correctly
//         perpath.e[v].transmittance = 0.0f;
//         if(path_project(&perpath, v, s_propagate_mutate) ||
//             perpath.v[v].flags != t->nee_cache.v[v].flags) continue;
//         if(!(perpath.v[v].flags & s_environment))
//         { // all non envmap vertices need to be sufficiently close:
//           const float *oldx = t->nee_cache.v[v].hit.x;
//           const float eps = MAX(MAX(1.0f, fabsf(oldx[0])), MAX(fabsf(oldx[1]), fabsf(oldx[2])));
//           if(fabsf(perpath.v[v].hit.x[0] - oldx[0]) > 1e-3f*eps ||
//              fabsf(perpath.v[v].hit.x[1] - oldx[1]) > 1e-3f*eps ||
//              fabsf(perpath.v[v].hit.x[2] - oldx[2]) > 1e-3f*eps)
//             continue;
//         }
//       }
//       shader_prepare(&perpath, v); // update for new wavelength (XXX may not be necessary after project()?)
//       perpath.length = v+1;
//       const double fw = measurement[v-1] * shader_brdf(&perpath, v-1) * path_G(&perpath, v) *
//         shader_vol_transmittance(&perpath, v) * lights_eval_vertex(&perpath, v) * sampler_mis_weight(&perpath);
//       const double pdf = perpdf[v-1] * path_tech_vertex_pdf_as_sampled(&perpath, v);
//       value[i][v] = fw / jacobian[v];

//       if(!(value[i][v] > 0.0)) value[i][v] = 0.0;
//       else pdfsum[v] += pdf / jacobian[v];
//       if(!(pdfsum[v] < DBL_MAX) || !(pdfsum[v] >= 0.0)) pdfsum[v] = 0.0;
//     }
//     if(0)
//     { // exception handling: kill all contributions of this lattice point
// fail:;
//       for(int j=0; j<=path->length; j++)
//         value[i][j] = jacobian[j] = measurement[j] = perpdf[j] = 0.0;
//     }
//   }

//   // do MIS and splat
//   for(int i=0; i<r1_n; i++)
//   {
//     double val = 0.0;
//     for(int v=0;v<=path->length;v++) // collect all fake splats
//       if(value[i][v] > 0.0 && pdfsum[v] > 0.0)
//         val += value[i][v] / pdfsum[v];
//     if(!(val < DBL_MAX)) continue; // happens in pdf=1e-36 cases
//     perpath.sensor.pixel_i = pixel_i[i];
//     perpath.sensor.pixel_j = pixel_j[i];
//     perpath.lambda = lambda[i];
//     perpath.length = path->length; // oh, that's a lie.
//     view_splat(&perpath, val);
//   }

//   // swap back rand buffers
//   tmp = t->curr_rand;
//   t->curr_rand = t->tent_rand;
//   t->tent_rand = tmp;
  return 0;
#endif
}

void pointsampler_mutate(path_t *curr, path_t *tent)
{
  pointsampler_mutate_with_pixel(curr, tent, -1, -1);
}

void pointsampler_mutate_with_pixel(path_t *curr_path, path_t *tent_path, float x, float y)
{
  pointsampler_t *s = rt.pointsampler;
  const int tid = common_get_threadid();
  pointsampler_thr_t *t = s->t + tid;
  t->mutated_samples = 0;

  path_init(tent_path, tent_path->index, tent_path->sensor.camid);
  path_init(&t->nee_cache, tent_path->index, tent_path->sensor.camid);
  memset(t->nee_valid, 0, sizeof(t->nee_valid));
  sampler_create_path(tent_path);
}

void pointsampler_prepare_frame(pointsampler_t *s)
{
  for(int tid=0;tid<rt.num_threads;tid++)
  {
    pointsampler_thr_t *t = s->t + tid;
    // enable/disable tile randomisation:
  #if randomtiles
    for(int k=0;k<POINTSAMPLER_NUM_LATTICE_DIM;k++) t->cp_shift[k] = points_rand(rt.points, tid);    
  #else
    for(int k=0;k<POINTSAMPLER_NUM_LATTICE_DIM;k++) t->cp_shift[k] = 0.0f;
  #endif
  }
}

#undef POINTSAMPLER_MAX_RANDS_PER_VERT
#undef POINTSAMPLER_NUM_DIMS
#undef POINTSAMPLER_NUM_MAX_MUTATIONS
