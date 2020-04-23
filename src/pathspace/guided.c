#include "accel.h"
#include "dbor.h"
#include "svd2.h"
#include "fakegaussian.h"
#include "heap.h"
#include "lbvh.h"
#include "lights.h"
#include "matrix2.h"
#include "matrix2d.h"
#include "matrix3.h"
#include "matrix3d.h"
#include "pathspace/guided.h"
#include "pathspace/raydifferentials.h"
#include "pointsampler.h"
#include "points.h"
#include "sampler_common.h"
#include "shader.h"
#include "spectrum.h"
#include "threads.h"
#include "view.h"
#include "color_palette.h"
#include <float.h>

// configure constant kernel cutoff in multiples of sigma: 2.0f is generally
// lower variance, 4.0f is conservative (i.e. leaves no gaps between gaussians)
#define CUTOFF 2.f

#define USE_PATH_CDF 1
#define USE_DBOR 1
#define RAY_DIFFERENTIALS 1
#define RAY_DIFF_MIN 10
#define RAY_DIFF_MAX 90
#ifndef RAY_DIFF_ADAPTIVE
#define RAY_DIFF_ADAPTIVE 1
#endif
//#define TRUST_RECORD_THR 0.25     // for living room test
#define TRUST_RECORD_THR 0.25     // up this to accept more paths to the cache
#define NUM_NB 10
#define SIMPLE_GAUSSIAN_SAMPLING 0
#ifndef RENDER_WITH_GAUSS
#define RENDER_WITH_GAUSS 1
#endif

typedef struct gpath_t
{
  uint32_t vi;      // start vertex index
  uint32_t len;     // path length
  float lambda;     // wavelength as constructed
  uint32_t nb_cnt;  // number of neighbour paths

  // throughput = mc / pdf
  double measurement_contribution;
  double guided_pdf;
  double uniform_pdf;
  float pdf_connect;
  float pixel_x;
  float pixel_y;
  double weight;
  uint32_t meta;    // meta info
  atomic_uint_fast32_t sampled;
  atomic_uint_fast32_t success;
  float pixel_step;
  int age;
}
gpath_t;

// type: 0 = sampled uniformly 
//       1 = sampled by guiding
static void gpath_set_sampling_type(gpath_t *path, sampling_type_t st)
{
  if (st) path->meta |=  (1 << 0);
  else    path->meta &= ~(1 << 0);
}

static sampling_type_t gpath_get_sampling_type(const gpath_t *path)
{
  if ((path->meta >> 0) & 1) return st_guided;
  return st_uniform;
}

typedef struct gvert_t
{
  float x[3];       // world space position
  uint16_t flags;   // inside or envmap
  uint16_t mode;    // scatter mode (volume, reflect, transmit)
  float sigma_dist; // marginalised sigma for distances
  float pdf_bsdf;
  float p_bsdf;
  float jacobian;
  float roughness;

  float S12S22i[9];
  float E[9];
  float w[3];
  float mean[6];
}
gvert_t;

typedef struct guided_stats_t
{
  uint64_t num_samples;
  uint64_t num_paths;
  uint64_t err_project[16];
  uint64_t err_extend[16];
  uint64_t num_uniform_paths_recorded;
  uint64_t num_guiding_paths_recorded;
  uint64_t num_uniform_paths_in_cache;
  uint64_t num_guiding_paths_in_cache;
}
guided_stats_t;

typedef struct guided_cache_tls_t
{
  heap_t *heap;
  uint32_t *path_list;
  float *path_pdf;
  uint32_t size;
  guided_stats_t stats;
}
guided_cache_tls_t;

typedef struct guided_cache_t
{
  gpath_t *guide_paths;          // cached paths
  gvert_t *guide_verts;          // vertex pool used by cached paths
  uint32_t num_paths, num_verts; // allocation sizes
  uint64_t num_pathverts;        // path + vert counts, for reasons of sync in one reg. these can actually be used to compute pdfs and such
  uint64_t num_new_pathverts;    // updated counts while adding more paths during an iteration

  double *path_block_cdf;        // simple global path sampling
  double *path_cdf;              // stage1 and stage2
  uint32_t path_block_cdf_cnt;
  uint32_t path_cdf_cnt;

  guided_cache_tls_t *tls;
  uint64_t job;
  uint64_t num_jobs;

  lbvh_t bvh;

  // TODO: clean up before release!
  dbor_t *dbor;                  // buffers for density based outlier rejection

  int gauss;                     // which gaussian shape to use: 1 b-spline 2 fakegauss 3 const
  render_mode_t render_mode;
  int* nearest_paths;
  float* nearest_paths_dists;
  int learn_iteration;
}
guided_cache_t;
  
// ================================================================================
//  helpers to access the path cache
// ================================================================================

uint32_t guided_num_guide_paths(const guided_cache_t *c)
{
  return c->num_pathverts >> 32;
}

static gvert_t *get_vert(const guided_cache_t *c, const int pi, const int vi)
{
  // vertex numbers are real vertex numbers like on path_t
  return c->guide_verts + c->guide_paths[pi].vi + vi;
}

static int path_len(const guided_cache_t *c, const int pi)
{
  return c->guide_paths[pi].len;
  // assert(pi >= 0);
  // assert(pi < guided_num_guide_paths(c));
  // if(pi == guided_num_guide_paths(c) - 1)
  //   return (c->num_pathverts&0xffffffffu) - c->guide_paths[pi].vi;
  // return c->guide_paths[pi+1].vi - c->guide_paths[pi].vi;
}

// ================================================================================
//  init and cleanup
// ================================================================================
guided_cache_t *guided_init(
    int num_paths)
{
  guided_cache_t *c = malloc(sizeof(*c));
  memset(c, 0, sizeof(*c));
  c->gauss = 2; // fake gaussian by default and for learning
  c->tls = calloc(rt.num_threads, sizeof(guided_cache_tls_t));
  const uint32_t num_cached_paths = num_paths;
  const uint32_t sqrt_num = sqrtf(num_cached_paths)+1;
  const uint32_t size_block = sqrt_num*sqrt_num;
  c->path_block_cdf = malloc(sizeof(*c->path_block_cdf)*sqrt_num);
  c->path_cdf = malloc(sizeof(*c->path_cdf)*size_block);
  for(int t=0;t<rt.num_threads;t++)
  {
    c->tls[t].size = num_cached_paths;
    c->tls[t].path_pdf = malloc(sizeof(float)*c->tls[t].size);
    c->tls[t].path_list = malloc(sizeof(uint32_t)*c->tls[t].size);
  }

  c->num_paths = num_paths;
  c->num_verts = num_paths * 5; // speculate on low average vertex count
  c->guide_paths = malloc(sizeof(*c->guide_paths)*c->num_paths);
  c->guide_verts = malloc(sizeof(*c->guide_verts)*c->num_verts);
  c->num_pathverts = 0;
  c->num_new_pathverts = 0;

  c->nearest_paths = calloc(NUM_NB*c->num_paths, sizeof(int));
  c->nearest_paths_dists = calloc(NUM_NB*c->num_paths, sizeof(float));

  c->dbor = dbor_init(view_width(), view_height(), 1, 10);
  c->render_mode = rm_learn;
  fprintf(stderr, "[guided] allocating %.02f MB for %d/%d cache paths\n", 
      ((sizeof(gvert_t))*c->num_verts + sizeof(gpath_t)*c->num_paths)/1024.0f/1024.0f, num_cached_paths, num_paths);

  return c;
}

void guided_cleanup(guided_cache_t *c)
{
  free(c->guide_paths);
  free(c->guide_verts);
  lbvh_cleanup(&c->bvh);
  for(int t=0;t<rt.num_threads;t++)
  {
    heap_cleanup(c->tls[t].heap);
    free(c->tls[t].path_list);
    free(c->tls[t].path_pdf);
  }
  free(c->tls);
  free(c->path_block_cdf);
  free(c->path_cdf);
  free(c->nearest_paths);
  free(c->nearest_paths_dists);
  dbor_cleanup(c->dbor);
  free(c);
}

void guided_prepare_frame(
    guided_cache_t *c,
    int num_nb,
    int heap_size)
{
  for(int t=0;t<rt.num_threads;t++)
  {
    if (c->tls[t].heap) heap_cleanup(c->tls[t].heap);
    c->tls[t].heap = heap_init(heap_size);
  }
}

// ================================================================================
//  begin path from the eye (image space adaptive sampling + depth for volumes)
// ================================================================================

static int guided_sample_dist(
    const guided_cache_t *c,
    int pi,           // id of guide path
    int vi,           // vertex id (>0, we'll sample a target around this one)
    float rz,         // random number
    float *dist)
{
  assert(vi>0);
  const gvert_t* vn = get_vert(c, pi, vi);
  const gvert_t* vp = get_vert(c, pi, vi-1);
  assert(vn->mode & s_volume);
  const float sigma = vn->sigma_dist;
  const double d[] = {
    vn->x[0]-vp->x[0],
    vn->x[1]-vp->x[1],
    vn->x[2]-vp->x[2]};

  // sample gaussian:
  double g = 0.0;
  const int tid = common_get_threadid();
  if(c->gauss == 2)
    g = fakegaussian_sample(rz);
  else if(c->gauss == 1)
    g = sample_cubic_bspline(rz, points_rand(rt.points, tid),
        points_rand(rt.points, tid), points_rand(rt.points, tid));
  else
    g = CUTOFF*(1.0f-2.0f*rz); // uniform kernel

  *dist = sqrt(dotproduct(d,d)) + g * sigma;
  if(!(*dist > 0.0f)) return 1;
  return 0;
}

static double guided_volume_gauss_distance_pdf(
    const guided_cache_t *c,
    const float sigma,
    const float xp[3], 
    const float xn[3], 
    const float vp[3], 
    const float vn[3])
{
  // need to evaluate distance sampling from guide path
  const double dp[] = { xn[0]-xp[0], xn[1]-xp[1], xn[2]-xp[2]};
  const double dg[] = { vn[0]-vp[0], vn[1]-vp[1], vn[2]-vp[2]};
  const double lp = sqrt(dotproduct(dp, dp));
  const double lg = sqrt(dotproduct(dg, dg));
  const double dt = (lp-lg)/sigma;
  if(c->gauss == 1)
    return sample_cubic_bspline_pdf(dt) / sigma;
  else if(c->gauss == 2)
    return fakegaussian_pdf(dt) / sigma;
  else if(fabs(dt) <= CUTOFF)
    return 1.0f / (2.0f * CUTOFF * sigma);
  return 0.0;
}


// helper function. find point to shoot towards
// return 0 if all good and 1 if fail (gauss sampled out of pdf range)
static int guided_sample_pos(
    const guided_cache_t *c,
    int pi,           // id of guide path
    int vi,           // vertex id (>0, we'll sample a target around this one)
    float rx,         // random numbers
    float ry,
    float rz,
    const float *xp,  // previous scattering event on path
    float *sampled_x) // return target point
{
  assert(vi>0);
  const gvert_t* vn = get_vert(c, pi, vi);
  const gvert_t* vp = get_vert(c, pi, vi-1);
  const int vol = (vn->mode & s_volume);
  const int env = (vn->flags & s_environment);
  assert(!(env && vol)); // can be either env or volume

  // draw fake gaussian random numbers:
  float g1 = 0, g2 = 0, g3 = 0;
  const int tid = common_get_threadid();
  if(c->gauss == 1) // b-spline
  {
    g1 = sample_cubic_bspline(rx, ry,
        points_rand(rt.points, tid), points_rand(rt.points, tid));
    g2 = sample_cubic_bspline(
        points_rand(rt.points, tid), points_rand(rt.points, tid),
        points_rand(rt.points, tid), points_rand(rt.points, tid));
  }
  else if(c->gauss == 2)
  { // fake gaussian
    g1 = fakegaussian_sample(rx);
    g2 = fakegaussian_sample(ry);
  }
  else
  { // constant kernel
    if(vol)
    { // uniform point in sphere
      sample_sphere(&g1, &g2, &g3, rx, ry);
      const float su = CUTOFF*sqrtf(rz);
      g1 *= su; g2 *= su; g3 *= su;
    }
    else
    { // uniform point in disk
      const float su = sqrtf(rx);
      g1 = CUTOFF * su*cosf(2.f*M_PI*ry);
      g2 = CUTOFF * su*sinf(2.f*M_PI*ry);
      g3 = 0.0f;
    }
  }

  // construct tangent frame from vertex and v-1
  float n[3] = {0};
  for(int k=0;k<3;++k)
    n[k] = env ? vn->x[k] : (vn->x[k]-vp->x[k]);
  normalise(n);
  float u[3], v[3];
  get_onb(n, u, v);

  // take the conditional at xp
  // S' = S11 - S12 S22^-1 S21
  // u' = u1  + S12 S22^-1 (x - u2)
  float delta[3], uc[3];
  delta[0] = dotproduct(u, xp);
  delta[1] = dotproduct(v, xp);
  delta[2] = (vp->mode & s_volume) ? dotproduct(n, xp) : 0.f;
  for(int k=0;k<3;k++) delta[k] -= vn->mean[k+3];
  mat3_mulv(vn->S12S22i, delta, uc);
  for(int k=0;k<3;k++) uc[k] += vn->mean[k];

  float Bi[9] = {0};
  float D[9] = {0};
  D[0] = vn->w[0];
  D[4] = vn->w[1];
  D[8] = vn->w[2];
  mat3_mul(vn->E, D, Bi);

  if (vol) 
  {
    if(c->gauss == 1)
      g3 = sample_cubic_bspline(
          points_rand(rt.points, tid), points_rand(rt.points, tid),
          points_rand(rt.points, tid), points_rand(rt.points, tid));
    else if(c->gauss == 2)
      g3 = fakegaussian_sample(rz);

    float gt[3] = {0}, g[] = { g1, g2, g3 };
    mat3_mulv(Bi, g, gt);
    for(int i=0;i<3;i++)
    {
      sampled_x[i] =
        + (uc[0] + gt[0]) * u[i]
        + (uc[1] + gt[1]) * v[i]
        + (uc[2] + gt[2]) * n[i];
    }
  }
  else
  {
    float gt[3] = {0}, g[] = { g1, g2, 0.f };
    mat3_mulv(Bi, g, gt);
    assert(gt[2] == 0.f);
    for(int i=0;i<3;i++)
    {
      sampled_x[i] = dotproduct(vn->x, n) * n[i]
        + (uc[0] + gt[0]) * u[i]
        + (uc[1] + gt[1]) * v[i];
    }
    if(env) // environment maps need to explicitly add back the point position:
      for(int k=0;k<3;k++) sampled_x[k] += xp[k];
  }

  return 0;
}

// returns the pdf to sample the given next vertex, in solid angle measure (not projected).
// that is, it divides the far cosine and distance squared out of a vertex area measure,
// but leaves the cosine at the previous vertex (because there are no normals on guide paths).
static double guided_pos_pdf(
    const guided_cache_t *c,
    const gvert_t* vp,
    const gvert_t* vn,
    const float *xp,         // previous scattering event on path
    const float *sampled_x)  // for envmaps, pass the direction instead
{
  const int vol = (vn->mode & s_volume);
  const int env = (vn->flags & s_environment);
  assert(!(env && vol)); // can be either env or volume

  // construct tangent frame from vertex and v-1
  float n[3] = {0};
  for(int k=0;k<3;++k) // envmaps don't store world space x but direction:
    n[k] = env ? vn->x[k] : (vn->x[k]-vp->x[k]);
  normalise(n);
  float u[3], v[3];
  get_onb(n, u, v);

  // take the conditional at xp
  // S' = S11 - S12 S22^-1 S21
  // u' = u1  + S12 S22^-1 (x - u2)
  float delta[3], uc[3];
  delta[0] = dotproduct(u, xp);
  delta[1] = dotproduct(v, xp);
  delta[2] = (vp->mode & s_volume) ? dotproduct(n, xp) : 0.f;
  for(int k=0;k<3;k++) delta[k] -= vn->mean[k+3];
  mat3_mulv(vn->S12S22i, delta, uc);
  for(int k=0;k<3;k++) uc[k] += vn->mean[k];
  
  float B[9] = {0};
  float D[9] = {0};
  float Et[9] = {0};
  if (vol)
    mat3_transpose(vn->E, Et);
  else
    mat3_transpose_sub2(vn->E, Et);
  D[0] = 1.f/vn->w[0];
  D[4] = 1.f/vn->w[1];
  D[8] = 1.f/vn->w[2];
  mat3_mul(D, Et, B);

  float sqrtdetSc = vn->w[0] * vn->w[1] * vn->w[2];

  double pdf = 0;
  if (vol)
  {
    float d[3] = {
      dotproduct(sampled_x, u)-uc[0],
      dotproduct(sampled_x, v)-uc[1],
      dotproduct(sampled_x, n)-uc[2]};
    float dt[3] = {0};
    mat3_mulv(B, d, dt);
    if(c->gauss == 1)
    {
      pdf = sample_cubic_bspline_pdf(dt[0])
          * sample_cubic_bspline_pdf(dt[1])
          * sample_cubic_bspline_pdf(dt[2])
          * 1.0/sqrtdetSc;
    }
    else if(c->gauss == 2)
    {
      pdf = fakegaussian_pdf(dt[0])
          * fakegaussian_pdf(dt[1])
          * fakegaussian_pdf(dt[2])
          * 1.0/sqrtdetSc;
    }
    else
    {
      const float dist = dotproduct(dt, dt);
      if(dist <= CUTOFF*CUTOFF) pdf = 3.0f/(4.0f*M_PI*CUTOFF*CUTOFF*CUTOFF);
      pdf *= 1.f/sqrtdetSc;
    }
    // convert 3d vertex area measure to solid angle measure for outgoing direction at xp:
    const float e[3] = {
      sampled_x[0]-xp[0],
      sampled_x[1]-xp[1],
      sampled_x[2]-xp[2]};
    pdf *= dotproduct(e, e);
  }
  else
  {
    float wo[3] = {0}, offset[3], t;
    if(env)
    { // envmaps pass sampled_x as direction
      for(int k=0;k<3;++k) wo[k] = sampled_x[k];
      // the plane is at distance 1:
      t = 1.0f/dotproduct(n, wo);
      if(t <= 0) return 0.0;
      for(int k=0;k<3;k++) offset[k] = t * wo[k];
    }
    else
    {
      for(int k=0;k<3;++k)
        wo[k] = sampled_x[k]-xp[k];
      normalise(wo);
      t = (dotproduct(vn->x, n)-dotproduct(xp, n))/dotproduct(n, wo);
      if (t <= 0) return 0.0;
      for(int k=0;k<3;k++) offset[k] = xp[k] + t * wo[k];
    }
    float d[3] = {
      dotproduct(offset, u)-uc[0],
      dotproduct(offset, v)-uc[1],
      0.f};
    float dt[3] = {0};
    mat3_mulv(B, d, dt);
    if(c->gauss == 1)
    {
      pdf = sample_cubic_bspline_pdf(dt[0])
          * sample_cubic_bspline_pdf(dt[1])
          * 1.0/sqrtdetSc;
    }
    else if(c->gauss == 2)
    {
      pdf = fakegaussian_pdf(dt[0])
          * fakegaussian_pdf(dt[1])
          * 1.0/sqrtdetSc;
    }
    else
    {
      const float dist = dt[0]*dt[0]+dt[1]*dt[1];
      if(dist <= CUTOFF*CUTOFF)
      {
        pdf = 1.0/sqrtdetSc/(M_PI*CUTOFF*CUTOFF);
      }
    }
    pdf *= t*t/fabsf(dotproduct(n, wo));
  }
  return pdf;
}

static float guided_extend_bsdf(
    const path_t *path,
    int v)
{
  if(v <= 1) return 0;
  // roughness will be set to 0 for index matched materials in the dielectric's
  // prepare() shader, so these should return 1.0 here, too.
  const float r = (path->v[v-1].material_modes & s_volume) ?
    fabsf(path->v[v-1].interior.mean_cos) :
    (1.0f-path->v[v-1].shading.roughness);
  return (path->v[v-1].material_modes & s_fiber) ? 0.f : CLAMP(1.5*r-0.5f, 0.0f, 1.0f);
}

int guided_sample_guide_path(guided_cache_t *c)
{
  if(guided_num_guide_paths(c) <= 0) return -1;
#if USE_PATH_CDF==1
  uint32_t blockid = sample_cdfd(c->path_block_cdf, c->path_block_cdf_cnt, points_rand(rt.points, common_get_threadid()));
  return c->path_cdf_cnt*blockid +
    sample_cdfd(c->path_cdf+c->path_cdf_cnt*blockid, c->path_cdf_cnt, points_rand(rt.points, common_get_threadid()));
#else
  const double xi = points_rand(rt.points, common_get_threadid());
  return CLAMP(guided_num_guide_paths(c) * xi, 0, guided_num_guide_paths(c)-1);
#endif
}

#if SPECTRAL_IS==1
static const float mutation_max = 0.5f;
static const float mutation_sigma = 0.5f * 2.f / 8.f;//range of fake gaussian is -4 to 4 (=> divide by 8), range of mutation is [-mutation_max, mutation_max] (=> multiply by 2)

static float guided_mutate(const float x, const float r)
{
  const float g = fakegaussian_sample(r) * mutation_sigma;
  const float xm = x + g;
  const int xi = floorf(xm);
  assert(xm - xi >= 0.0);
  assert(xm - xi <= 1.0);
  return xm - xi;
}

static float guided_mutate_pdf(const float x, const float y)
{
  const float d1 = fabsf(x - y), d2 = 1.0f-d1;
	const float g1 = fakegaussian_pdf(d1/mutation_sigma);
  const float g2 = fakegaussian_pdf(d2/mutation_sigma);
  return (g1 + g2)/mutation_sigma;
}
#endif

int guided_sample(guided_cache_t *c, path_t *path, uint32_t pathid)
{
  guided_stats_t *stats = &c->tls[common_get_threadid()].stats;
  stats->num_samples++;
  stats->num_paths++;
  assert(!path->length);

  c->guide_paths[pathid].sampled++;

  while(1)
  {
#define ERR(e) {stats->err_extend[e]++; return e;}
    if(path->length && (path->length == path_len(c, pathid)))
    {
      c->guide_paths[pathid].success++;
      ERR(1); // length of guide path reached, this is not an error
    }

    if(path->length == PATHSPACE_MAX_VERTS)
      ERR(2); // max path length reached

    if(guided_num_guide_paths(c) <= 0)
      ERR(3); // no cache points

    // minimal init on next data points:
    memset(path->v+path->length, 0, sizeof(vertex_t));
    memset(path->e+path->length, 0, sizeof(edge_t));

    if(path->length == 0)
    { // start from scratch at camera:
      // init essentials on new path:
      path->tangent_frame_scrambling = 0.1f + points_rand(rt.points, common_get_threadid())*(0.9f-0.1f);
#if SPECTRAL_IS==1
      const float guide_lambda01 = (b->guide_paths[pathid].lambda - spectrum_sample_min)
        /(spectrum_sample_max - spectrum_sample_min);
      const float lambda01 = guided_mutate(guide_lambda01, pointsampler(path, s_dim_lambda));
      assert(lambda01 >= 0.0);
      assert(lambda01 <= 1.0);
      path->lambda = spectrum_sample_min + lambda01 * (spectrum_sample_max - spectrum_sample_min);
      if (path->lambda > spectrum_sample_max || path->lambda < spectrum_sample_min) ERR(4);
#else
      float lambda_pdf;
      if(path->lambda == 0.0f)
      { // if it's already set, don't sample a new one (for bdpt)
        path->lambda = spectrum_sample_lambda(pointsampler(path, s_dim_lambda), &lambda_pdf);
      }
#endif
      path->time = view_sample_time(path, pointsampler(path, s_dim_time)); // unit pdf, /= 1.0f to throughput.
      path->v[0].throughput = view_cam_sample(path);

      path->v[0].mode = s_sensor;

      // sample start point always in global exterior medium:
      shader_exterior_medium(path);
      path->v[0].throughput = 1.0f; // dummy

      // we're actually creating two vertices at this point. the first is on the camera:
      path->length ++;
      memset(path->v+path->length, 0, sizeof(vertex_t));
      memset(path->e+path->length, 0, sizeof(edge_t));
    }

    const int v = path->length;
    assert(!(path->v[v-1].flags & s_environment)); // we did test that before (err 13)
    if(!(path->v[v-1].throughput > 0.0f)) ERR(5);

    // tentatively use same bsdf reflection/transmission mode as cached path
    path->v[v-1].mode = get_vert(c, pathid, v-1)->mode;
    path->v[v].mode   = get_vert(c, pathid, v)->mode;
    path->v[v].flags  = get_vert(c, pathid, v)->flags;

    // randomly sample usage of bsdf instead of guiding path:
    float p_bsdf = guided_extend_bsdf(path, v);
    const float xi = points_rand(rt.points, common_get_threadid());

    // remember our random number offset
    path->v[v].rand_beg = path->v[v-1].rand_beg + path->v[v-1].rand_cnt;
    path->v[v].pdf = 1.0f; // everybody will just *= his sampling here.
    path->v[v].throughput = path->v[v-1].throughput;

    if(xi < p_bsdf)
    {
      // sample the bsdf at vertex v-1. also inits e[v] with volume information.
      path->v[v-1].culled_modes = ~path->v[v-1].mode;
      path->v[v].throughput *= shader_sample(path);
      path->v[v-1].culled_modes = 0;

      if(!(path->v[v].throughput > 0.0f))
        ERR(6); // return without incrementing path length

      if(path->v[v].mode & s_volume)
      {
        // sample distance from guide cache
        if(guided_sample_dist(
              c, pathid, v,
              pointsampler(path, s_dim_free_path),
              &path->e[v].dist)) ERR(7);
      }
      else path->e[v].dist = 1.0f; // dummy to generate some vertex position.
      // we'll overwrite this distance soon anyways.

      for(int k=0;k<3;k++)
        path->v[v].hit.x[k] = path->v[v-1].hit.x[k] + path->e[v].dist * path->e[v].omega[k];
    }
    else
    {
      // sample target postion:
      if(guided_sample_pos(
            c, pathid, v,
            pointsampler(path, s_dim_omega_x),
            pointsampler(path, s_dim_omega_y),
            pointsampler(path, s_dim_free_path),
            path->v[v-1].hit.x,
            path->v[v].hit.x)) ERR(8);
    }

    path->e[v].dist = FLT_MAX;
    path->v[v].hit.prim = INVALID_PRIMID;
    // mutation mode requesting volume distance without medium?
    // this may happen if mutated to the other side of a surface.
    if(path->e[v].vol.shader == -1 && (path->v[v].mode & s_volume)) ERR(9);

    int err = 0;
    if((err = path_project(path, v, s_propagate_mutate)))
    {
      stats->err_project[err]++;
      return 100+err;
    }

    // hit light source but didn't aim for it?
    if((path->v[v].mode & s_emit) != (get_vert(c, pathid, v)->mode & s_emit))
      ERR(11);

    if((path->length+1 == path_len(c, pathid)) &&
        ((path->v[v].mode & s_emit) != (get_vert(c, pathid, v)->mode & s_emit) ||
         (path->v[v].flags != get_vert(c, pathid, v)->flags)))
      ERR(15); // miss light source

    if(path->v[v].flags != get_vert(c, pathid, v)->flags) ERR(13);

    // update throughput, which here is just the full measurement contribution:
    // now we have camera connection info, too
    if(v == 1)
    {
      path->v[v].throughput = 1.0f; // cam connect throughput will end up at v[0]
    }
    else if(v > 1)
    {
      const float bsdf = shader_brdf(path, v-1); // sets the mode
      if(path->v[v-1].mode != get_vert(c, pathid, v-1)->mode) ERR(12);
      if(!(bsdf > 0.0f)) ERR(14);
      assert(path->v[v-1].mode != s_absorb);
      // TODO: numerical precision may dictate to store it without * prev throughput
      path->v[v].throughput = path->v[v-1].throughput * bsdf;
    }

    path->length++; // now a valid vertex
    path->v[v].rand_cnt = s_dim_num_extend;

    if((path->v[0].mode & s_sensor) && (path->v[v].mode & s_emit))
      path->throughput = path->v[v].throughput * lights_eval_vertex(path, v);

#undef ERR
  }
  return 0;
}

static uint32_t _lbvh_collect_surface(
    const guided_cache_t *cache,
    const lbvh_t *b,         // the bvh
    const float *x,          // origin of ray
    const float *w,          // ray direction
    uint32_t *out_list,
    const int list_alloc,
    const int plen)
{
  int list_size = 0;
  const lbvh_node_t *node = b->nodes;
  uint64_t stack[100];
  int sp = 0;
  uint64_t current;
  stack[0] = 0;
  const float iw[3] = {1.0f/w[0], 1.0f/w[1], 1.0f/w[2]};

  while(1)
  {
    float tmin[2] = {0, 0}, tmax[2] = {FLT_MAX, FLT_MAX};
    for(int k=0;k<3;k++)
    {
      for(int c=0;c<2;c++)
      {
        const float tm0 = (node->aabb[c][k]   - x[k]) * iw[k];
        const float tM0 = (node->aabb[c][k+3] - x[k]) * iw[k];
        const float tm = MIN(tm0, tM0);
        const float tM = MAX(tm0, tM0);
        tmin[c] = tm > tmin[c] ? tm : tmin[c];
        tmax[c] = tM < tmax[c] ? tM : tmax[c];
      }
    }

    if(tmin[0] <= tmax[0])
    {
      current = node->child[0];
      if(tmin[1] <= tmax[1]) stack[sp++] = node->child[1];
    }
    else if(tmin[1] <= tmax[1]) current = node->child[1];
    else
    { // pop stack
      if(sp == 0) return list_size;
      sp--;
      current = stack[sp];
    }
    while(current & (1ul<<63))
    {
      uint32_t prim = (current ^ (1ul<<63)) >> 5;
      const uint32_t cnt = current & ~-(1<<5);
      for(int i=0;i<cnt;i++)
      {
        // reject all but surface vertices
        const uint32_t pi = b->primid[prim];
        const gvert_t *vn = get_vert(cache, pi, 1);
        if(!(vn->mode & s_volume) && !(vn->flags & s_environment) &&
            (path_len(cache, pi) == plen) && (cache->guide_paths[pi].measurement_contribution > 0.0f))
        {
          if(list_size >= list_alloc) return list_size;
          out_list[list_size++] = pi;
        }
        prim++;
      }
      if(sp == 0) return list_size;
      sp--;
      current = stack[sp];
    }
    node = b->nodes + current;
  }
  assert(0);
  return 0.0;
}

// point query
static uint32_t _lbvh_pdf_volume(
    const guided_cache_t *cache,
    const lbvh_t *b,         // the bvh
    const float *x,          // query point
    uint32_t *out_list,
    const int list_alloc,
    const int plen)
{
  const lbvh_node_t *node = b->nodes;
  uint32_t stack[100];
  int sp = 0;

  int list_size = 0;

  while(1)
  {
    for(int c=0;c<2;c++)
    {
      if(node->child[c] & (1ul<<63))
      { // leaf node
        uint32_t prim = (node->child[c] ^ (1ul<<63)) >> 5;
        const uint32_t cnt = node->child[c] & ~-(1<<5);
        for(int i=0;i<cnt;i++)
        {
          // reject all but volume vertices
          const int pi = b->primid[prim];
          const gvert_t *vn = get_vert(cache, pi, 1);
          if((vn->mode & s_volume) && !(vn->flags & s_environment) &&
             (path_len(cache, pi) == plen) && (cache->guide_paths[pi].measurement_contribution > 0.0f))
          {
            if(list_size >= list_alloc) return list_size;
            out_list[list_size++] = pi;
          }
          prim++;
        }
      }
      else
      { // inner node
        if(x[0] >= node->aabb[c][0] && x[0] <= node->aabb[c][3] &&
           x[1] >= node->aabb[c][1] && x[1] <= node->aabb[c][4] &&
           x[2] >= node->aabb[c][2] && x[2] <= node->aabb[c][5])
          stack[sp++] = node->child[c];
      }
    }
    // pop stack:
    if(--sp < 0) return list_size;
    node = b->nodes + stack[sp];
  }
  return list_size;
}

static int guided_interesting_path(
    const path_t *path)
{
  return path->length > 2;
  // DEBUG caustics:
  //return path->length == 5 && (path->v[3].mode & s_transmit);
}

static int guided_interesting_gpath(
    const gpath_t *path,
    const gvert_t *verts)
{
  return path->len > 2;
  // DEBUG caustics:
  //return path->len == 5 && (verts[3].mode & s_transmit);
}

static void guided_path_to_gpath(
    path_t *path,
    gpath_t *gpath,
    gvert_t *gverts)
{
  gpath->lambda = path->lambda;
  gpath->len = path->length;
  gpath->pdf_connect = view_cam_pdf_connect(path, 0);

  for(int vi=0;vi<path->length;vi++)
  {
    gverts[vi].mode = path->v[vi].mode;
    gverts[vi].flags = path->v[vi].flags;
    if(gverts[vi].flags & s_environment)
      for(int k=0;k<3;k++) // envmaps store direction, not position
        gverts[vi].x[k] = path->e[vi].omega[k];
    else
      for(int k=0;k<3;k++)
        gverts[vi].x[k] = path->v[vi].hit.x[k];
    assert(dotproduct(gverts[vi].x, gverts[vi].x) >= 0.0);
  
    if (vi > 0)
    {
      // save scatter mode to restore later
      vertex_scattermode_t mode_store = path->v[vi-1].culled_modes;
      
      path->v[vi-1].culled_modes = ~path->v[vi-1].mode;
      gverts[vi-1].pdf_bsdf = vi > 1 ? shader_pdf(path, vi-1) : 0.0f;

      // restore scatter mode later
      path->v[vi-1].culled_modes = mode_store;
      gverts[vi-1].p_bsdf = guided_extend_bsdf(path, vi);

      gverts[vi-1].jacobian = 1.0f/
               path_lambert(path, vi-1, path->e[vi].omega);
    }
    else gverts[vi].jacobian = 1.0f;
  }
}

static int guided_pdf_collect(
    const guided_cache_t *c,
    const uint32_t vol, // second vertex in volume or surface?
    const float v0[3],   // first vertex location
    const float v1[3],   // second vertex location
    const int path_length,
    const int size,
    uint32_t *path_list)
{
  int path_cnt = 0;
  if(vol)
    path_cnt = _lbvh_pdf_volume(c,
        &c->bvh, v1,
        path_list, size, path_length);
  else
  {
    float omega[3] = {v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]}; 
    normalise(omega);
    path_cnt = _lbvh_collect_surface(c,
        &c->bvh, v0, omega,
        path_list, size, path_length);
  }
  return path_cnt;
}

// sum contributions of overlapping gaussians,
// multiply probability to choose these paths:
static double guided_pdf_accumulate(
    const guided_cache_t *c,
    const int path_cnt, 
    const uint32_t *path_list, 
    const float *path_pdf)
{
  double pdf = 0.0;
  for(int k=0;k<path_cnt;k++)
  {
#if USE_PATH_CDF==1
    const int block = path_list[k] / c->path_cdf_cnt;
    const int inner = path_list[k] - block * c->path_cdf_cnt;
    double block_p = block ?
      c->path_block_cdf[block] - c->path_block_cdf[block-1] : c->path_block_cdf[0];
    double inner_p = inner ?
      c->path_cdf[path_list[k]] - c->path_cdf[path_list[k]-1] : c->path_cdf[block * c->path_cdf_cnt];
    pdf += path_pdf[k] * inner_p * block_p;
#else
    pdf += path_pdf[k];
#endif
  }

  return pdf;
}


// evaluate pdf to sample this path with the current cache.
// only collect paths from bvh query at first point, assume these will be few points.
// then go through those w/o accel struct, thinning out as the pdf becomes 0 for higher bounces
static double guided_pdf_gpath(
    const guided_cache_t *c,
    gpath_t *gpath, 
    gvert_t *gverts,
    int *num_overlapping_gpaths)
{
  if(num_overlapping_gpaths) *num_overlapping_gpaths = 0;
  if(!guided_interesting_gpath(gpath, gverts)) return 0.0;
  // we are calling this while building the cache, before swapping over the numbers:
  if(c->num_pathverts == 0 && c->num_new_pathverts == 0) return 0.0;

  // find candidates at v[1]:
  uint32_t *path_list = c->tls[common_get_threadid()].path_list;
  float *path_pdf = c->tls[common_get_threadid()].path_pdf;
  // experiments indicate that if we could limit
  // this to say 10 we'd be no slower than ptdl:
  const uint32_t size = c->tls[common_get_threadid()].size;
  int path_cnt = guided_pdf_collect(c, gverts[1].mode & s_volume,
      gverts[0].x, gverts[1].x, gpath->len, size, path_list);
  if(path_cnt == 0) return 0.0;

  // do a quick pass of simple culling:
  for(int i=0;i<path_cnt;i++)
  {
    int vi = c->guide_paths[path_list[i]].vi+1;
    assert(c->guide_paths[path_list[i]].len == gpath->len);
    for(int j=1;j<gpath->len;j++)
    { // early out test for modes and flags
      if((c->guide_verts[vi-1].mode != gverts[j-1].mode) ||
         (c->guide_verts[vi].flags != gverts[j].flags) ||
         (c->guide_paths[path_list[i]].measurement_contribution <= 0.0f))
      { // cull path
        path_list[i--] = path_list[--path_cnt];
        break;
      }
      vi++;
    }
  }

  // for all path vertices starting at v[1]:
  // go through candidate list, thin out if pdf = 0.
  for(int vi=1;vi<gpath->len;vi++)
  {
    const float pdf_bsdf = gverts[vi-1].pdf_bsdf;
    const float jacobian = gverts[vi-1].jacobian;
    const int vol = gverts[vi].mode  & s_volume;

#if SPECTRAL_IS==1
    const float gpath_lambda01 = (gpath->lambda - spectrum_sample_min)
      /(spectrum_sample_max - spectrum_sample_min);
#endif

    for(int i=0;i<path_cnt;i++)
    {
      float pdf_extend_dist = 1.0f;
      const uint32_t pi = path_list[i];
      float p_bsdf = gverts[vi-1].p_bsdf;
      const gvert_t* vn = get_vert(c, pi, vi);
      const gvert_t* vp = get_vert(c, pi, vi-1);
      if(vi == 1)
      {
#if SPECTRAL_IS==1
        const float guide_lambda01 = (c->guide_paths[pi].lambda - spectrum_sample_min)
          /(spectrum_sample_max - spectrum_sample_min);
        const float lambda_pdf = guided_mutate_pdf(gpath_lambda01, guide_lambda01);
        //this should be
        //lambda_pdf *= spectrum_lambda_pdf to upscale the mutate_pdf from [0,1] to the wavelength range
        //lambda_pdf /= spectrum_lambda_pdf because the unguided sampling ignores the uniform wavelength pdf
        //which cancels out to 1
        path_pdf[i] = lambda_pdf; // init
#else
        path_pdf[i] = 1.f; // init
#endif
      }

      // pdf for path extension via bsdf
      if(p_bsdf > 0.0 && vol)
      {
        pdf_extend_dist = guided_volume_gauss_distance_pdf(
            c, vn->sigma_dist, gverts[vi-1].x, 
            gverts[vi].x, vp->x, vn->x);
      }

      float pdf_pos = 0;
      if (p_bsdf < 1)
        pdf_pos = jacobian * guided_pos_pdf(c, vp, vn, gverts[vi-1].x, gverts[vi].x);
      path_pdf[i] *= (1.f-p_bsdf) * pdf_pos
        + p_bsdf * pdf_bsdf * pdf_extend_dist;

      if(!(path_pdf[i] == path_pdf[i]) || isinf(path_pdf[i]))
        path_pdf[i] = 0.0f;

      if(path_pdf[i] == 0.0f)
      { // remove path in case pdf = 0
        path_pdf[i]    = path_pdf[--path_cnt];
        path_list[i--] = path_list[path_cnt];
      }
    }
  }

  if(path_cnt == 0) return 0.0;
  double pdf = guided_pdf_accumulate(c, path_cnt, path_list, path_pdf);
#if USE_PATH_CDF!=1
  pdf /= guided_num_guide_paths(c);
#endif
  assert(pdf == pdf);
  // multiply aperture pdf from camera for v[0]
  pdf *= gpath->pdf_connect;

  assert(pdf == pdf);
  if(num_overlapping_gpaths) *num_overlapping_gpaths = path_cnt;
  return pdf;
}

// evaluate pdf to sample this path with the current cache.
// only collect paths from bvh query at first point, assume these will be few points.
// then go through those w/o accel struct, thinning out as the pdf becomes 0 for higher bounces
double guided_pdf_dwp(
    const guided_cache_t *c,
    path_t *path)
{
  if(!guided_interesting_path(path)) return 0.0;
  if(c->num_pathverts == 0) return 0.0;

  uint32_t *path_list = c->tls[common_get_threadid()].path_list;
  float *path_pdf = c->tls[common_get_threadid()].path_pdf;
  const uint32_t size = c->tls[common_get_threadid()].size;
  // find candidates at v[1]:
  int path_cnt = guided_pdf_collect(c, path->v[1].mode & s_volume,
      path->v[0].hit.x, path->v[1].hit.x, path->length, size, path_list);
  if(path_cnt == 0) return 0.0;

  // do a quick pass of simple culling:
  for(int i=0;i<path_cnt;i++)
  {
    assert(c->guide_paths[path_list[i]].len == path->length);
    int vi = c->guide_paths[path_list[i]].vi+1;
    if (c->guide_paths[path_list[i]].len != path->length)
    {
      path_list[i--] = path_list[--path_cnt];
      continue;
    }

    for(int j=1;j<path->length;j++)
    { // early out test for modes and flags
      if((c->guide_verts[vi-1].mode != path->v[j-1].mode) ||
         (c->guide_verts[vi].flags != path->v[j].flags) ||
         (c->guide_paths[path_list[i]].measurement_contribution <= 0.0f))
      { // cull path
        path_list[i--] = path_list[--path_cnt];
        break;
      }
      vi++;
    }
  }
  
  // for all path vertices starting at v[1]:
  // go through candidate list, thin out if pdf = 0.
  for(int vi=1;vi<path->length;vi++)
  {
    // save scatter mode to restore later
    vertex_scattermode_t mode_store = path->v[vi-1].culled_modes;

    path->v[vi-1].culled_modes = ~path->v[vi-1].mode;
    const float pdf_bsdf = vi > 1 ? shader_pdf(path, vi-1) : 0.0f;
    const int vol = path->v[vi].mode & s_volume;

    // restore scatter mode
    path->v[vi-1].culled_modes = mode_store;

    // used to convert pos_pdf from solid angle to projected solid angle
    const float jacobian = 1.0f/path_lambert(path, vi-1, path->e[vi].omega);

#if SPECTRAL_IS==1
    const float path_lambda01 = (path->lambda - spectrum_sample_min)
      /(spectrum_sample_max - spectrum_sample_min);
#endif

    for(int i=0;i<path_cnt;i++)
    {
      float pdf_extend_dist = 1.0f;
      const uint32_t pi = path_list[i];
      float p_bsdf = guided_extend_bsdf(path, vi);
      const gvert_t* vn = get_vert(c, pi, vi);
      const gvert_t* vp = get_vert(c, pi, vi-1);
      if(vi == 1)
      {
#if SPECTRAL_IS==1
        const float guide_lambda01 = (c->guide_paths[pi].lambda - spectrum_sample_min)
          /(spectrum_sample_max - spectrum_sample_min);
        const float lambda_pdf = guided_mutate_pdf(path_lambda01, guide_lambda01);
        path_pdf[i] = lambda_pdf; // init
#else
        path_pdf[i] = 1.f; // init
#endif
      }

      // pdf for path extension via bsdf
      if(p_bsdf > 0.0 && vol)
      {
        pdf_extend_dist = guided_volume_gauss_distance_pdf(
            c, vn->sigma_dist, path->v[vi-1].hit.x, 
            path->v[vi].hit.x, vp->x, vn->x);
      }

      float pdf_pos = 0;
      if (p_bsdf < 1) // if envmap pass direction instead:
        pdf_pos = jacobian * guided_pos_pdf(c, vp, vn, path->v[vi-1].hit.x,
            (path->v[vi].flags & s_environment) ?
            path->e[vi].omega :
            path->v[vi].hit.x);
      path_pdf[i] *= ((1.f-p_bsdf) * pdf_pos
                   + p_bsdf * pdf_bsdf * pdf_extend_dist);

      // XXX: FIXME path_pdf[i] can be inf
      if(!(path_pdf[i] == path_pdf[i]) || isinf(path_pdf[i]))
        path_pdf[i] = 0.0f;

      if(path_pdf[i] == 0.0f)
      { // remove path in case pdf = 0
        path_pdf[i]    = path_pdf[--path_cnt];
        path_list[i--] = path_list[path_cnt];
      }
    }
  }

  if(path_cnt == 0) return 0.0;
  // sum contributions of overlapping gaussians,
  // multiply probability to choose these paths:
  double pdf = guided_pdf_accumulate(c, path_cnt, path_list, path_pdf);

#if USE_PATH_CDF!=1
  pdf /= guided_num_guide_paths(c);
#endif
  assert(pdf == pdf);
  // multiply aperture pdf from camera for v[0]
  pdf *= view_cam_pdf_connect(path, 0);

  assert(pdf == pdf);
  return pdf;
}

// ================================================================================
//  record a given path as guide path
// ================================================================================
int guided_record_path(
    guided_cache_t *c, 
    const path_t *path, 
    const uint32_t gpath_id,
    const double measurement_contribution, 
    const double guided_pdf,
    const double uniform_pdf,
    const float uniform_ratio,
    const sampling_type_t st)
{
  const double pdf = uniform_ratio * uniform_pdf + (1-uniform_ratio) * guided_pdf;
  const double throughput = measurement_contribution / pdf;

  if(!(throughput > 1e-10f)) return 0;
  if(!(throughput < FLT_MAX)) return 0;

  dbor_splat(c->dbor, path->sensor.pixel_i, path->sensor.pixel_j, throughput);
  
#if USE_DBOR==1
  {
    const double trust = dbor_trust(c->dbor, path->sensor.pixel_i, path->sensor.pixel_j, throughput) / (rt.frames+1);

    if(trust > TRUST_RECORD_THR) 
    {
      if (st == st_guided)
      {
        const gpath_t* gpath = c->guide_paths+gpath_id;
        // only reject if it has high dimensional neighbours
        if (gpath->nb_cnt > 9) return 0;
      }
      else return 0;
    }
  }
#endif

  const float contrib = throughput;
  if(!(contrib < FLT_MAX)) return 0;
  assert(contrib > 0);
  
  guided_stats_t *stats = &c->tls[common_get_threadid()].stats;
  if(guided_interesting_path(path))
  {
    uint64_t pathid = (uint64_t)(-1);
    heap_t* heap = c->tls[common_get_threadid()].heap;
    if (heap_full(heap)) 
    {
      if (1.f/contrib > heap_max_val(heap))
        return 0;

      float val;
      heap_remove(heap, &pathid, &val);
    }

    // use 64-bit atomics to increment path and vertex buffer at the same time:
    uint64_t inc = pathid == (uint64_t)(-1) ? (1ul<<32) | path->length : path->length;
    uint64_t pv = __sync_fetch_and_add(&c->num_new_pathverts, inc);
    uint32_t path_idx = pathid == (uint64_t)(-1) ? pv >> 32 : pathid;
    uint32_t vert_idx = pv & 0xffffffff;
    if(path_idx >= c->num_paths || vert_idx + path->length >= c->num_verts)
    {
      __sync_fetch_and_add(&c->num_new_pathverts, -inc);
      return 1; // path buffer full?
    }
    
    heap_insert(heap, path_idx, 1.f/contrib);

    gpath_t *p = c->guide_paths + path_idx;
    guided_path_to_gpath((path_t*)path, p, c->guide_verts + vert_idx);
    
    p->vi = vert_idx;
    p->measurement_contribution = measurement_contribution;
    p->guided_pdf = guided_pdf;
    p->uniform_pdf = uniform_pdf;
    p->lambda = path->lambda;
    p->len = path->length;
    p->pixel_x = path->sensor.pixel_i;
    p->pixel_y = path->sensor.pixel_j;
    p->sampled = 0;
    p->success = 0;
    p->pixel_step = RAY_DIFF_MAX;
    p->age = 0;

    gpath_set_sampling_type(p, st);
    assert(gpath_get_sampling_type(p) == st);
    if (st == st_uniform) stats->num_uniform_paths_recorded++;
    else                  stats->num_guiding_paths_recorded++;

    for(int j=0;j<path->length;j++)
    {
      c->guide_verts[vert_idx].mode = path->v[j].mode;
      c->guide_verts[vert_idx].flags = path->v[j].flags;

      if (path->v[j].mode & s_volume)
        c->guide_verts[vert_idx].roughness = 1.f;
      else
        c->guide_verts[vert_idx].roughness = path->v[j].shading.roughness;

      if(c->guide_verts[vert_idx].flags & s_environment)
        for(int k=0;k<3;k++) // envmaps store direction, not position
          c->guide_verts[vert_idx].x[k] = path->e[j].omega[k];
      else
        for(int k=0;k<3;k++)
          c->guide_verts[vert_idx].x[k] = path->v[j].hit.x[k];
      assert(dotproduct(
            c->guide_verts[vert_idx].x,
            c->guide_verts[vert_idx].x) >= 0.0);
      vert_idx++;
    }
  }
  return 0;
}

void guided_clear(guided_cache_t *c)
{
  c->learn_iteration = 0;
  c->num_pathverts = 0;
  c->num_new_pathverts = 0;
  c->gauss = 2; // fake gaussian by default and for learning
  dbor_clear(c->dbor);
}

void guided_collect_stats(guided_cache_t *s, FILE* f)
{
  static const char *err_extend_str[16] = {
    "no error",
    "length of guide path reached",
    "max path vertices reached",
    "no cache points",
    "spectral sampling failed",
    "no previous throughput",
    "bsdf sampling failed",
    "distance sampling failed",
    "position sampling failed",
    "volume stack inconsistent",
    "",
    "premature light source",
    "mode inconsistent",
    "flags inconsistent",
    "bsdf zero",
    "miss light source",
  };

  guided_stats_t stats = {0};
  for(int t=0;t<rt.num_threads;t++)
  {
    stats.num_samples += s->tls[t].stats.num_samples;
    stats.num_paths += s->tls[t].stats.num_paths;
    for(int k=0;k<16;k++) stats.err_project[k] += s->tls[t].stats.err_project[k];
    for(int k=0;k<16;k++) stats.err_extend[k] += s->tls[t].stats.err_extend[k];
    stats.num_uniform_paths_recorded += s->tls[t].stats.num_uniform_paths_recorded;
    stats.num_uniform_paths_in_cache += s->tls[t].stats.num_uniform_paths_in_cache;
    stats.num_guiding_paths_recorded += s->tls[t].stats.num_guiding_paths_recorded;
    stats.num_guiding_paths_in_cache += s->tls[t].stats.num_guiding_paths_in_cache;
    memset(&s->tls[t].stats, 0, sizeof(guided_stats_t));
  }
  uint64_t fails = 0;
  for(int k=1;k<16;k++) fails += stats.err_project[k];
  for(int k=2;k<16;k++) fails += stats.err_extend[k];

  fprintf(f, "guided stats:\n");
  fprintf(f, "  total successful extend calls %.2f%% (%lu/%lu)\n",
      (stats.num_samples-fails)*100.0/stats.num_samples,
      (stats.num_samples-fails), stats.num_samples);
  // extend error 01 is no error.
  for(int k=2;k<16;k++)
    if(stats.err_extend[k])
      fprintf(f, "  extend fail  %02d : %.2f%% %lu (%s)\n", k, stats.err_extend[k]*100.0/stats.num_paths, stats.err_extend[k], err_extend_str[k]);
  for(int k=0;k<16;k++)
    if(stats.err_project[k])
      fprintf(f, "  project fail %02d : %.2f%% %lu\n", k, stats.err_project[k]*100.0/stats.num_paths, stats.err_project[k]);

  const uint64_t num_paths_in_cache = stats.num_guiding_paths_in_cache + stats.num_uniform_paths_in_cache;
  fprintf(f, "  num paths in cache: %lu\n", num_paths_in_cache);
  fprintf(f, "  ratio of guide paths in cache (rest uniform): %.2f%%\n", stats.num_guiding_paths_in_cache*100.0/num_paths_in_cache);
  fprintf(f, "    uniform path survival: %.2f%% (%lu/%lu)\n", stats.num_uniform_paths_in_cache*100.0/stats.num_uniform_paths_recorded,
      stats.num_uniform_paths_in_cache, stats.num_uniform_paths_recorded);
  fprintf(f, "    guiding path survival: %.2f%% (%lu/%lu)\n", stats.num_guiding_paths_in_cache*100.0/stats.num_guiding_paths_recorded,
      stats.num_guiding_paths_in_cache, stats.num_guiding_paths_recorded);
}

void guided_export_cache_info(guided_cache_t *c)
{
  dbor_export(c->dbor, "dbor", 1);
  guided_debug_path_sampling(c, 100);

  const int old_num_verts = c->num_pathverts & 0xffffffff;
  FILE *f = fopen("cache.obj", "w");
  for(int v=0;v<old_num_verts;v++)
  {
    if(c->guide_verts[v].flags & s_environment)
    { // envmap vertices encode direction, not position
      fprintf(f,  "v %g %g %g\n",
          c->guide_verts[v-1].x[0]+c->guide_verts[v].x[0],
          c->guide_verts[v-1].x[1]+c->guide_verts[v].x[1],
          c->guide_verts[v-1].x[2]+c->guide_verts[v].x[2]);
    }
    else
    {
      fprintf(f,  "v %g %g %g\n",
          c->guide_verts[v].x[0],
          c->guide_verts[v].x[1],
          c->guide_verts[v].x[2]);
    }
  }

  fprintf(f, "o cache-vertices\n");
  for(int v=0;v<old_num_verts;v++)
    fprintf(f,  "f %d\n", v+1);

  fclose(f);
}

static inline void raydiff(
    const guided_cache_t *c,
    const int pi,
    const float pixel_step,
    float *min,
    float *max)
{
  path_t tmp;
  tmp.v[0].mode = s_sensor;
  tmp.sensor.aperture_set = 1;
  tmp.sensor.aperture_x = 0.0f;
  tmp.sensor.aperture_y = 0.0f;
  for(int k=0;k<3;k++) tmp.v[1].hit.x[k] = get_vert(c, pi, 1)->x[k];
  tmp.time = 0.0f;
  // XXX TODO: fix this for stereo
  tmp.sensor.camid = 0;
  tmp.lambda = 550.0f;
  tmp.index = 0;
  view_cam_connect(&tmp);
  tmp.sensor.pixel_i = CLAMP(tmp.sensor.pixel_i, 0.f, ((float)view_width())-1e-4f);
  tmp.sensor.pixel_j = CLAMP(tmp.sensor.pixel_j, 0.f, ((float)view_height())-1e-4f);
  for(int k=0;k<3;k++) tmp.v[1].hit.n[k] = -tmp.e[1].omega[k];
  float rd_i[3], rd_j[3];
  if(!raydifferentials_v1(&tmp, 1, 1, rd_i, rd_j))
    *min = MAX(dotproduct(rd_i, rd_i), dotproduct(rd_j, rd_j));
  if(!raydifferentials_v1(&tmp, pixel_step, pixel_step, rd_i, rd_j))
    *max = MAX(dotproduct(rd_i, rd_i), dotproduct(rd_j, rd_j));
}

static inline void init_vertex(
    gvert_t *vn,
    double *Sc,      // 3x3 conditional covariance block
    double *mean,    // 6D mean of gaussian
    double *S12S22i, // 3x3 block to be used when computing the conditional mean
    double gauss_min_dist)
{
  const int vn_vol = vn->mode  & s_volume;
  const int vn_env = vn->flags & s_environment;
  double E[9] = {0}, Et[9] = {0};
  double w[3] = {0};
  if (vn_vol)
  {
    getEVDSymmetric3x3(w, Et, Sc);
    mat3d_transpose(Et, E);
  }
  else
  {
    double w2d[2] = {0};
    double Et2d[4] = {0};
    double Sc2d[4] = { Sc[0], Sc[1], Sc[3], Sc[4] };
    getEVDSymmetric2x2(w2d, Et2d, Sc2d);
    w[0] = w2d[0];
    w[1] = w2d[1];
    w[2] = 1;
    Et[0] = Et2d[0];
    Et[1] = Et2d[1];
    Et[3] = Et2d[2];
    Et[4] = Et2d[3];
    Et[8] = 1;
    mat3d_transpose_sub2(Et, E);
  }

  // clamp eigenvalues
  if(vn_env)
    for (int k=0;k<3;++k)
      w[k] = MAX(0.00001f, w[k]);
  else
    for (int k=0;k<3;++k)
      w[k] = MAX(gauss_min_dist*gauss_min_dist, w[k]);

  for(int k=0;k<9;++k)
    vn->E[k] = E[k];
  for (int k=0;k<3;++k)
    vn->w[k] = sqrtf(w[k]);
  for(int k=0;k<9;k++)
    vn->S12S22i[k] = S12S22i[k];
  for (int k=0;k<6;k++)
    vn->mean[k] = mean[k];

  if (vn_vol) 
  {
    // use regularised version here, too:
    double Sr[9], Tmp[9];
    for(int j=0;j<3;j++)
    for(int i=0;i<3;i++) Tmp[3*j+i] = w[j] * Et[3*j+i];
    mat3d_mul(E, Tmp, Sr);
    const double S11[4] = {Sr[0], Sr[1], Sr[3], Sr[4]};
    const double S21[2] = {Sr[6], Sr[7]};
    double S11i[4] = {0}, tmp[2] = {0};
    mat2d_invert(S11, S11i);
    mat2d_mulv(S11i, S21, tmp);
    const double sigma = sqrt(MAX(1e-5, Sr[8] - tmp[0]*S21[0] - tmp[1]*S21[1]));
    vn->sigma_dist = sigma;
  }
}

static void *par_compute_pdf(void *arg)
{
  guided_cache_t *c = arg;
  while(1)
  {
    uint64_t p = __sync_fetch_and_add(&c->job, 1);
    if(p >= c->num_jobs) return 0;
    gpath_t* gpath = c->guide_paths+p;
    gvert_t* gverts = c->guide_verts+gpath->vi;
    gpath->guided_pdf = guided_pdf_gpath(c, gpath, gverts, 0);
  }
}

static void *par_compute_raydiff(void *arg)
{
  guided_cache_t *c = arg;
  while(1)
  {
    uint64_t p = __sync_fetch_and_add(&c->job, 1);
    if(p >= c->num_jobs) return 0;

    // init with zero bounds around points:
    gvert_t* vn = get_vert(c, p, 1);
    const gvert_t* vp = get_vert(c, p, 0);
    const int vn_vol = vn->mode & s_volume;
    // compute onb for vertex vn
    float u[3], v[3], n[3];
    for(int k=0;k<3;k++) n[k] = vn->x[k] - vp->x[k];
    normalise(n);
    get_onb(n, u, v);

    gpath_t* gpath = &c->guide_paths[p];
    float ps_max = MAX(1.f, ((float)RAY_DIFF_MAX)/(c->learn_iteration) + RAY_DIFF_MIN);
    float ps_min = MAX(1.f, ((float)RAY_DIFF_MAX)/(c->learn_iteration) + 1);
#if RAY_DIFF_ADAPTIVE==1
    gpath->pixel_step = vn->roughness * ps_max + (1.f - vn->roughness) * ps_min;
#else
    gpath->pixel_step = ps_max;
#endif

    float gauss_dist2 = 1, gauss_min_dist = 0;
    raydiff(c, p, gpath->pixel_step, &gauss_min_dist, &gauss_dist2);
    gauss_min_dist = sqrtf(gauss_min_dist);

    double Sc[9] = {0}, mean[6] = {0}, S12S22i[9] = {0};
    mat3d_set_identity(Sc);
    for(int k=0;k<9;k++)
      Sc[k] *= gauss_dist2;
    mean[0] = dotproduct(u, vn->x);
    mean[1] = dotproduct(v, vn->x);
    mean[2] = vn_vol ? dotproduct(n, vn->x) : 0.f;

    init_vertex(vn, Sc, mean, S12S22i, gauss_min_dist);

    // grow bounding boxes by gaussian support:
    // pretend to sample really far gaussian distances and bound that.
    float Bi[9] = {0};
    float D[9] = {0};
    D[0] = vn->w[0];
    D[4] = vn->w[1];
    D[8] = vn->w[2];
    mat3_mul(vn->E, D, Bi);
    for(int k=0;k<3;k++)
      c->bvh.aabb[6*p+k] = c->bvh.aabb[6*p+k+3] = vn->x[k];
    for(int i=0;i<(vn_vol ? 8 : 4);i++)
    {
      const float cutoff = (c->gauss == 3) ? CUTOFF : 4.0f;
      float gt[3], g[] = { (i&1)?cutoff:-cutoff, (i&2)?cutoff:-cutoff, vn_vol ? (i&4)?cutoff:-cutoff : 0.0};
      mat3_mulv(Bi, g, gt);
      float corner[3];
      for(int k=0;k<3;k++) 
      {
        if (vn_vol)
        {
          corner[k] =
            + (vn->mean[0] + gt[0]) * u[k]
            + (vn->mean[1] + gt[1]) * v[k]
            + (vn->mean[2] + gt[2]) * n[k];
        }
        else
        {
          corner[k] = dotproduct(vn->x, n) * n[k]
            + (vn->mean[0] + gt[0]) * u[k]
            + (vn->mean[1] + gt[1]) * v[k];
        }
      }
      for(int k=0;k<3;k++) c->bvh.aabb[6*p+k] = MIN(c->bvh.aabb[6*p+k], corner[k]);
      for(int k=3;k<6;k++) c->bvh.aabb[6*p+k] = MAX(c->bvh.aabb[6*p+k], corner[k-3]);
    }
  }
}

static void *par_compute_neighbours(void *arg)
{
  guided_cache_t *c = arg;
  while(1)
  {
    uint64_t pi = __sync_fetch_and_add(&c->job, 1);
    if(pi >= c->num_jobs) return 0;

    gpath_t* gpath_i = c->guide_paths+pi;
    const gvert_t* gverts_i = c->guide_verts+gpath_i->vi;

    // find candidates by searching bvh:
    uint32_t *path_list = c->tls[common_get_threadid()].path_list;
    const uint32_t size = c->tls[common_get_threadid()].size;
    int path_cnt = guided_pdf_collect(c, gverts_i[1].mode & s_volume,
        gverts_i[0].x, gverts_i[1].x, gpath_i->len, size, path_list);
    for(int pjj=0;pjj<path_cnt;pjj++)
    {
      const int pj = path_list[pjj];
      double dist_path = 0;
      if(c->guide_paths[pi].len != c->guide_paths[pj].len)
        continue;

      int valid_path = 1;
      const gpath_t* gpath_j = c->guide_paths+pj;
      const gvert_t* gverts_j = c->guide_verts+gpath_j->vi;
      for(int vi=1;vi<gpath_i->len;vi++)
      { // early out test for modes and flags
        if((gverts_i[vi-1].mode != gverts_j[vi-1].mode) ||
            (gverts_i[vi].flags != gverts_j[vi].flags) ||
            (gpath_j->measurement_contribution <= 0.0f))
        { // cull path
          valid_path = 0;
          break;
        }

        // need to discard vertices [vi==1] if outside of ray diff gaussian!
        if(vi == 1)
        {
          double v1_gauss =
            guided_pos_pdf(c, gverts_i + vi-1, gverts_i + vi, gverts_j[vi-1].x, gverts_j[vi].x);
          if(v1_gauss <= 0.0)
          {
            valid_path = 0;
            break;
          }
        }

        float dx[3] = { 
          gverts_i[vi].x[0] - gverts_j[vi].x[0], 
          gverts_i[vi].x[1] - gverts_j[vi].x[1],
          gverts_i[vi].x[2] - gverts_j[vi].x[2] };
        dist_path += dotproduct(dx, dx);
      }
      if (!valid_path)
        continue;

      float dist_max = dist_path;
      int dist_max_idx = -1;
      for (int k=0;k<NUM_NB;++k)
      {
        if (c->nearest_paths_dists[pi*NUM_NB+k] > dist_max)
        {
          dist_max = c->nearest_paths_dists[pi*NUM_NB+k];
          dist_max_idx = k;
        }
      }
      if (dist_max_idx >= 0)
      {
        c->nearest_paths[pi*NUM_NB+dist_max_idx] = pj;
        c->nearest_paths_dists[pi*NUM_NB+dist_max_idx] = dist_path;
      }
    }

    gpath_i->nb_cnt = 0;
    for (int k=0;k<NUM_NB;++k)
    {
      if (c->nearest_paths[pi*NUM_NB+k] != -1)
        gpath_i->nb_cnt++;
    }
  }
}

static void *par_compute_gaussians(void *arg)
{
  guided_cache_t *c = arg;
  while(1)
  {
    uint64_t pi = __sync_fetch_and_add(&c->job, 1);
    if(pi >= c->num_jobs) return 0;
      
    int nb_cnt = c->guide_paths[pi].nb_cnt;
    float gauss_min_dist = 0;
    float gauss_dist2 = 1;
    {
      // min for eigenvalue clamping is always 1 pixel step
      float tmp;
      raydiff(c, pi, 1, &gauss_min_dist, &tmp);
      gauss_min_dist = sqrtf(gauss_min_dist);
    
      // for lonely samples use the age based step_size with
      // progressive shrinking that the change for sampling sucess
      // increases over time
      int path_age = c->guide_paths[pi].age + 1;
      float ps = MAX(1.f, ((float)RAY_DIFF_MAX/(path_age+1)));
      raydiff(c, pi, ps, &tmp, &gauss_dist2);
    }

    for(int vi=2;vi<path_len(c, pi);vi++)
    {
      // compute bounds for all paths
      gvert_t* vn = get_vert(c, pi, vi);

      // compute onb for vertex vn
      float u[3], v[3], n[3];
      const gvert_t* vp = get_vert(c, pi, vi-1);
      if(vn->flags & s_environment)
        for(int k=0;k<3;k++) n[k] = vn->x[k]; // envmaps store direction, not position
      else
        for(int k=0;k<3;k++) n[k] = vn->x[k] - vp->x[k];
      float vn_dist = sqrtf(dotproduct(n, n));
      for(int k=0;k<3;k++) n[k] /= vn_dist;
      get_onb(n, u, v);

      const int vn_vol = vn->mode & s_volume;
      const int vp_vol = vp->mode & s_volume;

      double mean[6] = {0};    // mean of the 3d or 6d gaussian
      double Sc[9] = {0};      // covariance matrix for which the eigenvalue decomposition should be computed
      double S12S22i[9] = {0}; // conditionalzation matrix if available

      if (nb_cnt < 10 || SIMPLE_GAUSSIAN_SAMPLING)
      {
        mean[0] = dotproduct(u, vn->x);
        mean[1] = dotproduct(v, vn->x);
        mean[2] = vn_vol ? dotproduct(n, vn->x) : 0.f;
        mat3d_set_identity(Sc);
        for(int k=0;k<9;k++)
          Sc[k] *= gauss_dist2;
      }
      else  // regular case compute 6x6 covariance matrix
      {
        // compute mean in uvn space
        for (int k=0; k<NUM_NB; k++)
        {
          const int pj = c->nearest_paths[pi*NUM_NB+k];
          if (pj == -1)
            continue;

          gvert_t* vpn = get_vert(c, pj, vi-1);
          gvert_t* vnn = get_vert(c, pj, vi);

          mean[0] += dotproduct(u, vnn->x);
          mean[1] += dotproduct(v, vnn->x);
          mean[2] += vn_vol ? dotproduct(n, vnn->x) : 0.f;
          mean[3] += dotproduct(u, vpn->x);
          mean[4] += dotproduct(v, vpn->x);
          mean[5] += vp_vol ? dotproduct(n, vpn->x) : 0.f;
        }
        assert(nb_cnt > 2);
        for(int k=0;k<6;++k)
          mean[k] /= nb_cnt;

        // we compute the following 4d-6d covariance matrix in blocks:
        //     | S11 S12 |
        // S = | S21 S22 | 
        double S11[9] = {0}, S12[9] = {0}, S21[9] = {0}, S22[9] = {0};
        // compute covariance matrix
        for (int k=0; k<NUM_NB; k++)
        {
          const int pj = c->nearest_paths[pi*NUM_NB+k];
          if (pj == -1)
            continue;

          gvert_t* vpn = get_vert(c, pj, vi-1);
          gvert_t* vnn = get_vert(c, pj, vi);

          for(int j=0;j<6;j++)
          {
            if(!vn_vol && j==2) continue;
            if(!vp_vol && j==5) continue;
            double dy = 0;
            if(j == 0) dy = dotproduct(u, vnn->x);
            if(j == 1) dy = dotproduct(v, vnn->x);
            if(j == 2) dy = dotproduct(n, vnn->x);
            if(j == 3) dy = dotproduct(u, vpn->x);
            if(j == 4) dy = dotproduct(v, vpn->x);
            if(j == 5) dy = dotproduct(n, vpn->x);
            for(int i=j;i<6;i++)
            {
              if(!vn_vol && i==2) continue;
              if(!vp_vol && i==5) continue;
              double dx = 0;
              if(i == 0) dx = dotproduct(u, vnn->x);
              if(i == 1) dx = dotproduct(v, vnn->x);
              if(i == 2) dx = dotproduct(n, vnn->x);
              if(i == 3) dx = dotproduct(u, vpn->x);
              if(i == 4) dx = dotproduct(v, vpn->x);
              if(i == 5) dx = dotproduct(n, vpn->x);
              assert(dx == dx);
              assert(dy == dy);

              const double cov = (dx-mean[i])*(dy-mean[j])/ (nb_cnt - 1);
              assert(cov == cov);
              if     (i <  3 && j <  3) S11[3* j   +i  ] = S11[3* i   +j  ] += cov;
              else if(i >= 3 && j >= 3) S22[3*(j-3)+i-3] = S22[3*(i-3)+j-3] += cov;
              else if(j <  3 && i >= 3) S21[3*(i-3)+j  ] = S12[3* j   +i-3] += cov;
            }
          }
        }
        // check for nans
        for(int k=0;k<9;k++) assert(S11[k] == S11[k]);
        for(int k=0;k<9;k++) assert(S12[k] == S12[k]);
        for(int k=0;k<9;k++) assert(S21[k] == S21[k]);
        for(int k=0;k<9;k++) assert(S22[k] == S22[k]);

        double tmp[9];
        double dettmp = 0;
        if (vp_vol) dettmp = mat3d_invert(S22, tmp);
        else        
        {
          dettmp = mat3d_invert_sub2(S22, tmp);
          tmp[8]=0.f;
        }
        if(dettmp <= 0.0)
        {
          c->guide_paths[pi].measurement_contribution = 0.f;
          c->guide_paths[pi].weight = 0.f;
          continue;
        }

        mat3d_mul(S12, tmp, S12S22i);
        mat3d_mul(S12S22i, S21, tmp);
        mat3d_sub(S11, tmp, Sc);
      }
      // now that we have a covariance matrix for the different settings, init
      // the vertex:
      init_vertex(vn, Sc, mean, S12S22i, gauss_min_dist);
    }
  }
}

void guided_build_cdf(guided_cache_t *c, 
    const int num_nb,
    const float uniform_ratio,
    const render_mode_t rm)
{
  if (rm == rm_learn) {
    c->learn_iteration++;
  }
  guided_stats_t *stats = &c->tls[common_get_threadid()].stats;
  if(c->num_new_pathverts == c->num_pathverts && rm == rm_learn) return; // nothing to do

  if(c->num_new_pathverts == 0) return;
  
  const int num_paths = c->num_new_pathverts>>32;
  const int old_path_cnt = guided_num_guide_paths(c);

  /*
  // only rebuild if more then 50% new paths
  uint32_t diff = num_paths - old_path_cnt;
  if (diff < old_path_cnt / 4 && rt.frames > 8 && rm == rm_learn)
    return;
  */

  // have at least 75% of heap full
  uint32_t heap_size = c->tls[common_get_threadid()].heap->size;
  if ((old_path_cnt > 1000) && (num_paths - old_path_cnt < rt.num_threads*heap_size/2) && rm == rm_learn)
    return; // don't recompute cache. too few new paths.
  
  // switch to constant kernel. do this before re-evaluating all pdf for the
  // old guide paths.
  if(rm == rm_render) c->gauss = 3; // constant
  else c->gauss = 2; // reset to fake gauss, just to be sure.

#if RENDER_WITH_GAUSS==1
  c->gauss = 2;
#endif

  c->render_mode = rm;

  double time_build_cdf = common_time_wallclock();
  double time_build_cdf_recompute = common_time_wallclock();

#if 0 // debug: dump histogram of path lengths in cache:
  int cache_pl_hist[6] = {0};
  for (int p = 0; p < num_paths; ++p)
  {
    gpath_t* gpath = c->guide_paths+p;
    int bin = CLAMP(gpath->len-3, 0, 5);
    cache_pl_hist[bin]++;
  }
  printf("pl < 4: %g%%\n", (100.f*cache_pl_hist[0])/num_paths);
  for (int i = 1; i < 5; ++i)
    printf("pl  %d: %g%%\n", i+3, (100.f*cache_pl_hist[i])/num_paths);
  printf("pl > %d: %g%%\n", 4+3, (100.f*cache_pl_hist[5])/num_paths);
#endif

  
  threads_t *t = rt.threads;
#if USE_PATH_CDF==1
  // compute new guide cdf for old guide paths
  c->job = 0;
  c->num_jobs = old_path_cnt;
  for(int k=0;k<rt.num_threads;k++)
    pthread_pool_task_init(t->task + k, &t->pool, par_compute_pdf, c);
  pthread_pool_wait(&t->pool);
#endif
  
#if 0 // debug: dump count of newly discovered paths
  int cnt_learn_paths = 0;
  for (int p=old_path_cnt; p<num_paths; ++p)
  {
    gpath_t* gpath = c->guide_paths+p;
    if (!(gpath->guided_pdf > 0))
      cnt_learn_paths++;
  }
  fprintf(stderr, "learned %d new independent paths out of %d new total\n", cnt_learn_paths, num_paths - old_path_cnt);
#endif
  
  time_build_cdf_recompute = common_time_wallclock() 
                           - time_build_cdf_recompute;
  


  //===============================================================================================================
  // build bvh around first vertices in scene v[1], according to uniform ray diff gaussian
  // assume that no 2-vert paths on the envmap are in here, i.e. v[1] is always a valid surface or volume vertex
  if(c->bvh.num_alloced_prims < num_paths)
  {
    lbvh_cleanup(&c->bvh);
    lbvh_init(&c->bvh, num_paths);
  }
  else c->bvh.num_prims = num_paths;
  c->job = 0;
  c->num_jobs = num_paths;
  for(int k=0;k<rt.num_threads;k++)
    pthread_pool_task_init(t->task + k, &t->pool, par_compute_raydiff, c);
  pthread_pool_wait(&t->pool);
  lbvh_build(&c->bvh);


  double time_build_cdf_gaussians = common_time_wallclock();

  //===============================================================================================================
  // compute nearest neighbour paths
  if (1) //rm == rm_render)
  {
    //printf("compute nearest paths\n");
    for(int pi=0;pi<num_paths;pi++)
    for(int k=0;k<NUM_NB;++k)
    {
      c->nearest_paths[pi*NUM_NB+k] = -1;
      c->nearest_paths_dists[pi*NUM_NB+k] = FLT_MAX;
    }
    c->job = 0;
    c->num_jobs = num_paths;
    for(int k=0;k<rt.num_threads;k++)
      pthread_pool_task_init(t->task + k, &t->pool, par_compute_neighbours, c);
    pthread_pool_wait(&t->pool);
  }

  //===============================================================================================================
  // for each path vertex index (i.e. bounce) v[2..end]
  // compute the gaussian from neighbourhood information
  c->job = 0;
  c->num_jobs = num_paths;
  for(int k=0;k<rt.num_threads;k++)
    pthread_pool_task_init(t->task + k, &t->pool, par_compute_gaussians, c);
  pthread_pool_wait(&t->pool);
  
  time_build_cdf_gaussians = common_time_wallclock()
                           - time_build_cdf_gaussians;

  if (rm == rm_render)
  {
    printf("final build cdf. num_paths %d, old_path_cnt %d\n", num_paths, old_path_cnt);
#if USE_PATH_CDF==1
  // compute new guide cdf for old guide paths
    c->job = 0;
    c->num_jobs = num_paths;
    for(int k=0;k<rt.num_threads;k++)
      pthread_pool_task_init(t->task + k, &t->pool, par_compute_pdf, c);
    pthread_pool_wait(&t->pool);
#endif
  }

#if USE_PATH_CDF==1
  double time_build_cdf_cdf = common_time_wallclock();

  // build 2d path cdf:
  c->path_cdf_cnt = sqrtf(num_paths)+1;
  c->path_block_cdf_cnt = c->path_cdf_cnt;
  assert(c->path_block_cdf_cnt * c->path_cdf_cnt >= num_paths);
  for(int p=0;p<num_paths;p++)
  {
    const double updf = c->guide_paths[p].uniform_pdf;
    const double gpdf = c->guide_paths[p].guided_pdf;
    const double f = c->guide_paths[p].measurement_contribution;
    double weight = f / (uniform_ratio * updf + (1.f - uniform_ratio) * gpdf);
    if (rm == rm_learn)
      weight = 0.5f * weight + 0.5f * c->guide_paths[p].weight;
    if(!(weight > 0)) weight = 0.0;
    if(!(weight < FLT_MAX)) weight = 0.0;
    c->path_cdf[p] = weight;
    c->guide_paths[p].weight = weight;
    c->guide_paths[p].age++;
    assert(c->path_cdf[p] == c->path_cdf[p]);
  }
  for(int p=num_paths;p<c->path_block_cdf_cnt*c->path_cdf_cnt;p++)
    c->path_cdf[p] = 0.0;
  for(int b=0;b<c->path_block_cdf_cnt;b++)
  {
    for(int p=c->path_cdf_cnt*b+1;p<c->path_cdf_cnt*(b+1);p++)
      c->path_cdf[p] += c->path_cdf[p-1];
    if(c->path_cdf[c->path_cdf_cnt*(b+1)-1] > 0.0)
      for(int p=c->path_cdf_cnt*b;p<c->path_cdf_cnt*(b+1)-1;p++)
        c->path_cdf[p] /= c->path_cdf[c->path_cdf_cnt*(b+1)-1];
    c->path_block_cdf[b] = c->path_cdf[c->path_cdf_cnt*(b+1)-1];
    c->path_cdf[c->path_cdf_cnt*(b+1)-1] = 1.0;
  }
  for(int p=0;p<num_paths;p++)
    assert(c->path_cdf[p] == c->path_cdf[p]);
  for(int p=num_paths;p<c->path_block_cdf_cnt*c->path_cdf_cnt;p++)
    assert(c->path_cdf[p] == c->path_cdf[p]);
  for(int b=1;b<c->path_block_cdf_cnt;b++)
    c->path_block_cdf[b] += c->path_block_cdf[b-1];
  for(int b=0;b<c->path_block_cdf_cnt-1;b++)
    c->path_block_cdf[b] /= c->path_block_cdf[c->path_block_cdf_cnt-1];
  c->path_block_cdf[c->path_block_cdf_cnt-1] = 1.0;
  for(int b=0;b<c->path_block_cdf_cnt;b++)
    assert(c->path_block_cdf[b] == c->path_block_cdf[b]);
  time_build_cdf_cdf = common_time_wallclock()
                     - time_build_cdf_cdf;
#endif
  
#if 1 // debug: count paths based on uniform or guided sampling
  for(int p=old_path_cnt;p<num_paths;p++)
  {
    gpath_t *path = c->guide_paths+p;
    sampling_type_t st = gpath_get_sampling_type(path);

    if (path->measurement_contribution > 0.f) {
      if (st == st_uniform) stats->num_uniform_paths_in_cache++;
      else                  stats->num_guiding_paths_in_cache++;
    }
  }
#endif

  c->num_pathverts = c->num_new_pathverts;
  
  for(int t=0;t<rt.num_threads;t++)
    heap_clear(c->tls[t].heap);
  
  time_build_cdf = common_time_wallclock()
                 - time_build_cdf;
  
#if 0 // debug: output timing breakdown
  double time_build_cdf_rest = time_build_cdf
    - time_build_cdf_recompute
    - time_build_cdf_gaussians
#if USE_PATH_CDF == 1
    - time_build_cdf_cdf
#endif
;

  printf("\n time build cdf %g\n    recompute %g\n    gaussians %g\n    cdf %g\n    copy %g\n    rest %g\n",
      time_build_cdf, 
      time_build_cdf_recompute/time_build_cdf,
      time_build_cdf_gaussians/time_build_cdf,
      time_build_cdf_cdf/time_build_cdf,
      time_build_cdf_rest/time_build_cdf);
#endif
  // guided_collect_stats(s, stdout);
}
  
int guided_debug_sample_around_guide_path(guided_cache_t* c, int pathid, int num_samples, int vcnt, FILE* f)
{
  if (!(c->guide_paths[pathid].measurement_contribution > 0))
    return vcnt;
  const int pathlen = path_len(c, pathid);
  for(int vi=0;vi<pathlen;vi++)
  {
    gvert_t* v = get_vert(c, pathid, vi);
    assert(v);
    if(v->flags & s_environment)
      fprintf(f,  "v %g %g %g\n",
          get_vert(c, pathid, vi-1)->x[0] + v->x[0],
          get_vert(c, pathid, vi-1)->x[1] + v->x[1],
          get_vert(c, pathid, vi-1)->x[2] + v->x[2]);
    else
      fprintf(f,  "v %g %g %g\n", v->x[0], v->x[1], v->x[2]);
  }
  fprintf(f,  "o guide_path_%03d\n", pathid);
  for(int v=0;v<pathlen-1;v++)
    fprintf(f,  "f %d %d\n", vcnt+v, vcnt+v+1);
  vcnt += pathlen;
  
  fprintf(f,  "o guide_path_gaussians_%03d\n", pathid);
  for(int v=1;v<pathlen;v++)
  {
    gvert_t* vp = get_vert(c, pathid, v-1);
    gvert_t* vn = get_vert(c, pathid, v);
    float Bi[9] = {0};
    float D[9] = {0};
    D[0] = vn->w[0];
    D[4] = vn->w[1];
    D[8] = vn->w[2];
    mat3_mul(vn->E, D, Bi);
    float u[3], v[3], n[3] = {
      (vn->flags & s_environment) ? vn->x[0] : (vn->x[0]-vp->x[0]),
      (vn->flags & s_environment) ? vn->x[1] : (vn->x[1]-vp->x[1]),
      (vn->flags & s_environment) ? vn->x[2] : (vn->x[2]-vp->x[2])};
    normalise(n);
    get_onb(n, u, v);

    // marginalize at xp
    float delta[3], uc[3];
    delta[0] = dotproduct(u, vp->x);
    delta[1] = dotproduct(v, vp->x);
    delta[2] = dotproduct(n, vp->x);
    for(int k=0;k<3;k++) delta[k] -= vn->mean[k+3];
    mat3_mulv(vn->S12S22i, delta, uc);
    for(int k=0;k<3;k++) uc[k] += vn->mean[k];

    for(int i=0;i<8;i++)
    {
      const float cutoff = (c->gauss == 3) ? CUTOFF : 4.0f;
      float gt[3], g[] = { (i&1)?cutoff:-cutoff, (i&2)?cutoff:-cutoff, vn->mode & s_volume ? (i&4)?cutoff:-cutoff : 0.f };
      // transform corner to distorted gauss
      mat3_mulv(Bi, g, gt);
      float corner[3];
      // transform corner to world
      for(int k=0;k<3;k++) 
      {
        if(vn->mode & s_volume)
          corner[k] = u[k]*(uc[0]+gt[0]) + v[k]*(uc[1]+gt[1]) + (uc[2]+gt[2])*n[k];
        else
          corner[k] = dotproduct(vn->x, n) * n[k] + (uc[0]+gt[0])*u[k] + (uc[1]+gt[1])*v[k];
        if(vn->flags & s_environment)
          corner[k] += vp->x[k];
      }
      fprintf(f, "v %f %f %f\n", corner[0], corner[1], corner[2]);
    }
    fprintf(f, "f %d %d %d %d\n", vcnt    , vcnt + 1, vcnt + 5, vcnt + 4);//bottom
    fprintf(f, "f %d %d %d %d\n", vcnt + 2, vcnt + 3, vcnt + 7, vcnt + 6);//top
    fprintf(f, "f %d %d %d %d\n", vcnt    , vcnt + 1, vcnt + 3, vcnt + 2);//front
    fprintf(f, "f %d %d %d %d\n", vcnt + 4, vcnt + 5, vcnt + 7, vcnt + 6);//back
    fprintf(f, "f %d %d %d %d\n", vcnt    , vcnt + 4, vcnt + 6, vcnt + 2);//left
    fprintf(f, "f %d %d %d %d\n", vcnt + 1, vcnt + 5, vcnt + 7, vcnt + 3);//right
    vcnt += 8;
  }

  {
    fprintf(f,  "o nearest_paths_%03d\n", pathid);
    for(int k=0;k<NUM_NB;++k)
    {
      int pk = c->nearest_paths[pathid*NUM_NB+k];
      if (pk >= 0)
      {
        gpath_t* gpath = c->guide_paths+pk;
        gvert_t* gverts = c->guide_verts+gpath->vi;
        const int len = gpath->len;
        assert(pathlen == len);
        for(int vi=0;vi<len;vi++)
          if(gverts[vi].flags & s_environment)
            fprintf(f,  "v %g %g %g\n",
                gverts[vi-1].x[0] + gverts[vi].x[0],
                gverts[vi-1].x[1] + gverts[vi].x[1],
                gverts[vi-1].x[2] + gverts[vi].x[2]);
          else
            fprintf(f,  "v %g %g %g\n", gverts[vi].x[0], gverts[vi].x[1], gverts[vi].x[2]);
        for(int v=0;v<len-1;v++)
          fprintf(f,  "f %d %d\n", vcnt+v, vcnt+v+1);
        vcnt += pathlen;
      }
    }
  }

  fprintf(f,  "o sampled_paths_%03d\n", pathid);
  // sample several paths for this guide path
  for (int i = 0; i < num_samples; ++i)
  {
    // sample path
    path_t path;
    path_init(&path, -1ul, 0);
    assert(path.index == -1ul);
    guided_sample(c, &path, pathid);
    
    // write path to file
    const int pathlen = path.length;
    for(int vi=0;vi<pathlen;vi++)
      fprintf(f,  "v %g %g %g\n", path.v[vi].hit.x[0], path.v[vi].hit.x[1], path.v[vi].hit.x[2]);
    for(int v=0;v<pathlen-1;v++)
      fprintf(f,  "f %d %d\n", vcnt+v, vcnt+v+1);
    vcnt += pathlen;
  }
  return vcnt;
}

void guided_debug_path_sampling(guided_cache_t *c, int num_paths)
{
  char filename[512];
  snprintf(filename, sizeof(filename), "sampling_debug.obj");
  FILE* f = fopen(filename, "w");
  int vcnt = 1;

  printf("debug path sampling of best %d of %d paths\n", num_paths, guided_num_guide_paths(c));

  int step = MAX(1, guided_num_guide_paths(c)/num_paths);
  for (int pathid=0;pathid<guided_num_guide_paths(c);pathid+=step)
  {
    vcnt = guided_debug_sample_around_guide_path(c, pathid, 32, vcnt, f);
  }

  fclose(f);
  printf("exported !!! \n");
}

void guided_gpath_hist(
    guided_cache_t* c,
    float* hist,
    int w,
    int h)
{
  int num_paths = guided_num_guide_paths(c);
  
#if 1
  for (int p = 0; p < num_paths; ++p)
  {
    gpath_t* gpath = c->guide_paths+p;
    int px = (int)gpath->pixel_x;
    int py = (int)gpath->pixel_y;
    assert(px < w);
    assert(py < h);
    const uint32_t offset = (px+py*w);
    
    float w = gpath->weight/num_paths;
    hist[3*offset+0] += w;
    hist[3*offset+1] += w;
    hist[3*offset+2] += w;
  }
#else
  for (int p = 0; p < num_paths; ++p)
  {
    gpath_t* gpath = c->guide_paths+p;
    int px = (int)gpath->pixel_x;
    int py = (int)gpath->pixel_y;
    assert(px < w);
    assert(py < h);
    const uint32_t offset = (px+py*w);
    
    float w = ((float)gpath->success)/(gpath->sampled);
    hist[3*offset+0] += w;
    hist[3*offset+1] += 1.f;
  }

  for (int y = 0; y < h; ++y)
  for (int x = 0; x < w; ++x)
  {
    const uint32_t offset = (x+y*w);
    if (hist[3*offset+1] > 0)
    {
      const float w = hist[3*offset+0]/hist[3*offset+1];
      float color[3];
      color_palette(1.f - w, color);

      hist[3*offset+0] = color[0];
      hist[3*offset+1] = color[1];
      hist[3*offset+2] = color[2];
    }
    else
    {
      hist[3*offset+0] = 0.f;
      hist[3*offset+1] = 0.f;
      hist[3*offset+2] = 0.f;
    }
  }
#endif
}

#undef CUTOFF
#undef USE_PATH_CDF
#undef USE_DBOR
#undef RAY_DIFFERENTIALS
#undef RAY_DIFF_MIN
#undef RAY_DIFF_MAX
#undef SIMPLE_GAUSSIAN_SAMPLING
