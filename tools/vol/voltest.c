#include "vol/vol.h"
#include "vol/trace.h"
#include "pnoise.h" // for perlin noise
#include "cloud.h"

#include <stdlib.h>
#include <stdio.h>

// linkage for random number generator
rt_t rt;
__thread rt_tls_t rt_tls;

#ifndef VOLTEST_RES
#define VOLTEST_RES 1
#endif
#if VOLTEST_RES == 2
static const float voxel_size = 1./64.0f; // insane, like 12GB
#elif VOLTEST_RES == 1
static const float voxel_size = 1./8.0f; // 1.1 GB
#else
static const float voxel_size = 1.0f; // 2.1 MB
#endif

static inline float torus_density(
    const float *p,   // point to evaluate
    const float *c0,  // center of the torus
    const float r,    // small radius
    const float R,    // big radius
    const float *n0)  // normal/orientation
{
  float tx = p[0], ty = p[1], tz = p[2];

  const float xt0[3] = {tx-c0[0], ty-c0[1], tz-c0[2]};
  float a0[3], b0[3];
  get_onb(n0, a0, b0);
  const float x0[3] = {dotproduct(xt0, a0), dotproduct(xt0, b0), dotproduct(xt0, n0)};
  const float norm0 = sqrtf(x0[0]*x0[0] + x0[1]*x0[1]);
  const float dd0[3] = {x0[0]/norm0 * R - x0[0], x0[1]/norm0 * R - x0[1], 0.0f - x0[2]};
  return r - sqrtf(dotproduct(dd0, dd0));
}

static inline float fake_rand(uint32_t *seed)
{
  *seed *= 0x15A4E35u;
  return (float) (*seed >> 16) / (0xFFFF + 1.0f);
}

static inline float test_density(const float *p, float *temperature, int motion_sample)
{
#if 1 // cloud
  *temperature = 0;
  return cloud_density(p, motion_sample);
  // return canonical_density(p, motion_sample);
#endif
#if 0 // looping rain on the sun
  const float scale = 3.0f;
  const int octaves = 4;
  const float q0[3] = {p[0] + 70000, p[1] + 10000, p[2] + 23000};
  const float q1[3] = {p[0] + 10000, p[1] + 70000, p[2] + 10000};
  const float q2[3] = {p[0] + 13000, p[1] + 10000, p[2] + 70000};
  float ng[3] = { pnoise(q0, octaves, scale), pnoise(q1, octaves, scale), pnoise(q2, octaves, scale)};
  const float radius = 90.0f;
  const float thickness = 6.f;
  const float ns = .5f;
  const float c[3] = {0.f, 0.f, -100.f};
  for(int k=0;k<3;k++) ng[k] *= ns;
  const float pp[3] =
  {
    p[0] - 0.5f*ng[0] - c[0],
    p[1] - 0.5f*ng[1] - c[1],
    p[2] - 0.5f*ng[2] - c[2]
  };
  const float ppp[3] =
  {
    p[0] - c[0],
    p[1] - c[1],
    p[2] - c[2]
  };
  const float add = fabsf(20.0*pnoise(ng, octaves, 1.0f));

  int outside = sqrtf(ppp[0]*ppp[0] + ppp[1]*ppp[1] + ppp[2]*ppp[2]) > radius + 0.9*thickness ? 1.f : 0.f;
  int cull_inside = sqrtf(ppp[0]*ppp[0] + ppp[1]*ppp[1] + ppp[2]*ppp[2]) < radius ? 1.f : 0.f;
  if(cull_inside)
  {
    *temperature = 0.0f;
    return 0.0f;
  }
  float val = 0.0f;
  // constant emission:
  // if(fabsf(radius - sqrtf(pp[0]*pp[0] + pp[1]*pp[1] + pp[2]*pp[2])) > thickness) val = 1.0f;
  // more emission at border:
  float dist = 
  CLAMP((MAX(0.0, sqrtf(pp[0]*pp[0] + pp[1]*pp[1] + pp[2]*pp[2])-radius))/thickness, 0.0, 2.0);
  *temperature = 0.0f;
  if(dist < .96)
  {
    val = dist*dist*dist*dist;
    const float norm = val/(val + 15.0f);
    *temperature = MAX(180.0*val, 1600.0+norm * 1800.0f);
    float noise = 30*pnoise(p, 5, 10);
    *temperature = MAX(0, *temperature + noise);
    dist = 0.96f-dist;
    val = dist*dist;
  }

  // float val = MAX(0.0, thickness-fabsf(radius - sqrtf(pp[0]*pp[0] + pp[1]*pp[1] + pp[2]*pp[2])))/thickness;
  // val = MAX(0.0, radius - sqrtf(pp[0]*pp[0] + pp[1]*pp[1] + pp[2]*pp[2]));
  // loops:
  const int num = 70;
  uint32_t seed = 0x15A4E35;
  float c0[3] = {10.0f, 0.0f, -4.0f};
  float r = 0.5f;
  float tval = 0.0f;
  for(int k=0;k<num;k++)
  {
    const float f = k/(float)num;
    const float f2 = f*f;
    const float f4 = f2*f2;
    float R = 6.0f+4.0f*f4 + fake_rand(&seed)*1.2f;
    c0[1] = 4.0*f2 + fake_rand(&seed)*1.4f;
    c0[2] = -1.0 + f4*5.0 + fake_rand(&seed)*0.6f;
    const float tilt = 0.5f*k/num + fake_rand(&seed)*1.5;
    const float nn = sqrtf(1.f + tilt*tilt);
    const float n0[3] = {1.0f/nn, -tilt/nn, 0.0f};
    float tdens = MAX(0.0, torus_density(p, c0, r, R, n0)/r);
    tdens *= tdens;
    tval += tdens;
  }
  // TODO: also a couple orthogonal for the straight rain?

  // mix and set temperature (loops are hotter)
  val *= 20.0f;
  tval *= 2.4f;
  if(val > 0.0f) val += 0.8*add;
  if(tval > 0.0f) tval += 0.1*add;

#if 0
  { // sun
    const float norm = val/(val + 15.0f);
    *temperature = MAX(180.0*val, 1600.0+norm * 1800.0f);
    float noise = 50*pnoise(p, 2, 10);
    *temperature = MAX(0, *temperature + noise);
  }
#endif
  if (outside) { // loop
    const float norm = tval/(tval + 15.0f);
    *temperature = MAX(700*tval, 2000.0+norm * 3200.0f);
  }
  val += tval;

  // TODO: add some scattering atmosphere?
  return val;
#else
#if 1
#if 1
  if(p[0] < -16 || p[0] > 16 ||
     p[1] < -16 || p[1] > 16 ||
     p[2] < -16 || p[2] > 16) return 0.0f; // early out
  const float scale = 5;
  const int octaves = 3;
  const float q0[3] = {p[0] + 70000, p[1] + 10000, p[2] + 23000};
  const float q1[3] = {p[0] + 10000, p[1] + 70000, p[2] + 10000};
  const float q2[3] = {p[0] + 13000, p[1] + 10000, p[2] + 70000};
  const float ng[3] = { pnoise(q0, octaves, scale), pnoise(q1, octaves, scale), pnoise(q2, octaves, scale)};
  const float pp[3] = {p[0] - ng[0], p[1] - ng[1], p[2] - ng[2]};
  const float add = 1e-3f*(ng[0]*ng[0]+ng[1]*ng[1]+ng[2]*ng[2]);
#else
  if(p[0] < - 6 || p[0] > 26 ||
     p[1] < -16 || p[1] > 16 ||
     p[2] < -16 || p[2] > 16) return 0.0f; // early out
  const float pp[3] = {p[0]-10, p[1]-0, p[2]-0};
  const float add = 0.0f;
#endif
  const float val = fmaxf(0.0, 4.0f-fabsf(10.0f - sqrtf(pp[0]*pp[0] + pp[1]*pp[1] + pp[2]*pp[2])));
  *temperature = 800.0*(add+val);
  assert(*temperature >= 0.f);
  if(val > 0.0f) return add+val;
  return 0.0f;
#else
  if(p[0] < -4 || p[0] > 4 ||
     p[1] < -32 || p[1] > 32 ||
     p[2] < -4 || p[2] > 4) return 0.0f; // early out
  *temperature = 2000.f;
  float v = 3.f*(p[1]+32.f)/64.f;
  v = expf(v) - 1;
  assert(v > 0.f);
  return v;
#endif
#endif
}

static inline void motion(float t, float *off)
{
#if 0
#if 0 // about the linear motion the eulerian version had:
  for(int k=0;k<3;k++) off[k] = 0.0f;
  off[1] = 30.0f * (t-.5f);
#else // interesting spline
  const float x0[3] = {-16,  16, -16},
              x1[3] = {-16,  16,  16},
              x2[3] = { 16,  16,  16},
              x3[3] = { 16, -16,  16};
  for(int k=0;k<3;k++)
    off[k] = (1.0-t)*(1.0-t)*(1.0-t)*x0[k] + 3.0*(1.0-t)*(1.0-t)*t*x1[k] +
             3.0*(1.0-t)*t*t*x2[k] + t*t*t*x3[k];
#endif
#else
  // no motion blur:
  for(int k=0;k<3;k++) off[k] = 0.0f;
#endif
}

vol_payload_type_t payload_fill(
    void *data, vol_payload_uncompressed_t *payload,
    const float *aabb, const int force_static)
{ // struct comes in memset(0), only write what's needed.
  vol_payload_type_t fill = s_vol_empty;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) default(shared)
#endif
  for(int i=0;i<512;i++)
  {
    vol_index_t vx = {0};
    vx.idx = i;
    // sample input, fill leaf struct
    float x[3];
    x[0] = aabb[0] + (aabb[0+3]-aabb[0])*(vx.i+.5f)/8.0f;
    x[1] = aabb[1] + (aabb[1+3]-aabb[1])*(vx.j+.5f)/8.0f;
    x[2] = aabb[2] + (aabb[2+3]-aabb[2])*(vx.k+.5f)/8.0f;

    float off[3];
    float temperature = 0.0f, density = 0.0f;
    if(force_static)
    {
      density = test_density(x, &temperature, 0);
      if(density > 0.0) fill = s_vol_static;
      assert(density == density && density < FLT_MAX && density >= 0);
      payload->d[i][0] = density;
      payload->t[i][0] = temperature;
    }
    else
    {
      for(int s=0;s<VOL_MOTION_SAMPLES;s++)
      {
        float t = s/(VOL_MOTION_SAMPLES-1.0f);
        motion(t, off);
        for(int l=0;l<3;l++) off[l] += x[l];
        density = test_density(off, &temperature, s);
        if(density > 0.0) fill = s_vol_full;
        assert(density == density && density < FLT_MAX && density >= 0);
        payload->d[i][s] = density;
        payload->t[i][s] = temperature;
      }
    }
  }
  return fill;
}

int main (int argc, char *argv[])
{
  if(argc < 2)
  {
    fprintf(stderr, "usage: voltest output.vol [-s,--static] [--shader-id x]\n");
    exit(1);
  }

  int force_static = 0;
  int shader_id = 0;
  if(argc > 2) for(int a=2;a<argc;a++)
  {
    if(argv[a][0] == '-' && (argv[a][1] == 's' || !strcmp(argv[a], "--static")))
      force_static = 1;
    if(argv[a][0] == '-' && !strcmp(argv[a], "--shader-id") && (++a < argc))
      shader_id = atol(argv[a]);
  }

  // load helper textures
  cloud_init();
  
  // init random numbers:
  rt.points = points_init(1, 0);

  float aabb[6];
  aabb[0] = -32;
  aabb[1] = -32;
  aabb[2] = -32;
  aabb[3] = 32;
  aabb[4] = 32;
  aabb[5] = 32;
  // XXX debug: comment out if large tree already constructed:
  if(vol_create_tree(argv[1], aabb, voxel_size, 0, 0,
        &payload_fill, 0, shader_id, force_static))
    exit(1);

  vol_tree_t *tree = vol_open(argv[1]);
  assert(tree);
  FILE *f;

  exit(0);

  const float mu_t = 0.01f;
  const int lod = 1;


  // test ray tracing transmittance:
  char filename[512];
  int num = 32;
  float *buf = malloc(sizeof(float)*3*512*512);
  // for(int interpolation=0;interpolation<4;interpolation++)
  for(int interpolation=0;interpolation<1;interpolation++)
  {
    memset(buf, 0, sizeof(float)*3*512*512);
    double start = common_time_wallclock();
    int idx = 0;
    for(int j=0;j<512;j++)
    for(int i=0;i<512;i++)
    {
      for(int s=0;s<num;s++)
      {
        const float pos[3] = {0, 0, -400}; // TODO: anti aliasing?
        float dir[3] = {-31.5+i*64.0/512.0, -31.5+j*64.0/512.0, 400};
        normalise(dir);
        const float time = drand48();
        const float t = vol_trace_transmittance(tree, pos, dir, 500, mu_t, interpolation, lod, time, 550, 0, 0);
        for(int k=0;k<3;k++) buf[idx + k] += t/num;
      }
      idx += 3;
    }
    double end = common_time_wallclock();
    fprintf(stderr, "time to process transmittance %d: %g seconds.\n", interpolation, end-start);
    snprintf(filename, sizeof(filename), "transmittance_%02d.pfm", interpolation);
    f = fopen(filename, "wb");
    fprintf(f, "PF\n%d %d\n-1.0\n", 512, 512);
    fwrite(buf, 3*sizeof(float), 512*512, f);
    fclose(f);
  }

  // test ray tracing distance sampling
  num = 64;
  for(int interpolation=0;interpolation<2;interpolation++)
  {
    memset(buf, 0, sizeof(float)*3*512*512);
    double start = common_time_wallclock();
    int idx = 0;
    for(int j=0;j<512;j++)
    for(int i=0;i<512;i++)
    {
      for(int s=0;s<num;s++)
      {
        const float pos[3] = {0, 0, -400}; // TODO: anti aliasing?
        float dir[3] = {-31.5+i*64.0/512.0, -31.5+j*64.0/512.0, 400};
        normalise(dir);
        const float time = drand48();
        const float t = vol_trace_sample(tree, pos, dir, 500, mu_t, drand48(), interpolation, lod, time, 550.0, 0,0,0,0);
        if(t == FLT_MAX)
          for(int k=0;k<3;k++) buf[idx + k] += 1.0f/num;
      }
      idx += 3;
    }
    double end = common_time_wallclock();
    fprintf(stderr, "time to process sampling %d: %g seconds.\n", interpolation, end-start);
    snprintf(filename, sizeof(filename), "sampling_%02d.pfm", interpolation);
    f = fopen(filename, "wb");
    fprintf(f, "PF\n%d %d\n-1.0\n", 512, 512);
    fwrite(buf, 3*sizeof(float), 512*512, f);
    fclose(f);
  }

  // test ray marching + world space sampling
  num = 64;
  for(int interpolation=0;interpolation<2;interpolation++)
  {
    memset(buf, 0, sizeof(float)*3*512*512);
    double start = common_time_wallclock();
    int idx = 0;
    for(int j=0;j<512;j++)
    for(int i=0;i<512;i++)
    {
      const float pos[3] = {0, 0, -400}; // TODO: anti aliasing?
      float dir[3] = {-31.5+i*64.0/512.0, -31.5+j*64.0/512.0, 400};
      normalise(dir);
      const float tmin = 400.0/dir[2], tmax = tmin + 80.0;
      float thickness = 0.0f;
      for(int s=0;s<num;s++)
      {
        const float time = drand48();
        float x[3];
        for(int k=0;k<3;k++) x[k] = pos[k] + (tmin + (tmax-tmin)*s/(num-1.0f)) * dir[k];
        float res[3];
        int empty = vol_lookup(tree, x[0], x[1], x[2], s_vol_density, interpolation, lod, time, res);
        if(!empty) thickness += res[0] * (tmax-tmin)/(num-1.0f);
      }
      for(int k=0;k<3;k++) buf[idx + k] = expf(-thickness * mu_t);
      idx += 3;
    }
    double end = common_time_wallclock();
    fprintf(stderr, "time to process marching %d: %g seconds.\n", interpolation, end-start);
    snprintf(filename, sizeof(filename), "marching_%02d.pfm", interpolation);
    f = fopen(filename, "wb");
    fprintf(f, "PF\n%d %d\n-1.0\n", 512, 512);
    fwrite(buf, 3*sizeof(float), 512*512, f);
    fclose(f);
  }

  free(buf);

  vol_close(tree);
  exit(0);
}
