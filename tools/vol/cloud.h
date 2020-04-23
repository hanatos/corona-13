#pragma once

// some glsl-like convenience functions
static inline float clamp(float x, float lowerlimit, float upperlimit)
{
  if (x < lowerlimit) x = lowerlimit;
  if (x > upperlimit) x = upperlimit;
  return x;
}

static inline float smoothstep(float edge0, float edge1, float x)
{
  x = clamp((x - edge0)/(edge1 - edge0), 0.0, 1.0); 
  return x*x*(3 - 2*x);
}

static inline float fract(float f)
{
  return f-floorf(f);
}

static inline float mix(float f1, float f2, float t)
{
  return f1*(1.0-t) + f2*t;
}

static inline float remap(float val, float m, float M, float nm, float nM)
{
  return clamp(nm + (val-m)/(M-m) * (nM-nm), nm, nM);
}

static inline int mfloor(float f)
{
  return f < 0 ? (int)(f-1.0f) : (int)f;
}

static inline void *read_file(const char *filename, size_t *len)
{
  FILE *f = fopen(filename, "rb");
  if(!f) return 0;
  fseek(f, 0, SEEK_END);
  const size_t filesize = ftell(f);
  fseek(f, 0, SEEK_SET);
  char *file = (char *)malloc(filesize+1);

  size_t rd = fread(file, sizeof(char), filesize, f);
  file[filesize] = 0;
  if(rd != filesize)
  {
    free(file);
    file = 0;
    fclose(f);
    return 0;
  }
  if(len) *len = filesize;
  fclose(f);
  return file;
}


static uint8_t *cloud_tex0 = 0;
static uint8_t *cloud_tex1 = 0;
static uint8_t *cloud_tex2 = 0;

// load static textures
static inline void cloud_init()
{
  cloud_tex0 = read_file("worley.tex", 0);
  cloud_tex1 = read_file("erosion.tex", 0);
  cloud_tex2 = read_file("curl.tex", 0);
  assert(cloud_tex0);
  assert(cloud_tex1);
  assert(cloud_tex2);
}

static inline void texture3d(const uint8_t *const restrict tex, float u, float v, float w, int size, float *restrict res)
{
  float uuf = (u - mfloor(u))*size;
  float vvf = (v - mfloor(v))*size;
  float wwf = (w - mfloor(w))*size;
  int uu = clamp((int)uuf, 0, size-1);
  int vv = clamp((int)vvf, 0, size-1);
  int ww = clamp((int)wwf, 0, size-1);
  int uu1 = uu == size-1 ? 0 : uu+1;
  int vv1 = vv == size-1 ? 0 : vv+1;
  int ww1 = ww == size-1 ? 0 : ww+1;
  float fu = uuf-uu;
  float fv = vvf-vv;
  float fw = wwf-ww;
  for(int k=0;k<4;k++)
    res[k] =
      (1.0-fu)*(1-fv)*(1-fw)*tex[4*(size*(size*ww +vv )+uu )+k]/255.0f+
      (    fu)*(1-fv)*(1-fw)*tex[4*(size*(size*ww +vv )+uu1)+k]/255.0f+
      (1.0-fu)*(  fv)*(1-fw)*tex[4*(size*(size*ww +vv1)+uu )+k]/255.0f+
      (    fu)*(  fv)*(1-fw)*tex[4*(size*(size*ww +vv1)+uu1)+k]/255.0f+
      (1.0-fu)*(1-fv)*(  fw)*tex[4*(size*(size*ww1+vv )+uu )+k]/255.0f+
      (    fu)*(1-fv)*(  fw)*tex[4*(size*(size*ww1+vv )+uu1)+k]/255.0f+
      (1.0-fu)*(  fv)*(  fw)*tex[4*(size*(size*ww1+vv1)+uu )+k]/255.0f+
      (    fu)*(  fv)*(  fw)*tex[4*(size*(size*ww1+vv1)+uu1)+k]/255.0f;
}

static inline void texture2d(const uint8_t *const restrict tex, float u, float v, int size, float *restrict res)
{
  float uuf = (u - mfloor(u))*size;
  float vvf = (v - mfloor(v))*size;
  int uu = clamp((int)uuf, 0, size-1);
  int vv = clamp((int)vvf, 0, size-1);
  int uu1 = (uu+1)%size;
  int vv1 = (vv+1)%size;
  float fu = uuf-uu;
  float fv = vvf-vv;
  for(int k=0;k<4;k++)
    res[k] =
    (1.0-fu)*(1.0-fv)*tex[4*(size*vv+uu)+k]/255.0f+
    (    fu)*(1.0-fv)*tex[4*(size*vv+uu1)+k]/255.0f+
    (1.0-fu)*(    fv)*tex[4*(size*vv1+uu)+k]/255.0f+
    (    fu)*(    fv)*tex[4*(size*vv1+uu1)+k]/255.0f;
}

static inline float cloud_density(const float *const restrict p, const int motion_sample)
{
#ifdef CLOUD_CUT
  if(p[0] < 0) return 0.0f;
#endif
  const float time = 10.0f + 10.0*motion_sample/(VOL_MOTION_SAMPLES-1.0f);
  // hard coded outer shell sphere
  const float c[3] = {0.0, 0.0, 0.0};
  float radius1 = 30.0f;
  float radius0 = 18.0f;
  float feather = 0.1f;

  float dist[3]; // = {p[0]-c[0],p[1]-c[1],p[2]-c[2]};
  for(int k=0;k<3;k++) dist[k] = p[k] - c[k];
  float dist2 = dotproduct(dist, dist);

#if 0 // sphere
  if(dist2 < radius1*radius1) return blend;
  return 0.0;
#endif

  const float cloud_lo = -30, cloud_mid = -15, cloud_hi = 32;
  const float cloud_top_offset = 2;
  const float wind[3] = {8, 0, 1};
  const float cloud_anvil_bias = 0.0;

  // blend using distance from center of earth instead:
  float height = p[2];
  float res[4];
  float envelope = 1;//smoothstep(cloud_lo, cloud_mid, height);// *
                   // smoothstep(6*cloud_hi, cloud_mid, height);
  const float scale = 1/(cloud_hi-cloud_lo);
  const float height_fraction = clamp((height - cloud_lo)/(cloud_hi-cloud_lo), 0.0f, 1.0f);
  float pp[3];
  for(int k=0;k<3;k++)
    pp[k] = (p[k] + height_fraction * wind[k] * cloud_top_offset);
  // float coverage = dist2 < radius1*radius1 ? 0.80 : 0.0; // clamp(fbm_4(0.005*(vec3(2,1,1)*p+wind*u_time)), 0.f, .4f);
  float coverage = clamp((radius1*radius1-dist2)/(radius1*radius1), 0.0, 1.0);
  texture3d(cloud_tex0,
      (pp[0]+0.9*wind[0]*time)*.3*scale,
      (pp[1]+0.9*wind[1]*time)*.3*scale,
      (pp[2]+0.9*wind[2]*time)*.3*scale,
      128,
      res);

  float lowfreq = res[1]*0.625 + res[2]*0.25 + res[3]*0.125;
  float base = res[0];
  float d = remap(base, -(1.0 - lowfreq), 1.0, 0.0, 1.0);
  // d = remap(d, 0.8, 1.0, 0.0, 1.0);

  d *= envelope;

  // coverage = pow(coverage, remap(height_fraction, 0.7, 0.8, 1.0, mix(1.0, 0.5, cloud_anvil_bias)));
  // d = remap(d, coverage, 1.0, 0.0, 1.0);
  d = remap(d, 1-coverage*mix(0.9, 0.2, height_fraction), 1.0, 0.0, 1.0);
  // d = remap(d, mix(0.7, 0.90, height_fraction), 1.0, 0.0, 1.0);
  d *= coverage;
#if 1
#if 1
  texture2d(cloud_tex2,
      100*scale*(pp[0]+.3*wind[0]*time),
      100*scale*(pp[2]+.3*wind[1]*time),
      128,
      res);
  pp[0] += 0.05*res[0];
  pp[1] += 0.05*res[1];
#endif
  texture3d(cloud_tex1,
      (p[0]+.4*wind[0]*time)*scale*10.,
      (p[1]+.4*wind[1]*time)*scale*10.,
      (p[2]+.4*wind[2]*time)*scale*10.,
      32,
      res);
  float hfreqn = res[3];//mix(res[3], 1.0-res[3], clamp(height_fraction*10, 0, 1));
  d = remap(d, hfreqn * 0.2, 1.0, 0.0, 1.0);
#endif

  assert(d == d && d < FLT_MAX && d>= 0);

#if !defined CLOUD_FULL && !defined CLOUD_INNER && !defined CLOUD_OUTER
#define CLOUD_FULL
#endif

#ifdef CLOUD_FULL
  const float blend = 1.0f;
#endif
  // const float blo = 0.45, bhi = 0.50; // small inner core
  const float blo = 0.15, bhi = 0.20; // large inner core
#ifdef CLOUD_INNER
  const float blend = smoothstep(blo, bhi, d);
#endif
#ifdef CLOUD_OUTER
  const float blend = 1.0f-smoothstep(blo, bhi, d);
#endif
  // float blend = smoothstep(radius0*radius0, (radius0+feather)*(radius0+feather), dist2);
  return blend * d;
}

static inline float canonical_density(const float *const restrict p, const int motion_sample)
{
  if(p[1] > -24) return 0;
  // ramp from left to right: increase density: T(s) = scale(s) * exp(-sum(densities)) => density + log(scale(s))
  // ramp top down: increase added noise
  const float x = (p[0]+32)/64.0;
  const float y = (p[2]+32)/64.0;

  float d = 24*(0.01 + logf(1.0f+x));
  float yy = clamp(-0.1 + logf(1.0f+y), 0, 100);
  float q[3] = {1000+p[0], 1000+p[1], 1000+p[2]};
  float n = .5 * d * pnoise(q, 9, 3.331f);

  return clamp(d + yy*n, 0, 10000.0f);
}
