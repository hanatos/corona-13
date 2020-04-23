#pragma once
#include "corona_common.h" // for fast log and aligned alloc
#include "screenshot.h" // for fast log and aligned alloc
#include "filter.h"

#include <assert.h>
#include <string.h>

#define DBOR_DS 4

// fast implementation of density based outlier rejection
typedef struct dbor_t
{
  int stops;           // how many stops difference between buffers
  int num_buffers;     // how many buffers
  uint64_t buf_width;  // dimensions of buffer
  uint64_t buf_height;
  float *data;         // consecutive data backing
  float **throughputs; // pointer to the buffer begins
}
dbor_t;

// clear samples from buffers
static inline void dbor_clear(dbor_t *d)
{
  memset(d->data, 0, sizeof(float)*d->num_buffers*d->buf_width*d->buf_height);
}

// init a new dbor cascaded frame buffer
static inline dbor_t *dbor_init(int wd, int ht, int stops, int bufs)
{
  dbor_t *d = malloc(sizeof(*d));
  d->stops = stops;
  d->num_buffers = bufs;
  d->buf_width  = wd/DBOR_DS;
  d->buf_height = ht/DBOR_DS;
  d->data = common_alloc(16, sizeof(float)*d->num_buffers*d->buf_width*d->buf_height);
  d->throughputs = malloc(sizeof(float*)*d->num_buffers);
  for(int k=0;k<d->num_buffers;k++)
    d->throughputs[k] = d->data + d->buf_width*d->buf_height*k;
  dbor_clear(d);
  return d;
}

static inline void dbor_cleanup(dbor_t *d)
{
  free(d->data);
  free(d->throughputs);
  free(d);
}

// evaluate the trust value/duplicity of a sample
static inline float dbor_trust(const dbor_t *d, float x, float y, float throughput)
{
#if 1
  // const float logval = common_fasterlog2(throughput+1.0f)*.5f;
  const float logval = MAX(0, log2f(throughput));
  const int l = CLAMP((int)logval, 0, d->num_buffers-1);
  const int ll = CLAMP(l-1, 0, d->num_buffers-1);
  const int u = CLAMP(l+1, 0, d->num_buffers-1);
  int i = CLAMP((int)(x/DBOR_DS), 0, d->buf_width-1);
  int j = CLAMP((int)(y/DBOR_DS), 0, d->buf_height-1);
  if (ll == l || l == u)
  {
    return d->throughputs[ll][i+j*d->buf_width]
         + d->throughputs[ u][i+j*d->buf_width];
  }
  return d->throughputs[ll][i+j*d->buf_width]
       + d->throughputs[ l][i+j*d->buf_width]
       + d->throughputs[ u][i+j*d->buf_width];
#else
  if (!(throughput > 0.f))
    return 0.f;

  const float logval = MAX(0, log2f(throughput));
  const int i = CLAMP((int)logval, 0, d->num_buffers-1);
  float n = 0, n_avg = 0;
  for (int dh = -1; dh < 2; ++dh)
  for (int dw = -1; dw < 2; ++dw)
  {
    const int offset = CLAMP(y+dh, 0, d->buf_height-1) * d->buf_width + CLAMP(x+dw, 0, d->buf_width-1);
    float n_tmp = 0;
    if (i == 0) 
    {
      n_tmp = d->throughputs[  i][offset] / (1 << i)
            + d->throughputs[i+1][offset] / (1 << (i+1));
    }
    else if (i == (d->num_buffers - 1)) 
    {
      n_tmp = d->throughputs[i-1][offset] / (1 << (i-1))
            + d->throughputs[  i][offset] / (1 << i);
    }
    else 
    {
      n_tmp = d->throughputs[i-1][offset] / (1 << (i-1))
            + d->throughputs[  i][offset] / (1 << i)
            + d->throughputs[i+1][offset] / (1 << (i+1));
    }
    
    n_avg += n_tmp;
    if (dh == 0 && dw == 0)
      n = n_tmp;
  }
  n_avg /= 9.f;
  return n_avg;
#endif
}

static inline void dbor_filter_splat(
    float *buf,
    const float i,
    const float j,
    const float val,
    const int wd,
    const int ht)
{
  const int ii = (int)i;
  const int jj = (int)j;
  const float fi = i - ii;
  const float fj = j - jj;
  const int iip = ii+1, jjp = jj+1;
  float *p = buf + ii+wd*jj;
  if(ii  >= 0 && jj  >= 0 && ii  < wd && jj  < ht)
    common_atomic_add(p     , val * (1.0f-fi)*(1.0f-fj));
  if(iip >= 0 && jj  >= 0 && iip < wd && jj  < ht)
    common_atomic_add(p+1   , val * (     fi)*(1.0f-fj));
  if(ii  >= 0 && jjp >= 0 && ii  < wd && jjp < ht)
    common_atomic_add(p  +wd, val * (1.0f-fi)*(     fj));
  if(iip >= 0 && jjp >= 0 && iip < wd && jjp < ht)
    common_atomic_add(p+1+wd, val * (     fi)*(     fj));
}

// this is thread-safe since it uses atomics in filter_splat (also it depends
// on the pixel filter you chose in config.mk)
static inline float dbor_splat(dbor_t *d, float x, float y, float throughput)
{
#if 1
  if (throughput > (1<<d->num_buffers)) return 0.f;
  assert(throughput > 0.f);
  const float logval = MAX(0, log2f(throughput));
  const int l = CLAMP((int)logval, 0, d->num_buffers-1);
  const int u = CLAMP(l+1, 0, d->num_buffers-1);
  const float lv = (l == d->num_buffers-1) || throughput < 1.f ? 1.f :
    (((1<<l)/throughput)-0.5f)/0.5f;
  const float uv = 1.f - lv;
  dbor_filter_splat(d->throughputs[l], x/DBOR_DS, y/DBOR_DS, lv,
      d->buf_width, d->buf_height);
  dbor_filter_splat(d->throughputs[u], x/DBOR_DS, y/DBOR_DS, uv,
      d->buf_width, d->buf_height);
  return dbor_trust(d, x, y, throughput);
#else
  if (throughput > 0)
  {
    const float logval = MAX(0, log2f(throughput));
    const int l = CLAMP((int)logval, 0, d->num_buffers-1);
    const int u = l+1;
    const float v = (((1<<l)/throughput)-0.5f)/0.5f;
    const float lv = CLAMP(throughput < 1.f ? 1.f : v, 0.0f, 1.0f);
    const float uv = 1.f - lv;
    const float lt = lv * throughput;
    const float ut = uv * throughput;

    filter_splat(x, y, &lt, d->throughputs[l], 1, 0, 1, d->buf_width, d->buf_height);
    if (u < d->num_buffers)
      filter_splat(x, y, &ut, d->throughputs[u], 1, 0, 1, d->buf_width, d->buf_height);
  }
  return dbor_trust(d, x, y, throughput);
#endif
}

static inline void dbor_export(const dbor_t *d, const char* str, int num_samples)
{
  // write dbor buffer to file
  for(int k=0;k<d->num_buffers;k++)
  {
    char filename[1024];
    snprintf(filename, sizeof(filename), "%s_%d", str, k);
    screenshot_write(filename, d->throughputs[k], 1, 0, 1, d->buf_width, d->buf_height, 1.f / num_samples);
  }
}

#undef DBOR_DS
