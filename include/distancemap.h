#ifndef _DISTANCE_MAP_H
#define _DISTANCE_MAP_H

#include "corona_common.h"
#include <stdint.h>

// temporary struct
typedef struct distancemap_t
{
  uint16_t *data;
  int wd, ht, dp;
}
distancemap_t;

static inline uint16_t distancemap_get(
    distancemap_t *m,
    int x,
    int y,
    int z,
    int c)
{
  return m->data[c + 3*(x + m->wd*(y + m->ht*z))];
}

static inline void distancemap_set(
    distancemap_t *m,
    int x,
    int y,
    int z,
    int c,
    uint16_t val)
{
  m->data[c + 3*(x + m->wd*(y + m->ht*z))] = val;
}

static inline int distancemap_inside(
    distancemap_t *m,
    int x,
    int y,
    int z)
{
  return (0 <= x && x < m->wd &&
      0 <= y && y < m->ht &&
      0 <= z && z < m->dp); 
}

static inline void distancemap_combine(
    distancemap_t *m,
    int x, int y, int z,
    int dx, int dy, int dz)
{
  // compute on torus in x/y:
  int xx = (x + dx)%m->wd;
  int yy = (y + dy)%m->ht;
  int zz = z + dz;
  if(distancemap_inside(m, x, y, z) && distancemap_inside(m, xx, yy, zz))
  {
    uint16_t d[3]; d[0] = abs(dx); d[1] = abs(dy); d[2] = abs(dz);

    uint32_t v1[3], v2[3];
    for (int i = 0; i < 3; i++) v1[i] = distancemap_get(m, x, y, z, i);
    for (int i = 0; i < 3; i++) v2[i] = distancemap_get(m, xx, yy, zz, i) + d[i];

    if(dotproduct(v1, v1) > dotproduct(v2, v2))
      for (int i = 0; i < 3; i++)
        distancemap_set(m, x, y, z, i, v2[i]);
  }
}

#if 0
static inline void distancemap_combine(
    distancemap_t *m,
    int dx, int dy, int dz,
    int cx, int cy, int cz,
    int x, int y, int z)
{
  while (distancemap_inside(m, x, y, z) && distancemap_inside(m, x + dx, y + dy, z + dz))
  {
    uint16_t d[3]; d[0] = abs(dx); d[1] = abs(dy); d[2] = abs(dz);

    uint32_t v1[3], v2[3];
    for (int i = 0; i < 3; i++) v1[i] = distancemap_get(m, x, y, z, i);
    for (int i = 0; i < 3; i++) v2[i] = distancemap_get(m, x+dx, y+dy, z+dz, i) + d[i];

    if(dotproduct(v1, v1) > dotproduct(v2, v2))
      for (int i = 0; i < 3; i++)
        distancemap_set(m, x, y, z, i, v2[i]);

    x += cx; y += cy; z += cz;
  }
}
#endif

static inline void _distancemap_progress(float p)
{
  char spinner[] = {'|', '/', '-', '\\'};
  static int state = 0;
  const int sp_size = sizeof(spinner)/sizeof(char);
  fprintf(stderr, "\r[dmap] %c %04.1f%%", spinner[state], p);
  state = (state + 1)%sp_size;
}

// only public function here.
static inline float *distancemap_init(
    int wd,
    int ht,
    int dp,
    float depth,
    float *heightmap)
{
  distancemap_t m;
  m.wd = wd; m.ht = ht; m.dp = dp;
  m.data = (uint16_t *)malloc(sizeof(uint16_t)*wd*ht*dp*3);
  for (int y = 0; y < ht; y++)
  {
    for (int x = 0; x < wd; x++)
    {
      for (int z = 0; z < dp; z++)
      {
        // init everything below surface to 0 and cells
        // above the surface to their distance straight down (it'll only be made shorter later on)
        const int above = z - heightmap[x + wd*y] * dp;
        if (above <= 0 || z <= 1)
          for (int i = 0; i < 3; i++)
            distancemap_set(&m, x, y, z, i, 0u);
        else
        {
          distancemap_set(&m, x, y, z, 0, 0);
          distancemap_set(&m, x, y, z, 1, 0);
          distancemap_set(&m, x, y, z, 2, above);
        }
      }
    }
  }

#if 0
  // search bottom to top:
  for(int z=1;z<dp;z++)
  {
    fprintf(stderr, ".");
    for(int y=0;y<ht;y++)
    {
      for(int x=0;x<wd;x++)
      {
        // walk all possible cells which might contain something
        // potentially closer than what we have already
        int max = MAX(MAX(distancemap_get(&m, x, y, z, 0),
              distancemap_get(&m, x, y, z, 1)),
            distancemap_get(&m, x, y, z, 2));
        for(int yy=MAX(0, y-max);yy<MIN(ht-1, y+max);yy++)
        for(int xx=MAX(0, x-max);xx<MIN(wd-1, x+max);xx++)
        {
          distancemap_combine(&m, x, y, z, xx-x, yy-x, -1);
          // make sure we clamp that as we go:
          max = MAX(MAX(distancemap_get(&m, x, y, z, 0),
                distancemap_get(&m, x, y, z, 1)),
              distancemap_get(&m, x, y, z, 2));
          if(yy < y-max) yy = y-max;
          if(xx < x-max) xx = x-max;
        }
      }
    }
  }
  fprintf(stderr, "X");
#endif

    // TODO: modulo/torus?
#if 1
  // plane 0254 seems to jump up by one pixel:
  for(int z=1;z<dp;z++)
  {
    _distancemap_progress(100*z/(dp-1.0f));
    const int dz = -1; // previous plane was alredy perfectly filled (nearest hit is always below or on same level)
    // collect everything in this plane:
    for(int y=0;y<ht;y++)
      for(int x=0;x<wd;x++)
      {
        for(int ii=-1;ii<=1;ii++) for(int jj=-1;jj<=1;jj++)
          distancemap_combine(&m, x, y, z, ii, jj, dz);
       distancemap_combine(&m, x, y, z, -1, -1, 0);
       distancemap_combine(&m, x, y, z, -1,  0, 0);
       distancemap_combine(&m, x, y, z,  0, -1, 0);
       distancemap_combine(&m, x, y, z,  1, -1, 0);
      }
    // and propagate back again:
    for(int y=ht-1;y>=0;y--)
      for(int x=wd-1;x>=0;x--)
      {
        for(int ii=-1;ii<=1;ii++) for(int jj=-1;jj<=1;jj++)
          distancemap_combine(&m, x, y, z, ii, jj, dz);
        distancemap_combine(&m, x, y, z, 1, 1, 0);
        distancemap_combine(&m, x, y, z, 1, 0, 0);
        distancemap_combine(&m, x, y, z, 0, 1, 0);
        distancemap_combine(&m, x, y, z, -1, 1, 0);
      }
  }
  fprintf(stderr, "\n");
#endif
#if 0
  for(int z=1;z<dp;z++)
  {
    _distancemap_progress(100*z/(dp-1.0f));
    // do some diffusion in this z plane, get info from below:
    for(int k=0;k<4;k++) // k==4 seems to result in reasonably isotropic distance maps. k==2 still shows some ringing.
    {
    for(int y=0;y<ht;y++)
      for(int x=0;x<wd;x++)
       distancemap_combine(&m, x, y, z, -1, 0, -1);
    for(int x=0;x<wd;x++)
      for(int y=0;y<ht;y++)
        distancemap_combine(&m, x, y, z, 0, -1, -1);
    // and back
    for(int y=0;y<ht;y++)
      for(int x=wd-1;x>=0;x--)
        distancemap_combine(&m, x, y, z, 1, 0, -1);
    for(int x=0;x<wd;x++)
      for(int y=ht-1;y>=0;y--)
        distancemap_combine(&m, x, y, z, 0, 1, -1);
    }
  }
  fprintf(stderr, "\n");
#if 0
  // not needed
  for(int z=1;z<dp;z++)
  {
    fprintf(stderr, ".");
    // do some diffusion in between z planes.
    // lower minimum distances are only to be expected from z <<
    for(int y=0;y<ht;y++)
      for(int x=0;x<wd;x++)
        distancemap_combine(&m, x, y, z, 0, 0, -1);
  }
  fprintf(stderr, "X\n");
#endif
#endif

#if 0
  fprintf(stderr, "[distance map] pass one");
  for (int z = 1; z < dp; z++)
  {
    fprintf(stderr, ".");

    // combine with everything with dz = -1
    for (int y = 0; y < ht; y++)
    {
      distancemap_combine(&m,
          0, 0, -1,
          1, 0, 0,
          0, y, z);
    }
    
    for (int y = 1; y < ht; y++)
    {
      distancemap_combine(&m,
          0, -1, 0,
          1, 0, 0,
          0, y, z);
      distancemap_combine(&m,
          -1, 0, 0,
          1, 0, 0,
          1, y, z);
      distancemap_combine(&m,
          +1, 0, 0,
          -1, 0, 0,
          wd - 2, y, z);
    }

    for (int y = ht - 2; y >= 0; y--)
    {
      distancemap_combine(&m,
          0, +1, 0,
          1, 0, 0,
          0, y, z);
      distancemap_combine(&m,
          -1, 0, 0,
          1, 0, 0,
          1, y, z);
      distancemap_combine(&m,
          +1, 0, 0,
          -1, 0, 0,
          wd - 2, y, z);
    }
  }
  fprintf(stderr, "done\n");
  fprintf(stderr, "[distance map] pass two");

  for (int z = dp - 2; z >= 0; z--)
  {
    fprintf(stderr, ".");

    // combine with everything with dz = +1
    for (int y = 0; y < ht; y++)
    {
      distancemap_combine(&m,
          0, 0, +1,
          1, 0, 0,
          0, y, z);
    }
    
    for (int y = 1; y < ht; y++)
    {
      distancemap_combine(&m,
          0, -1, 0,
          1, 0, 0,
          0, y, z);
      distancemap_combine(&m,
          -1, 0, 0,
          1, 0, 0,
          1, y, z);
      distancemap_combine(&m,
          +1, 0, 0,
          -1, 0, 0,
          wd - 2, y, z);
    }
    for (int y = ht - 2; y >= 0; y--)
    {
      distancemap_combine(&m,
          0, +1, 0,
          1, 0, 0,
          0, y, z);
      distancemap_combine(&m,
          -1, 0, 0,
          1, 0, 0,
          1, y, z);
      distancemap_combine(&m,
          +1, 0, 0,
          -1, 0, 0,
          wd - 2, y, z);
    }
  }
  
  fprintf(stderr, "done\n");
#endif

  float *dmap = (float *)common_alloc(16, sizeof(float)*wd*ht*dp);
  for (int z = 0; z < dp; z++)
  {
    for (int y = 0; y < ht; y++)
    {
      for (int x = 0; x < wd; x++)
      {
        float value = 0;
        value += distancemap_get(&m, x, y, z, 0)*distancemap_get(&m, x, y, z, 0)/((float)wd*wd);
        value += distancemap_get(&m, x, y, z, 1)*distancemap_get(&m, x, y, z, 1)/((float)ht*ht);
        value += distancemap_get(&m, x, y, z, 2)*distancemap_get(&m, x, y, z, 2)/((float)dp*dp);
        value = sqrtf(value);
        dmap[x + wd*(y + ht*z)] = value;
      }
    }
  }
  free(m.data);
  return dmap;
}

#endif
