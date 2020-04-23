/*
    This file is part of corona-6: radiata.

    corona-6: radiata is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    corona-6: radiata is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with corona-6: radiata.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "prim.h"
#include "accel.h"

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#ifdef _OPENMP
  #include <omp.h>
#endif


#define ACCEL_ADD_TO_CELL(c, p) \
      if(cnt)                                                             \
      {                                                                   \
        if(c->prim != primid+1)                                             \
        {                                                                 \
          c->prim = primid+1;                                               \
          c->num_prims++;                                                 \
        }                                                                 \
      }                                                                   \
      else                                                                \
      {                                                                   \
        int k=0;                                                          \
        for(;k<c->num_prims;k++) if(g->index[c->prim + k] == primid) break;\
        if(k == c->num_prims)\
        {\
          assert(g->index_size > c->prim + c->num_prims);\
          g->index[c->prim + c->num_prims++] = primid;\
          accel_populate(g, p);\
        }\
      }

void rasterize(grid_t *g, int primid, const int cnt)
{
  // 20*4: float has 23 bits, last 3 are always garbled.
  prim_t stack[80];
  int top = 0;
  stack[0] = g->prim[primid];
  float aabb[6];
  int p1[3], p2[3];

  const float term = 0.1f*g->grid_size;

  while (top >= 0)
  {
    prim_get_aabb(stack+top, aabb);

    accel_gridpos(g, aabb,   p1);
    accel_gridpos(g, aabb+3, p2);

    const int l1norm = p2[2] + p2[1] + p2[0] - p1[2] - p1[1] - p1[0];
    // TODO: replace by only one coordinate > 0? for i=min..max add.
    if(l1norm == 1)
    {
      // 2 voxels or subdivided enough.
      // add primid to both hashes, remove prim, continue
      hashcell_t *c1 = accel_hashcell(g, aabb);
      hashcell_t *c2 = accel_hashcell(g, aabb+3);
      ACCEL_ADD_TO_CELL(c1, p1)
      ACCEL_ADD_TO_CELL(c2, p2)
      top --;
    }
    else if(l1norm == 0)
    {
      // only one voxel
      hashcell_t *c = accel_hashcell(g, aabb);
      ACCEL_ADD_TO_CELL(c, p1)
      top--;
    }
    else if(((aabb[3]-aabb[0] < term) && (aabb[4]-aabb[1] < term) && (aabb[5]-aabb[2] < term)) || top >= 15)
    {
      int corner[3], ipos[3];
      float fc[3];
      for(corner[0]=0;corner[0]<4;corner[0]+=3)
      for(corner[1]=0;corner[1]<4;corner[1]+=3)
      for(corner[2]=0;corner[2]<4;corner[2]+=3)
      {
        for(int k=0;k<3;k++) fc[k] = aabb[k+corner[k]];
        accel_gridpos(g, fc, ipos);
        hashcell_t *c = accel_hashcell(g, fc);
        ACCEL_ADD_TO_CELL(c, ipos);
      }
      top--;
    }
    else
    {
      // complex overlap
      // TODO: possible optimization using voxel cache: find all candidate voxels (rasterize aabb) and check against
      // voxels already added to.
      prim_subdivide(stack+top, stack+top, stack+top+1, stack+top+2, stack+top+3);
      top += 3;
    }
  }
}

accel_t* accel_init(const int num_prims, prim_t *prim)
{
  accel_t *g = (accel_t *)malloc(sizeof(accel_t));
  g->prim = prim;
  g->num_prims = num_prims;
  g->index = (int *)malloc(sizeof(int)*num_prims);

  // round hashtable size up to next power of 2.
  int i = (int)(1.3*g->num_prims) + 1;
  int k = 0x40000000;
  for(;k>i;k>>=1);
  g->num_entries = k<<1;
  //g->num_entries >>= 2;
  printf("[accel_init] hashtable has %d entries\n", g->num_entries);
  g->hashtable = (hashcell_t *)malloc(g->num_entries*sizeof(hashcell_t));
  return g;
}

void accel_cleanup(struct accel_t *g)
{
  free(g->hashtable);
  free(g->gridmask);
  free(g->index);
  free(g);
}

void accel_debug(struct accel_t *g)
{
  int cnt = 0;
  int maxc = 0;
  int nonempty = 0;
  for(int k=0;k<g->num_entries;k++) 
  {
    if(g->hashtable[k].num_prims) nonempty++;
    cnt += g->hashtable[k].num_prims;
    if(g->hashtable[k].num_prims > maxc) maxc = g->hashtable[k].num_prims;
  }
  printf("hashtable with %d/%d = %f fill, %d max\n", cnt, g->num_entries, cnt/(float)nonempty, maxc);
  ray_t ray;
  rayhit_t hit;
  hit.threadid = 0;
  hit.points = rt.points;
  FILE *f = fopen("dreggn0.ppm", "wb");
  fprintf(f, "P6\n%d %d\n255\n", g->num_cells[1], g->num_cells[2]);
  ray.dir[2] = ray.dir[1] = 0.0f;
  ray.dir[0] = 1.0f;
  for(int k=0;k<g->num_cells[2];k++)
  for(int j=0;j<g->num_cells[1];j++)
  {
    ray.pos[2] = (k)*g->grid_size + g->aabb[2];
    ray.pos[1] = (j)*g->grid_size + g->aabb[1];
    ray.pos[0] = g->aabb[0];
    unsigned char outc[3] = {0, 0, 0};
    for(int i=0;i<g->num_cells[0];i++)
    {
      int p[3] = {i, j, k};
      if(accel_populated(g, p))
      {
        hit.dist = INFINITY;
        hit.prim = -1;
        hashcell_t *c = accel_hashcelli(g, p);
        for(int l=c->prim;l<c->prim+c->num_prims;l++)
        //for(int l=0;l<g->num_prims;l++)
          prim_intersect(g->prim + g->index[l], &ray, &hit, g->index[l]);
        if(hit.prim != -1)
        {
          prim_get_normal(g->prim + hit.prim, &ray, &hit);
          normalise(hit.normal);
          for(int l=0;l<3;l++) outc[l] = fmaxf(0.0f, fminf(255.0f, 128.0f*(hit.normal[l] + 1.0f)));
          break;
        }
      }
    }
    fwrite(&outc, 3, sizeof(char), f);
  }
  fclose(f);

  f = fopen("dreggn1.pgm", "wb");
  fprintf(f, "P5\n%d %d\n255\n", g->num_cells[0], g->num_cells[2]);
  for(int k=0;k<g->num_cells[2];k++)
  for(int j=0;j<g->num_cells[0];j++)
  {
    int out = 0;
    for(int i=0;i<g->num_cells[1];i++)
    {
      int p[3] = {j, i, k};
      if(accel_populated(g, p)) out += 20;
    }
    out = out > 255 ? 255 : out;
    unsigned char outc = out;
    fwrite(&outc, 1, sizeof(char), f);
  }
  fclose(f);

  f = fopen("dreggn2.pgm", "wb");
  fprintf(f, "P5\n%d %d\n255\n", g->num_cells[0], g->num_cells[1]);
  for(int k=0;k<g->num_cells[1];k++)
  for(int j=0;j<g->num_cells[0];j++)
  {
    int out = 0;
    for(int i=0;i<g->num_cells[2];i++)
    {
      int p[3] = {j, k, i};
      if(accel_populated(g, p)) out += 20;
    }
    out = out > 255 ? 255 : out;
    unsigned char outc = out;
    fwrite(&outc, 1, sizeof(char), f);
  }
  fclose(f);
}

void accel_build(struct accel_t *g, const char* filename)
{
  // create aabb
  float aabb[6*rt.num_threads];
  for(int k=0;k<3;k++)
  {
    g->aabb[k] = INFINITY;
    g->aabb[k+3] = -INFINITY;
    for(int i=0;i<rt.num_threads;i++)
    {
      aabb[6*i + k] = INFINITY;
      aabb[6*i + k + 3] = -INFINITY;
    }
  }
//#pragma omp parallel for default(none) shared(g, aabb)
  for(int k=0;k<g->num_prims;k++)
  {
    int threadid = 0;
#ifdef _OPENMP
    threadid = omp_get_thread_num();
#endif
    float tmp[6];
    prim_get_aabb(g->prim + k, tmp);
    for(int k=0;k<3;k++) if(tmp[k] < aabb[6*threadid+k]) aabb[6*threadid+k] = tmp[k];
    for(int k=0;k<3;k++) if(tmp[k+3] > aabb[6*threadid+k+3]) aabb[6*threadid+k+3] = tmp[k+3];
  }
  for(int i=0;i<rt.num_threads;i++)
  {
    for(int k=0;k<3;k++) if(g->aabb[k] > aabb[6*i+k]) g->aabb[k] = aabb[6*i+k];
    for(int k=0;k<3;k++) if(g->aabb[k+3] < aabb[6*i+k+3]) g->aabb[k+3] = aabb[6*i+k+3];
  }

  // FIXME: dynamic models do not allow for dynamic mem alloc here!!!!

  // create grid: O(num_tris) cells, O(cells) hash table entries
  // thiago ize uses N_i = d_i sqrt3f( 5* N / V), d the diagonal and V the volume
  const float lambda = 3.f;
  const float volume = (g->aabb[0+3]-g->aabb[0])*(g->aabb[1+3]-g->aabb[1])*(g->aabb[2+3]-g->aabb[2]);
  g->grid_size = 0.0f;
  for(int k=0;k<3;k++)
  {
    g->num_cells[k] = (int)((g->aabb[k+3]-g->aabb[k])*powf(lambda*g->num_prims/volume, 0.3333f));
    g->grid_size = fmaxf(g->grid_size, (g->aabb[k+3]-g->aabb[k])/g->num_cells[k]);
  }
  for(int k=0;k<3;k++)
  {
    g->num_cells[k]++;
    g->grid_factor = 1.0f/g->grid_size;
  }

  printf("[accel_build] initing a %dx%dx%d grid\n", g->num_cells[0], g->num_cells[1], g->num_cells[2]);
  printf("[accel_build] grid size: %f \n", g->grid_size);

  for(int k=0;k<3;k++) assert(g->grid_size*g->num_cells[k] >= g->aabb[k+3] - g->aabb[k]);

  const int masksize = MAX(4, (sizeof(unsigned int)*g->num_cells[0]*g->num_cells[1]*g->num_cells[2])>>5);
  g->gridmask = (unsigned int *)malloc(masksize);
  memset(g->hashtable, 0, g->num_entries*sizeof(hashcell_t));
  memset(g->gridmask, 0, masksize);

  // construct array-hashtable (haines99)
#if 0
  accel_t gt[rt.num_threads];
  for(int k=0;k<rt.num_threads;k++)
  {
    gt[k].prim = g->prim;
    memcpy(gt[k].aabb, g->aabb, 6*sizeof(float));
    gt[k].grid_factor = g->grid_factor;
    gt[k].grid_size = g->grid_size;
    gt[k].num_entries = g->num_entries;
    gt[k].hashtable = (hashcell_t *)malloc(g->num_entries*sizeof(hashcell_t));
    memset(gt[k].hashtable, 0, g->num_entries*sizeof(hashcell_t));
  }
  // rasterize prims to grid, count refs: O(n/p)
  // TODO: avoid this step using memory predictor models?
//#pragma omp parallel for default(none) shared(g, gt) schedule(dynamic)
  for(int k=0;k<g->num_prims;k++)
  {
    int threadid = 0;
#ifdef _OPENMP
    threadid = omp_get_thread_num();
#endif
    // TODO: lower precision/ use only aabb for mem-conservative faster check?
    rasterize(gt+threadid, k, 1);
  }

  // O(n/p)
  int num_prim_refs[rt.num_threads+1];
  memset(num_prim_refs, 0, sizeof(int)*(rt.num_threads+1));
//#pragma omp parallel for default(none) shared(num_prim_refs, rt, gt) shared(g) schedule(static)
  for(int k=0;k<g->num_entries;k++)
  
    int threadid = 0;
#if 0//def _OPENMP
    threadid = omp_get_thread_num();
#endif
    g->hashtable[k].prim = num_prim_refs[threadid+1];
    for(int i=0;i<rt.num_threads;i++) num_prim_refs[threadid+1] += gt[i].hashtable[k].num_prims;
    g->hashtable[k].num_prims = 0;
  }

  // update hashtable.prim (+= num_prim_refs[t<=t0])
  /*
#pragma omp parallel for default(none) shared(num_prim_refs) shared(g) schedule(static)
  for(int k=0;k<g->num_entries;k++)
  {
    int threadid = 0;
#ifdef _OPENMP
    threadid = omp_get_thread_num();
#endif
    for(int k=0;k<=threadid;k++) g->hashtable[k].prim += num_prim_refs[k];
  }*/

  // TODO: find elegant reduction for this:
  g->index_size = 0;
  //for(int k=0;k<rt.num_threads;k++) g->index_size += num_prim_refs[k+1];
  // malloc ref array (num_prim_refs)
  g->index = (int *)malloc(sizeof(int)*g->index_size);
  printf("[accel_build] constructing index array with %d entries\n", g->index_size);
#endif

  for(int k=0;k<g->num_prims;k++)
    rasterize(g, k, 1);

  int pos = 0;
  for(int k=0;k<g->num_entries;k++)
  {
    g->hashtable[k].prim = pos;
    pos += g->hashtable[k].num_prims;
    g->hashtable[k].num_prims = 0;
  }

  g->index_size = pos;
  g->index = (int *)malloc(sizeof(int)*g->index_size);
  printf("[accel_build] constructing index array with %d entries\n", g->index_size);

  // O(n/p)
//#pragma omp parallel for default(none) shared(g) schedule(static)
  for(int k=0;k<g->num_prims;k++)
  {
    rasterize(g, k, 0);
  }

  //for(int k=0;k<rt.num_threads;k++) free(gt[k].hashtable);
  accel_debug(g);
}

void accel_intersect(const struct accel_t *g, const ray_t *ray, rayhit_t *hit)
{
  // small cache (hit rate for sponza ~ 0.002)
  // const int mailboxsize = 16;
  // int mailbox[mailboxsize];
  // int mailboxpos = 0;

  // TODO: SIMD this. perhaps like that:
  // init t[3], t[3] + 1,2,3 cells.
  // step t+4, collect tri indices (accel_populated)
  // if > 4, simd intersect primitives, early out if < t[3], fill mailbox cache

  // clip ray to aabb and rasterize inside
  int p[3], mind = 0;
  const int step[3] = {ray->dir[0] > 0 ? 1 : -1, ray->dir[1] > 0 ? 1 : -1, ray->dir[2] > 0 ? 1 : -1 };
  float t[3];
  float invdir[3];

  float tmin = 0.0f;
  float tmax = FLT_MAX;
  for(int k=0;k<3;k++)
  {
    invdir[k] = 1.0f/ray->dir[k];
    tmin = fmaxf(tmin, (g->aabb[k+(ray->dir[k] > 0 ? 0 : 3)] - ray->pos[k])*invdir[k]);
    tmax = fminf(tmax, (g->aabb[k+(ray->dir[k] < 0 ? 0 : 3)] - ray->pos[k])*invdir[k]);
  }
  if(tmin > tmax) return;
  hit->dist = tmax + 0.001f;
  //if(tmin > 0.0f) tmin += 0.5f*g->grid_size;
  for(int k=0;k<3;k++) t[k] = ray->pos[k] + tmin*ray->dir[k];
  accel_gridpos(g, t, p);
  //for(int k=0;k<3;k++) if(p[k] < 0) p[k] = 0;
  const float min[3] = {g->aabb[0] + g->grid_size*(p[0]-.5f), g->aabb[1] + g->grid_size*(p[1]-.5f), g->aabb[2] + g->grid_size*(p[2]-.5f)};
  const float max[3] = {g->grid_size+min[0], g->grid_size+min[1], g->grid_size+min[2]};
  for(int k=0;k<3;k++) t[k] = (ray->dir[k] > 0.0f ? max[k] - ray->pos[k] : min[k] - ray->pos[k])*invdir[k];

  const float delta[3] = {g->grid_size*fabsf(invdir[0]), g->grid_size*fabsf(invdir[1]), g->grid_size*fabsf(invdir[2]) };
  // bad idea: empty space sponza: 8.6->7.6 drop
  //int hash[3] = {p[0]*73856093, p[1]*19349663, p[2]*83492791};
  //const int hashstep[3] = {step[0]*73856093, step[1]*19349663, step[2]*83492791};
  //printf("\nnew ray\n");
  while(1)
  {
    // rather cheap texture lookup..
    if(accel_populated(g, p))
    {
      // rather expensive hash function calculation
      //const hashcell_t *c = g->hashtable + ((hash[0] ^ hash[1] ^ hash[2]) & (g->num_entries - 1));
      const hashcell_t *c = accel_hashcelli(g, p);
        //hit->dist = t[mind];
        //hit->prim = 0;//c->prim;
        //hit->u = hit->v = 0.3f;
        //return;
        int tmp = c->prim;
        for(int k=0;k<c->num_prims;k++)
        {
          /*for(int k=0;k<mailboxpos && k<mailboxsize;k++)
          {
            rt.accel->cachelookups++;
            if(mailbox[k] == g->index[tmp]) { rt.accel->cachehits++; goto next_prim;}
          }
          mailbox[mailboxpos++ & (mailboxsize-1)] = g->index[tmp];*/
          prim_intersect(g->prim + g->index[tmp], ray, hit, g->index[tmp]);
//next_prim:
          tmp++;
        }
        //return;
      }
    mind = 0;
    if(t[2] < t[1] && t[2] < t[0]) mind = 2;
    else if(t[1] < t[0]) mind = 1;

    if(t[mind] > hit->dist) return;
    t[mind] += delta[mind];
    p[mind] += step[mind];
    //hash[mind] += hashstep[mind];
      //if(p[0] == g->num_cells[0] - 1 || p[1] == g->num_cells[1] - 1 || p[2] == g->num_cells[2] - 1)
      //printf("leaf at (%d %d %d) with %d prims at %d\n", p[0], p[1], p[2], c->num_prims, c->prim);
  }
}

int accel_visible(const struct accel_t *g, const ray_t *ray)
{
  //tmin = 0, tmax = 1 then aabbtest
  return 1;
}

