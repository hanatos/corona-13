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

#include "corona_common.h"
#include "shader.h"
#include "prims.h"

int init(FILE *s, void **data)
{
  *data = 0;
  // read sample shader number
  if(fscanf(s, "%*d%*[^\n]\n"))
    fprintf(stderr, "[normals] could not read host number!\n");
  return 0;
}

float prepare(path_t *p, int v, void *data)
{
#if 0
  // bezier like anti-terminator-hack:
  // adjust hit point to avoid terminator problem (hooray, dirty hacks rock...)
  float hitu[3], hitv[3], hitw[3];
  for(int k=0;k<3;k++)
  {
    hitu[k] = hit->hit[k] - rt->prim[hit->prim].v[2][k];
    hitv[k] = hit->hit[k] - rt->prim[hit->prim].v[1][k];
    hitw[k] = hit->hit[k] - rt->prim[hit->prim].v[0][k];
  }
  const float dotu = fminf(0.0f, dotproduct(hitu, n+6)),  dotv = fminf(0.0f, dotproduct(hitv, n+2)),  dotw = fminf(0.0f, dotproduct(hitw, n));
  for(int k=0;k<3;k++)
  {
    hitu[k] -= dotu*n[6+k];
    hitv[k] -= dotv*n[3+k];
    hitw[k] -= dotw*n[k];
  }
  for(int k=0;k<3;k++) hit->hit2[k] = (1-hit->u-hit->v)*(rt->prim[hit->prim].v[0][k] + hitw[k]) + hit->v*(rt->prim[hit->prim].v[1][k] + hitv[k]) + hit->u*(rt->prim[hit->prim].v[2][k] + hitu[k]);
#endif
  return 1.0f;
}

