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
#include "shader.h"
#include <assert.h>

// aggregate shader for surface + interior (with potentially heterogeneous medium).
// usage: `interior <surface id> <interior id>'
// does not pass on any volume callbacks, these have to be passed on in mult shaders
// for the interior shader.
// could potentially be flattened into one fat mult shader instead.

typedef struct
{
  int surface;
  int interior;
}
interior_t;

mf_t sample(path_t *p, void *data)
{
  interior_t *s = (interior_t *)data;
  shader_so_t *surface = rt.shader->shader + s->surface;
  return surface->sample(p, surface->data);
}

mf_t brdf(path_t *p, int v, void *data)
{
  interior_t *s = (interior_t *)data;
  shader_so_t *surface = rt.shader->shader + s->surface;
  return surface->brdf(p, v, surface->data);
}

int init(FILE *s, void **data)
{
  interior_t *t = (interior_t *)malloc(sizeof(interior_t));
  int dreggn = 0;
  *data = t;
  if(fscanf(s, "%d %d", &t->surface, &t->interior) < 2)
  {
    fprintf(stderr, "[interior] could not read all parameters! expecting: <surface id> <interior id>\n");
    return 1;
  }
  shader_so_t *self = (shader_so_t *)data;
  assert(self->data == *data);
  int us = self - rt.shader->shader;
  if(t->surface  < 0) t->surface  = us + t->surface;
  if(t->interior < 0) t->interior = us + t->interior;
  dreggn = fscanf(s, "%*[^\n]\n");
  return dreggn == -1;
}

int shape_init(uint32_t shapeid, shader_so_t *self)
{
  interior_t *s = (interior_t *)self->data;
  shader_so_t *surface = rt.shader->shader + s->surface;
  if(surface->shape_init)
    if(surface->shape_init(shapeid, surface)) return 1;
  shader_so_t *interior = rt.shader->shader + s->interior;
  if(interior->shape_init)
    if(interior->shape_init(shapeid, interior)) return 1;
  return 0;
}

mf_t pdf(path_t *p, int e1, int v, int e2, void *data)
{
  interior_t *s = (interior_t *)data;
  shader_so_t *surface = rt.shader->shader + s->surface;
  return surface->pdf(p, e1, v, e2, surface->data);
}

void cleanup(void *data)
{
  free(data);
}

float prepare(path_t *p, int v, void *data)
{
  interior_t *s = (interior_t *)data;
  shader_so_t *surface = rt.shader->shader + s->surface;
  shader_so_t *interior = rt.shader->shader + s->interior;
  float att = 1.0;
  // act like a surface. the medium case will overwrite the shader
  // from the interior shader we set below and not go through this
  // wrapper any more.
  // need to init v[v].interior, but take care not to interfere with the rest.
  if(interior->prepare) interior->prepare(p, v, interior->data);
  if(surface->prepare) att = surface->prepare(p, v, surface->data);
  // this is only called for the surface case, prepare() on volumes is done on the interior shader directly.
  // we just need to set it, but handle the special case that the volume actually dropped it:
  if(p->v[v].interior.shader != -1)
    p->v[v].interior.shader = s->interior;
  return att;
}

