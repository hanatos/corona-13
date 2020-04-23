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

typedef struct
{
  int num;
  int *pre;
  int host;
}
mult_t;

mf_t sample(path_t *p, void *data)
{
  mult_t *s = (mult_t*) data;
  shader_so_t *host = rt.shader->shader + s->host;
  return host->sample(p, host->data);
}

mf_t brdf(path_t *p, int v, void *data)
{
  mult_t *s = (mult_t*) data;
  shader_so_t *host = rt.shader->shader + s->host;
  return host->brdf(p, v, host->data);
}

int volume_enabled(void *data)
{
  mult_t *s = (mult_t*) data;
  shader_so_t *host = rt.shader->shader + s->host;
  if(host->volume_enabled)
    return host->volume_enabled(host->data);
  else return 0;
}

mf_t volume_transmittance(path_t *p, int e, void *data)
{
  mult_t *s = (mult_t*) data;
  shader_so_t *host = rt.shader->shader + s->host;
  return host->volume_transmittance(p, e, host->data);
}

mf_t volume_pdf(const path_t *p, int e, void *data)
{
  mult_t *s = (mult_t*) data;
  shader_so_t *host = rt.shader->shader + s->host;
  return host->volume_pdf(p, e, host->data);
}

mf_t volume_pdf_adjoint(const path_t *p, int e, void *data)
{
  mult_t *s = (mult_t*) data;
  shader_so_t *host = rt.shader->shader + s->host;
  return host->volume_pdf_adj(p, e, host->data);
}

float volume_sample(path_t *p, int e, void *data)
{
  mult_t *s = (mult_t*) data;
  shader_so_t *host = rt.shader->shader + s->host;
  return host->volume_sample(p, e, host->data);
}

int init(FILE *s, void **data)
{
  mult_t *t = (mult_t *)malloc(sizeof(mult_t));
  int dreggn = 0;
  *data = t;
  // read sample shader number
  if(fscanf(s, "%d", &(t->num)) < 1)
  {
    fprintf(stderr, "[mult] shader could not read all parameters! expecting: <num>\n");
    return 1;
  }
  t->pre = (int*)malloc(t->num * sizeof(int));
  for(int k=0;k<t->num;k++)
  {
    if(fscanf(s, "%d", t->pre + k) < 1)
    {
      fprintf(stderr, "[mult] shader failed to read %d-th shader num!\n", k);
      return 1;
    }
  }
  if(fscanf(s, "%d", &t->host) < 1)
  {
    fprintf(stderr, "[mult] shader failed to read host num!\n");
    return 1;
  }
  dreggn = fscanf(s, "%*[^\n]\n");
  shader_so_t *self = (shader_so_t *)data;
  assert(self->data == *data);

  // get our shader number and offset the others if negative numbers:
  int us = self - rt.shader->shader;
  if(t->host < 0) t->host = us + t->host;
  for(int k=0;k<t->num;k++) if(t->pre[k] < 0) t->pre[k] = us + t->pre[k];

  shader_so_t *host = rt.shader->shader + t->host;
  if(!host->volume_enabled || !host->volume_enabled(host->data))
  {
    // kill our callbacks
    self->volume_transmittance = 0;
    self->volume_sample = 0;
    self->volume_pdf = 0;
    self->volume_pdf_adj = 0;
  }
  return dreggn == -1;
}

int shape_init(uint32_t shapeid, shader_so_t *self)
{
  mult_t *s = self->data;
  shader_so_t *host = rt.shader->shader + s->host;
  for(int k=0;k<s->num;k++)
  {
    shader_so_t *pre = rt.shader->shader + s->pre[k];
    if(pre->shape_init && pre->shape_init(shapeid, pre))
      return 1;
  }
  if(host->shape_init && host->shape_init(shapeid, host)) return 1;
  return 0;
}

mf_t pdf(path_t *p, int e1, int v, int e2, void *data)
{
  mult_t *s = (mult_t*) data;
  shader_so_t *host = rt.shader->shader + s->host;
  return host->pdf(p, e1, v, e2, host->data);
}

void cleanup(void *data)
{
  mult_t *s = (mult_t *)data;
  free(s->pre);
}

float prepare(path_t *p, int v, void *data)
{
  mult_t *s = (mult_t *)data;
  float att = 1.0f;
  for(int k=0;k<s->num;k++)
  {
    shader_so_t *pre = rt.shader->shader + s->pre[k];
    if(pre->prepare)
      att *= pre->prepare(p, v, pre->data);
  }
  shader_so_t *host = rt.shader->shader + s->host;
  if(host->prepare) att *= host->prepare(p, v, host->data);
  return att;
}

