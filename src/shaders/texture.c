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
#include "texture.h"
#include "framebuffer.h"

#include <math.h>
#include <assert.h>

typedef struct tex_t
{
  int num;
  tex_slot_t slot;
  framebuffer_t fb;
  float mul;
}
tex_t;

int init(FILE *s, void **data)
{
  tex_t *t = malloc(sizeof(*t));
  memset(t, 0, sizeof(*t));
  *data = t;
  char filename[1024];
  char c;
  t->mul = 1.0f;
  if(fscanf(s, " %c %s %f", &c, filename, &t->mul) < 2)
  {
    fprintf(stderr, "[texture] shader could not read all parameters! expecting: <dsevgrt> [filename] [mul]\n");
    return 1;
  }
  int dreggn = fscanf(s, "%*[^\n]\n");
  if(filename[0] == '#') filename[0] = '\0';
  t->slot = tex_parse_slot(c);

  if(fb_map(&t->fb, filename))
  {
    fprintf(stderr, "[texture] could not load framebuffer `%s'!\n", filename);
    return 1;
  }
  t->fb.retain = 1;

  return dreggn == -1;
}

void cleanup(void *data)
{
  tex_t *t = data;
  fb_cleanup(&t->fb);
}

float prepare(path_t *p, int v, void *data)
{
  tex_t *s = data;
  // float px[4];
  // fb_fetch(&s->fb, p->v[v].hit.s, p->v[v].hit.t, px);
  const float *px = fb_fetch(&s->fb, p->v[v].hit.s, p->v[v].hit.t);
  if(s->fb.header->channels == 4) p->v[v].shading.t = px[3];
  tex_set_slot_coeff(p, v, s->slot, s->mul, px);
  if(s->fb.header->channels == 4 && p->v[v].shading.t < 0.5f)
  {
    p->v[v].material_modes = p->v[v].mode = s_transmit | s_specular;
    return -1.0f;
  }
  return 1.0f;
}

