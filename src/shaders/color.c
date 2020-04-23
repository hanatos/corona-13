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
#include "matrix3.h"
#include "colourspaces.h"
#include "rgb2spec.h"
#include "shader.h"
#include "spectrum.h"
#include "texture.h"
#include "lights.h"

typedef struct
{
  float coeff[3];
  float mul;
  float roughness;
  tex_slot_t slot;
}
color_shader_t;

int init(FILE *s, void **data)
{
  color_shader_t *t = malloc(sizeof(*t));
  *data = t;
  float col[3];
  char c;
  t->slot = s_slot_diffuse;
  t->roughness = 1.0f; // cos^1
  if(fscanf(s, " %c %f %f %f %f", &c, col, col+1, col+2, &(t->roughness)) < 4)
  {
    fprintf(stderr, "[color] could not read all parameters! expecting: [dgsevr] r g b [roughness]\n");
    free(t);
    *data = 0;
    return 1;
  }
  t->mul = spectrum_rgb_to_coeff(col, t->coeff);
  int dreggn = 0;
  dreggn = fscanf(s, "%*[^\n]\n");
  t->slot = tex_parse_slot(c);
  if(t->slot == s_slot_count)
  {
    fprintf(stderr, "[color] could not parse texture slot %c!\n", c);
    free(t);
    *data = 0;
    return 1;
  }
  return dreggn == -1;
}

int shape_init(uint32_t shapeid, shader_so_t *self)
{
  color_shader_t *t = self->data;
  const __m128 l4 = _mm_set_ps(400.0f, 480.0f, 560.0f, 660.0f);
  const __m128 eval = rgb2spec_eval_sse(t->coeff, l4);
  if(t->slot == s_slot_emission)
    lights_init_light(shapeid, t->mul*(eval[0]+eval[1]+eval[2]+eval[3])/4.0f);
  return 0;
}

float prepare(path_t *p, int v, void *data)
{
  color_shader_t *s = (color_shader_t *)data;
  p->v[v].shading.roughness = s->roughness;
  tex_set_slot_coeff(p, v, s->slot, s->mul, s->coeff);
  return 1.0f;
}

