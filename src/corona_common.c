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
#include "prims.h"
#include "points.h"
#include "accel.h"
#include "camera.h"
#include "view.h"
#include "render.h"
#include "sampler.h"
#include "spectrum.h"
#include "shader.h"
#include "version.h"

int common_load_scene(FILE *f)
{
  int num_shapes = 0;
  if(!fscanf(f, "%d\n", &num_shapes))
  {
    fprintf(stderr, "[common_load_scene] corrupt model file: could not read number of shapes!\n");
    return 1;
  }
  if(num_shapes == 0) return 0;
  char line[BUFSIZ];

  int shader;
  char filename[512];
  char texture[512];

  prims_allocate(rt.prims, num_shapes);
  for(int shape=0;shape<num_shapes;shape++)
  {
    int items_read = fscanf(f, "%[^\n]\n", line);
    if(items_read == -1) fprintf(stderr, "\n");
    snprintf(texture, 512, "none");
    if(sscanf(line, "%d %s %s\n", &shader, filename, texture) < 2)
    {
      fprintf(stderr, "[common_load_scene] WARN: malformed line (%d): %s\n", shape + rt.shader->num_shaders + 1, line);
      continue;
    }
    if(shader >= rt.shader->num_shaders)
      fprintf(stderr, "[common_load_scene] WARN: shader %d in line %d (%s) out of bounds!\n", shader, shape, filename);
    if(shader < 0 || shader >= rt.shader->num_shaders) shader = 0;
    prims_load(rt.prims, filename, texture, shader);

    // load uv/normals etc:
    int discard = shader_shape_init(shape, rt.shader->shader + shader);
    if(discard)
      prims_discard_shape(rt.prims, shape);
  }
  prims_allocate_index(rt.prims);
  return 0;
}

void common_write_sidecar(const char *tempname)
{
  char filename[1024];
  sprintf(filename, "%s.txt", tempname);
  FILE *f = fopen(filename, "wb");
  fprintf(f, "corona-13: %s\n", VERSION);
  fprintf(f, "file     : %s\n", rt.basename);
  const float *aabb = accel_aabb(rt.accel);
  if(aabb[3] < aabb[0])
    fprintf(f, "aabb     : empty\n");
  else
    fprintf(f, "aabb     : (%.3f, %.3f)x(%.3f, %.3f)x(%.3f, %.3f) dm^3\n",
      aabb[0], aabb[3], aabb[1],
      aabb[4], aabb[2], aabb[5]);
  points_print_info(f);
  prims_print_info(f);
  accel_print_info(f);
  shader_print_info(f);
  render_print_info(f);
  view_print_info(f);
  sampler_print_info(f);
  pointsampler_print_info(f);
  fprintf(f, "camera   : ");
  colour_camera_print_info(f);
  fprintf(f, "input    : ");
  colour_input_print_info(f);
  fclose(f);
}

