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

#ifndef SHADER_H
#define SHADER_H

#include "corona_common.h"
#include "quaternion.h"
#include "pathspace.h"
#include "pointsampler.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <float.h>
#include <math.h>

/* shader callbacks for brdf */
typedef mf_t  (*sample_t)(path_t *p, void *data);
typedef int   (*inverse_sample_t)(const path_t *p, const int v, float *r_omega_x, float *r_omega_y, float *r_scatter_mode, void *data);
typedef mf_t  (*brdf_t)(path_t *p, int v, void *data);
typedef int   (*init_t)(FILE *s, void**);
typedef mf_t  (*pdf_t)(const path_t *p, int e1, int v, int e2, void*);
typedef void  (*cleanup_t)(void *);

/* volume shaders */
typedef mf_t  (*volume_transmittance_t)(path_t *p, int e, void *data);
typedef float (*volume_sample_t)(path_t *p, int e, void *data);
typedef mf_t  (*volume_pdf_t)(const path_t *p, int e, void *data);
typedef int   (*volume_enabled_t)(void *data);

/* sample end point v for next event estimation */
typedef mf_t  (*volume_sample_nee_t)(path_t *p, int v, void *data);
typedef mf_t  (*volume_pdf_nee_t)(const path_t *p, int v, void *data);

/* sample direction p->e[e].{omega,pdf} for forward next event estimation */
typedef mf_t  (*volume_sample_fnee_direction_t)(path_t *p, int e, void *data);
typedef mf_t  (*volume_pdf_fnee_direction_t)(const path_t *p, int e, void *data);

/* sample start point v[0] for light tracing */
typedef mf_t  (*volume_sample_light_t)(path_t *p, void *data);
typedef mf_t  (*volume_pdf_light_t)(const path_t *p, int v, void *data);

struct shader_so_t;
/* initialise shader for given shape */
typedef int   (*shape_init_t)(uint32_t shapeid, struct shader_so_t *self);
typedef float (*prepare_t)(path_t *p, int v, void *data);

/* envmap specific */
typedef mf_t  (*sky_eval_t)(const path_t *p, int v, void *);
typedef mf_t  (*sky_sample_t)(path_t *p, void *);
typedef mf_t  (*sky_pdf_t)(const path_t *p, int v, void *);


typedef struct
{
  void* data;
  sky_eval_t eval;
  init_t init;
  cleanup_t cleanup;
  sky_sample_t sample;
  sky_pdf_t pdf;
  int black;
  int sun;
}
shader_sky_t;

typedef struct shader_so_t
{
  void *data;
  prepare_t prepare;
  shape_init_t shape_init;
  sample_t sample;
  inverse_sample_t inverse_sample;
  brdf_t brdf;
  init_t init;
  cleanup_t cleanup;
  pdf_t pdf;

  volume_transmittance_t    volume_transmittance;
  volume_sample_t           volume_sample;
  volume_pdf_t              volume_pdf;
  volume_pdf_t              volume_pdf_adj;
  volume_enabled_t          volume_enabled;

  volume_sample_nee_t       volume_sample_nee;
  volume_pdf_nee_t          volume_pdf_nee;
  
  volume_sample_fnee_direction_t volume_sample_fnee_direction;
  volume_pdf_fnee_direction_t    volume_pdf_fnee_direction;
  volume_pdf_fnee_direction_t    volume_pdf_fnee_direction_adj;

  volume_sample_light_t     volume_sample_light;
  volume_pdf_light_t        volume_pdf_light;
}
shader_so_t;

typedef struct shader_t
{
  void* handle[64];
  int num_handles;
  char dlname[64][20];
  shader_sky_t skyshader;
  int num_shaders;
  shader_so_t *shader;
  int exterior_medium_shader; // shader id which holds the exterior medium
}
shader_t;

shader_t *shader_init(FILE* fd);
int shader_shape_init(uint32_t shapeid, struct shader_so_t *self);
void shader_cleanup(shader_t *s);

// surface material or phase function callbacks:

// sample new outgoing direction past last vertex in path:
mf_t shader_sample(path_t *p);
// recompute the random numbers that would produce the outgoing direction p->e[v+1].omega.
// if not supported by the shader, returns -1.0 in all dimensions
// returns 0 if successful, 1 if the direction could not have been produced by sampling at all.
int shader_inverse_sample(const path_t *p, int v, float *r_omega_x, float *r_omega_y, float *r_scatter_mode);
mf_t  shader_brdf(path_t *p, int v);
float shader_prepare(path_t *p, int v);
mf_t  shader_pdf(const path_t *p, int v);
mf_t  shader_pdf_adj(const path_t *p, int v);


// environment map shader:

mf_t shader_sky_eval(const path_t *p, int v);
mf_t shader_sky_pdf(const path_t *p, int v);
mf_t shader_sky_sample(path_t *p);
mf_t shader_sky_pdf_next_event(const path_t *p, int v);
mf_t shader_sky_sample_next_event(path_t *p);


// related to volumes:

// sets v[0].interior and e[0].vol:
void shader_exterior_medium(path_t *p);
// returns 1 if there is a heterogeneous medium on the edge, 0 otherwise
int shader_vol_hete(const path_t *p, int e);
// compute transmittance
mf_t shader_vol_transmittance(path_t *p, int e);
// sample free distance
float shader_vol_sample(path_t *p, int e);
// pdf for sampling method
mf_t shader_vol_pdf(const path_t *p, int e);
// adjoint pdf of free path sampling
mf_t shader_vol_pdf_adjoint(const path_t *p, int e);


static inline void shader_print_info(FILE *fd)
{
  fprintf(fd, "shader   : loadable shader support\n");
#ifdef SHADER_ROUGHENING
  fprintf(fd, "           using bsdf roughening as regularisation.\n");
#endif
}

#endif
