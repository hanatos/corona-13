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
#include "spectrum.h"

#include <stdio.h>
#include <math.h>
#include <assert.h>

typedef struct diff_t
{
  float eta1;
  float eta2;
  char filename[1024];
  float ref[90*400];
}
diff_t;

//Fresnel:
//Reflexionskoeffizient fuer die senkrechte Komponente
float reflexVert(float cosinangle, float cosoutangle, float n1, float n2)
{
  float rs1 = n1*cosinangle - n2*cosoutangle;
  float rs2 = n1*cosinangle + n2*cosoutangle;

  float rc = (rs1/rs2)*(rs1/rs2);

  return rc;
}

//Reflexionskoeffizient fuer die parallele Komponente
float reflexPar(float cosinangle, float cosoutangle, float n1, float n2)
{
  float rp1 = n2*cosinangle - n1*cosoutangle;
  float rp2 = n2*cosinangle + n1*cosoutangle;

  float rc = (rp1/rp2)*(rp1/rp2);

  return rc;
}

//Reflexionskoeffizient bei unpolarisiertem Licht
float reflexUnpo(float cosinangle, float n1, float n2)
{
  float sininangle = sqrtf(1.0f - cosinangle*cosinangle);

  //Ausfallswinkel nach Snellius in Abhaengigkeit von Einfallswinkel und Brechungsindizes der Medien	

  float sinoutangle = sininangle*(n1/n2);

  //Totalreflexion? 
  if(sinoutangle >= 1.0f){return 1.0f;}	  

  float cosoutangle = sqrtf(1.0f - sinoutangle*sinoutangle);

  float av = (reflexPar(cosinangle, cosoutangle, n1, n2) +
      reflexVert(cosinangle, cosoutangle, n1, n2))/2;

  return av;
}

extern float specularity(const float *in, const rayhit_t *hit, const float rr, void *data)
{
  return 1.0f;
}

extern int init(FILE *f, void **data)
{
  diff_t *d = (diff_t *)malloc(sizeof(diff_t));
  *data = d;
  int i = fscanf(f, "%f %f %s", &d->eta1, &d->eta2, d->filename);
  if(i < 3)
  {
    fprintf(stderr, "[diffractive] could not parse all arguments! expecting: <eta(400nm)> <eta(700nm)> <file>\n");
    d->eta1 = d->eta2 = 1.5f;
    return 1;
  }
  int dreggn = 0;
  dreggn = fscanf(f, "%*[^\n]\n");

  return 0;
}

extern void cleanup(void *data)
{
  free(data);
}

int shape_init(uint32_t shapeid, shader_so_t *self)
{
  diff_t *diff = (diff_t *)self->data;
  shader_data_t *data = shader_data(diff->filename, SDAT_RAW);
  if(!data) 
  {
    fprintf(stderr, "[diffractive] could not open %s!\n", diff->filename);
    return 1;
  }
  char *d = (char *)data->data;
  for(int j=0;j<90;j++) for(int i=0;i<400;i++)
    // sscanf(d, "%f", diff->ref + 400*j + i);
    diff->ref[400*j+i] = strtof(d, &d);
  for(int j=0;j<90;j++) for(int i=0;i<400;i++) diff->ref[400*j+i] *= 1./100.;
  /*for(int j=0;j<90;j++)
  {
    for(int i=0;i<400;i++) printf("%.1f ", diff->ref[400*j+i]);
    printf("\n\n");
  }*/
  return 0;
}

extern float sample(const float* omega_in, rayhit_t* hit, float* omega_out, const float qx, const float qy, const float rr, void* data)
{
  float eta;
  float cos_theta_i = omega_in[0]*hit->normal[0] + omega_in[1]*hit->normal[1] + omega_in[2]*hit->normal[2];

  // TODO: use fraunhofer number nd nf nc (for 589.2 nm, 486.1 nm and 656.3 nm) (=> abbe number)
  diff_t *diff = (diff_t *)data;
  float eta1 = diff->eta1;
  float eta2 = diff->eta2;
  hit->ior = eta = (700.0 - hit->lambda)*(700.0 - hit->lambda)*(eta1-eta2)/((700.0-400.0)*(700.0-400.0)) + eta2;
  // cauchy:
  // eta = eta0 + b /(lambda*lambda) //lambda in um or b in nm^2

  if(!hit->inside) eta = 1.0f/eta;

  const float cos_theta_t_2 = 1 - eta*eta * (1 - cos_theta_i*cos_theta_i);

#if 0
  // schlick's approximation of fresnel
  const float r_s = (eta - 1)*(eta - 1)/((eta + 1)*(eta + 1));
  const float R = r_s + (1 - r_s)*powf(1 + cos_theta_i, 5);
  // real fresnel for unpolarized light:
  // float R;
  // if(hit->inside) R = reflexUnpo(- cos_theta_i, hit->ior, 1.0);
  // else R = reflexUnpo(- cos_theta_i, 1.0, hit->ior);
#else
  const float angle = fminf(88.0, fmaxf(0.0, (180.0/M_PI*acosf(-cos_theta_i))));
  const float lambda = hit->lambda - 400.0;
  const float R00 = diff->ref[(int)(400*(int)(angle)+(int)lambda)];
  const float R01 = diff->ref[(int)(400*(int)(angle)+(int)lambda+1)];
  const float R10 = diff->ref[(int)(400*(int)(angle+1)+(int)lambda)];
  const float R11 = diff->ref[(int)(400*(int)(angle+1)+(int)lambda+1)];
  const float f1 = hit->lambda - (int)hit->lambda;
  const float f2 = angle - (int)angle;
  const float R = (1.0-f1)*(1.0-f2)*R00 + f1*(1.0-f2)*R01 + (1.0-f1)*f2*R10 + f1*f2*R11;
  // printf("R( %f %f ) = %f\n", (180.0/M_PI*acosf(-cos_theta_i)), lambda, R);
#endif

  // fresnel || total inner reflection -> reflection
  if((rr < R) || (cos_theta_t_2 <= 0.0f))
    for(int k=0;k<3;k++) omega_out[k] = omega_in[k] - 2*cos_theta_i*hit->normal[k];
  else
  {
    const float f = (- sqrtf(cos_theta_t_2) - eta*cos_theta_i);
    for(int k=0;k<3;k++) omega_out[k] = omega_in[k]*eta + f * hit->normal[k];
    const float il = 1.0f/sqrtf(omega_out[0]*omega_out[0] + omega_out[1]*omega_out[1] + omega_out[2]*omega_out[2]);
    for(int k=0;k<3;k++) omega_out[k] *= il;
  }
  // FIXME: eta1/eta2 for light rays!
  return .9f;
}


float pdf_rr(const float *omega_in, const rayhit_t *hit, const float *omega_out, const float rr, void *data)
{
  return 0.0f;
}

float pdf(const float *omega_in, const rayhit_t *hit, const float *omega_out, void *data)
{ // specularity will take care of dirac part.
  return 1.0f;
}

extern float brdf(float* in, rayhit_t* hit, float* out, void *data)
{
  return 0.0f;//1.0/M_PI;
}

