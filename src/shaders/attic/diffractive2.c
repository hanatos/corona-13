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
#include <strings.h>

typedef struct diff_ray_t
{
  float theta, phi, intensity;
}
diff_ray_t;

typedef struct diff_t
{
  float eta1;
  float eta2;
  char ref_filename[1024];
  char tra_filename[1024];
  int num_theta, num_phi;
  int *ref_dray_cnt;
  diff_ray_t **ref_dray;
  int *tra_dray_cnt;
  diff_ray_t **tra_dray;
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
  int i = fscanf(f, "%f %f %s %s", &d->eta1, &d->eta2, d->ref_filename, d->tra_filename);
  if(i < 3)
  {
    fprintf(stderr, "[diffractive] could not parse all arguments! expecting: <eta(400nm)> <eta(700nm)> <ref file> <trans file>\n");
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
  shader_data_t *rdata = shader_data(diff->ref_filename, SDAT_RAW);
  shader_data_t *tdata = shader_data(diff->tra_filename, SDAT_RAW);
  if(!rdata || !tdata)
  {
    fprintf(stderr, "[diffractive] could not open `%s' and `%s'!\n", diff->ref_filename, diff->tra_filename);
    return 1;
  }
  int32_t *d = (int32_t *)rdata->data;
  const int buflen = rdata->filesize/sizeof(float);
  diff->num_phi = d[0];
  diff->num_theta = d[1];
  diff->ref_dray_cnt = (int *)malloc(sizeof(int)*diff->num_theta*diff->num_phi);
  bzero(diff->ref_dray_cnt, sizeof(int)*diff->num_theta*diff->num_phi);
  diff->tra_dray_cnt = (int *)malloc(sizeof(int)*diff->num_theta*diff->num_phi);
  bzero(diff->tra_dray_cnt, sizeof(int)*diff->num_theta*diff->num_phi);

  d += 2;
  diff->ref_dray = (diff_ray_t **)malloc(sizeof(diff_ray_t*)*diff->num_theta*diff->num_phi);
  for(int i=0;i<diff->num_theta*diff->num_phi;i++)
  {
    assert(i-2 < buflen);
    diff->ref_dray_cnt[i] = d[i];
    diff->ref_dray[i] = (diff_ray_t *)malloc(sizeof(diff_ray_t)*diff->ref_dray_cnt[i]);
  }
  d = (int32_t *)tdata->data + 2;
  diff->tra_dray = (diff_ray_t **)malloc(sizeof(diff_ray_t*)*diff->num_theta*diff->num_phi);
  for(int i=0;i<diff->num_theta*diff->num_phi;i++)
  {
    assert(i-2 < buflen);
    diff->tra_dray_cnt[i] = d[i];
    diff->tra_dray[i] = (diff_ray_t *)malloc(sizeof(diff_ray_t)*diff->tra_dray_cnt[i]);
  }

  float *f = ((float *)rdata->data) + 2 + diff->num_theta*diff->num_phi;
  for(int i=0;i<diff->num_theta*diff->num_phi;i++)
  {
    assert(i-2-diff->num_theta*diff->num_phi < buflen);
    for(int k=0;k<diff->ref_dray_cnt[i];k++)
    {
      diff->ref_dray[i][k].phi       = M_PI/180.0 * f[3*i];
      diff->ref_dray[i][k].theta     = cosf(M_PI/180.0*f[3*i+1]);
      diff->ref_dray[i][k].intensity = f[3*i+2];
    }
  }
  f = ((float *)tdata->data) + 2 + diff->num_theta*diff->num_phi;
  for(int i=0;i<diff->num_theta*diff->num_phi;i++)
  {
    assert(i-2-diff->num_theta*diff->num_phi < buflen);
    for(int k=0;k<diff->tra_dray_cnt[i];k++)
    {
      diff->tra_dray[i][k].phi       = M_PI/180.0 * f[3*i];
      diff->tra_dray[i][k].theta     = cosf(M_PI/180.0*f[3*i+1]);
      diff->tra_dray[i][k].intensity = f[3*i+2];
    }
  }
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



#if 0
  if(!hit->inside) eta = 1.0f/eta;
  const float cos_theta_t_2 = 1 - eta*eta * (1 - cos_theta_i*cos_theta_i);
  // schlick's approximation of fresnel
  const float r_s = (eta - 1)*(eta - 1)/((eta + 1)*(eta + 1));
  const float R = r_s + (1 - r_s)*powf(1 + cos_theta_i, 5);
  // real fresnel for unpolarized light:
  // float R;
  // if(hit->inside) R = reflexUnpo(- cos_theta_i, hit->ior, 1.0);
  // else R = reflexUnpo(- cos_theta_i, 1.0, hit->ior);
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
#else
  const int theta = fminf(89, fmaxf(0.0, (180.0/M_PI*acosf(-cos_theta_i))));
  const int phi   = fminf(189, 180.0/M_PI * atan2f(dotproduct(omega_in, hit->b), dotproduct(omega_in, hit->a)));
  const int phi2  = abs(phi);
  int   pos   = diff->num_theta*phi2 + theta;
  if(pos < 0) pos = 0;
  if(pos >= diff->num_theta*diff->num_phi) pos = diff->num_theta*diff->num_phi;
  // TODO: rr based on ref => goto transmission
  if(diff->ref_dray_cnt[pos] == 0) return 0.0f;
  const int   d     = rr * diff->ref_dray_cnt[pos];
  const float R     = diff->ref_dray[pos][d].intensity;
  const float theta_out = diff->ref_dray[pos][d].theta;
  const float phi_out   = diff->ref_dray[pos][d].phi;
  const float phi_out2  = phi < 0 ? -phi_out : phi_out;
  for(int k=0;k<3;k++) omega_out[k] = hit->normal[k] * theta_out + sqrtf(1.0-theta_out*theta_out)*(hit->a[k]*cosf(phi_out2) + hit->b[k]*sinf(phi_out2));
  // for(int k=0;k<3;k++) omega_out[k] = hit->normal[k] * cosf(M_PI/180.0*theta) + sinf(M_PI/180.0*theta)*(hit->a[k]*cosf(M_PI/180.0*phi) + hit->b[k]*sinf(M_PI/180.0*phi));
  // for(int k=0;k<3;k++) omega_out[k] = hit->normal[k] * theta_out + sqrtf(1.0-theta_out*theta_out)*(hit->a[k]*cosf(M_PI/180.0*phi) + hit->b[k]*sinf(M_PI/180.0*phi));
  // printf("R( %f %f ) = %f\n", (180.0/M_PI*acosf(-cos_theta_i)), lambda, R);
  return R;
#endif
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

