// #include "corona_common.h"
#include "shader.h"
#include "spectrum_common.h"
#include "sampler_common.h"

// g++ -Wall -lm -lc -pipe  fluor_rgb.cpp multispecBRDF.C -I. -I../../include -I../../build -shared -o libfluor_rgb.so -fPIC
#include "multispecBRDF.hh"

extern "C"
{

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
double drand48();

typedef struct fluor_t
{
  multispecBRDF *s;
  char filename[512];
}
fluor_t;

#define ENERGY_SCALE (1./(300.0*20.0))
#define PREBAKED_COS(X) (1.f/(X))
// #define DIFFUSE

extern float specularity(const float *omega_in, const rayhit_t *hit, const float rr, void *data)
{
  return .4f;
}

extern float pdf(const float *omega_in, const rayhit_t *hit, const float *omega_out, void *data)
{
  return 1.0f/M_PI;
}

extern float pdf_rr(const float *omega_in, const rayhit_t *hit, const float *omega_out, const float rr, void *data)
{
  return dotproduct(hit->normal, omega_out)/M_PI;
}

float get_sum(const float theta_in, const float theta_out, const float phi_rel, multispecBRDF *s, const float lambda, int adjoint)
{
  // FIXME:
  adjoint = 1;
  // calculate color appearing towards camera (lambda_out spectrum)
  // when illuminated by constant (~white) spectrum (lambda_in spectrum, p(lambda_in) = 1.0)
  // float rgbout[3];
  // for(int k=0;k<3;k++) rgb[k] = 0.0f;
  float sum = 0.0f;
  int N = 15;
  // int M = 15;
  // sample lambda_out
  // for(int i=0;i<M;i++)
  {
    // and get rgb weighting factors
    // const float lambda_out = 400.0f + 20.0f*i;
    // const float lambda_out = 400.0f + (700.0-400.0)*drand48();// i/(float)M;
    // const float lambda_out = 400.0f + (700.0-400.0)*(i+.5f)/(float)M;
    // spectrum_p_to_rgb(lambda_out, 1.0f, rgbout);
    // printf("lambda_out %f, rgb %f %f %f\n", lambda_out, rgbout[0], rgbout[1], rgbout[2]);
    // for(int k=0;k<=i;k++)
    for(int k=0;k<N;k++)
    {
      // const float lambda_in = 400.0f + 20.0f*k;
      // const float lambda_in = 400.0f + (700.0-400.0)*drand48();// k/(float)N;
      const float lambda_in = 400.0f + (700.0-400.0)*(k+.5f)/(float)N;
      // spectrum_p_to_rgb(lambda_in, 1.0f, rgbin);
      float rd;
      if(adjoint) rd = s->readout(theta_out, theta_in,  phi_rel, lambda_in, lambda) *PREBAKED_COS(cosf(theta_out));
      else        rd = s->readout(theta_in,  theta_out, phi_rel, lambda_in, lambda)*PREBAKED_COS(cosf(theta_in));
      const float f = 300.0f*300.0f*rd*ENERGY_SCALE/((float)N);
      // printf("lambda_in, p(%f) = %f\n", lambda_in, f);
      // for(int n=0;n<3;n++) rgb[n] += rgbout[n]*f;
      sum += f;
    }
  }
  return sum;
  // for(int n=0;n<3;n++) if(rgb[n]<0.0f) rgb[n] = 0.0f;
#if 0
  printf("mean rgb appearance: %f %f %f\n", rgb[0], rgb[1], rgb[2]);

  N = 300;
  M = 700;
  // spectral version:
  for(int k=0;k<3;k++) rgb[k] = 0.0f;
  // sample lambda_out (ray from pt, for example)
  for(int i=0;i<M;i++)
  {
    // const float lambda_out = 400.0f + 20.0f*i;
    // const float lambda_out = spectrum_sample_lambda(drand48(), NULL);
    const float lambda_out = 400.f + (700.0f - 400.0f)*drand48();
    // const float lambda_out = 400.0f + (700.0-400.0)*i/(float)M;
    float sum = 0.0f;
    // sample (fewer) lambda_in from light source

#if 1
    // correct for loop
    for(int k=0;k<N;k++)
    {
      const float lambda_in = 400.f + (700.0f - 400.0f)*drand48();
      // const float lambda_in = 400.0f + (700.0-400.0)*k/(float)N;
      //                 p(lambda_out)   p(lambda_in)
      const float corr = 300.0f*(1./(float)M) * 300.0/(float)N;
      if(adjoint) rd = s->readout(theta_out, theta_in,  phi_rel, lambda_in, lambda_out) *PREBAKED_COS(cosf(theta_out));
      else        rd = s->readout(theta_in,  theta_out, phi_rel, lambda_in, lambda_out)*PREBAKED_COS(cosf(theta_in));
      const float f = corr * rd * ENERGY_SCALE;
      sum += f;
    }
#else
    // more efficient, incorrect sampling :(
    for(int k=0;k<N;k++)
    {
      const float lambda_in = 400.f + (lambda_out - 400.0f)*drand48();
      // const float lambda_in = 400.0f + (700.0-400.0)*k/(float)N;
      //                 p(lambda_out)   p(lambda_in)
      const float corr = 300.0f*(1./(float)M) * (lambda_out - 400.0)/(float)N;
      const float f = corr * s->readout(theta_out, theta_in, phi_rel, lambda_in, lambda_out)*ENERGY_SCALE*PREBAKED_COS(cosf(theta_out));
      sum += f;
    }
#endif
    // ray terminated, sum up in accum buf.
    spectrum_p_to_rgb(lambda_out, sum, rgbout);
    for(int n=0;n<3;n++) rgb[n] += rgbout[n];
  }
  printf("mean spec appearance: %f %f %f\n", rgb[0], rgb[1], rgb[2]);
  exit(1);
#endif
}

#if 0
void get_rgb(const float theta_in, const float theta_out, const float phi_rel, multispecBRDF *s, float *rgb)
{
  rgb[0] = rgb[1] = rgb[2] = 0.0f;
  for(int i=0;i<15;i++)
  {
    const float lambda_out = 400.0f + 20.0f*i;
    for(int k=0;k<=i;k++)
    {
      float tmprgb[3];
      spectrum_p_to_rgb(lambda_out, s->readout(theta_out, theta_in, phi_rel, 400.0f + k*20.0f, lambda_out)*PREBAKED_COS(cosf(theta_out))/15.0, tmprgb);
      for(int j=0;j<3;j++) rgb[j] += tmprgb[j];
    }
  }
}
#endif

extern float sample(const float *omega_in, rayhit_t *hit, float *omega_out, const float x1, const float x2, const float rr, void *data)
{
  float x, y, z;
  sample_cos(&x, &y, &z, x1, x2);
  for(int k=0;k<3;k++) omega_out[k] = hit->a[k]*x + hit->b[k]*y + hit->normal[k]*z;
  fluor_t *f = (fluor_t *)data;
  const float theta_in  = acosf(fminf(.999f, fmaxf(0.0f, - dotproduct(hit->normal, omega_in))));
  const float theta_out = acosf(z);
  float phi_rel = fabsf(2.f*M_PI*x2 - atan2f(- dotproduct(hit->b, omega_in), - dotproduct(hit->a, omega_in)));
  if(phi_rel > M_PI) phi_rel = 2.0f*M_PI - phi_rel;
  return M_PI*get_sum(theta_in, theta_out, phi_rel, f->s, hit->lambda, hit->adjoint);
}

extern float brdf(float *omega_in, rayhit_t *hit, float *omega_out, void *data)
{
  // only works with ptdl!
  const float dot_in = - dotproduct(hit->normal, omega_in);
  const float dot_out = dotproduct(hit->normal, omega_out);
  if(dot_in < 0.001f || dot_out < 0.001f) return 0.0f;
  fluor_t *f = (fluor_t *)data;
  const float theta_in  = acosf(fminf(.999f, dot_in));
  const float theta_out = acosf(fminf(.999f, dot_out));
  float phi_rel = fabsf(atan2f(dotproduct(hit->b, omega_out), dotproduct(hit->a, omega_out)) - atan2f(- dotproduct(hit->b, omega_in), - dotproduct(hit->a, omega_in)));
  if(phi_rel > M_PI) phi_rel = 2.0f*M_PI - phi_rel;
  return get_sum(theta_in, theta_out, phi_rel, f->s, hit->lambda, hit->adjoint);
}

extern int init(FILE *f, void **data)
{
  fluor_t *s = (fluor_t *)malloc(sizeof(fluor_t));
  *data = s;
  if(fscanf(f, "%s", s->filename) != 1)
  {
    printf("[fluor_rgb::init] ERROR: could not parse arguments!\n");
    printf(" expected: <data filename>\n");
    return 1;
  }
  // f->s = new multispecBRDF(string(filename));
  fscanf(f, "%*[^\n]\n");
  return 0;
}

int uvn_init(const char *fname, unsigned int tri, shader_so_t *self)
{
  fluor_t *f = (fluor_t *)self->data;
  shader_data_t *data = shader_data(f->filename, SDAT_RAW);
  if(!data)
  {
    fprintf(stderr, "[fluor_rgb]: could not open file `%s' !\n", f->filename);
    return 1;
  }

  f->s = new multispecBRDF((char *)data->data);
  f->s->smoothSpecularHighlight();
  return 0;
}

/*extern void cleanup(void *data)
{
  multispecBRDF *s = (multispecBRDF *)data;
  delete s;
}*/

}

