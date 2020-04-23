#include "corona_common.h"
#include "shader.h"
#include "spectrum_common.h"

// g++ -Wall -lm -lc -pipe fluor.cpp multispecBRDF.C -I. -I../../include -I../../build -shared -o libfluor.so -fPIC
#include "multispecBRDF.hh"

extern "C"
{

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

typedef struct fluor_t
{
  multispecBRDF *s;
  char filename[512];
}
fluor_t;


#define ENERGY_SCALE (1.0/20.0)  // no need to compensate for pdf(lambda_out), the render kernel will do this.
#define PREBAKED_COS(X) (1.f/(X))
//#define PREBAKED_COS(X) 1.f

extern float specularity(const float *omega_in, const rayhit_t *hit, const float rr, void *data)
{
  return .4f;
}

extern float prepare(const ray_t *ray, rayhit_t *hit, const float rr, void *data)
{ // sample new lambda (longer wavelength towards eye) and compensate pdf.
  hit->lambda_in = hit->lambda;

  float pdf;
  hit->lambda = spectrum_sample_lambda(rr, &pdf);
  return 1./pdf;
#if 0
  // :(
  if(hit->adjoint)
  {
    hit->lambda = 400.0f + (hit->lambda_in-400.0f)*rr;
    return (hit->lambda_in-400.0f);
  }
  else
  {
    hit->lambda = hit->lambda_in + (700.0f - hit->lambda_in)*rr;
    return (700.0f - hit->lambda_in);
  }
#endif
}


extern float pdf(const float *omega_in, const rayhit_t *hit, const float *omega_out, void *data)
{
  return 1.0f/((700.0f-400.0f)*M_PI);
}

extern float pdf_rr(const float *omega_in, const rayhit_t *hit, const float *omega_out, const float rr, void *data)
{
  return dotproduct(hit->normal, omega_out)/M_PI;
}

extern float sample_adj(const float *omega_in, rayhit_t *hit, float *omega_out, const float x1, const float x2, const float rr, void *data)
{
  float x, y, z;
  // multispecBRDF *s = (multispecBRDF *)data;
  fluor_t *f = (fluor_t *)data;
  const float su = sqrtf(x1);
  x = su*cosf(2.f*M_PI*x2);
  y = su*sinf(2.f*M_PI*x2);
  z = sqrtf(1.0 - x1);
  for(int k=0;k<3;k++) omega_out[k] = hit->a[k]*x + hit->b[k]*y + hit->normal[k]*z;
  const float theta_in  = acosf(fminf(.999f, fmaxf(0.0f, - dotproduct(hit->normal, omega_in))));
  const float theta_out = acosf(z);
  float phi_rel = fabsf(2.f*M_PI*x2 - atan2f(- dotproduct(hit->b, omega_in), - dotproduct(hit->a, omega_in)));
  if(phi_rel > M_PI) phi_rel = 2.0f*M_PI - phi_rel;
  return M_PI*f->s->readout(theta_out, theta_in, phi_rel, hit->lambda, hit->lambda_in)*ENERGY_SCALE*PREBAKED_COS(dotproduct(hit->normal, omega_out));
}

extern float sample(const float *omega_in, rayhit_t *hit, float *omega_out, const float x1, const float x2, const float rr, void *data)
{
  float x, y, z;
  // multispecBRDF *s = (multispecBRDF *)data;
  fluor_t *f = (fluor_t *)data;
  const float su = sqrtf(x1);
  x = su*cosf(2.f*M_PI*x2);
  y = su*sinf(2.f*M_PI*x2);
  z = sqrtf(1.0 - x1);
  for(int k=0;k<3;k++) omega_out[k] = hit->a[k]*x + hit->b[k]*y + hit->normal[k]*z;
  const float theta_in  = acosf(fminf(.999f, fmaxf(0.0f, - dotproduct(hit->normal, omega_in))));
  const float theta_out = acosf(z);
  float phi_rel = fabsf(2.f*M_PI*x2 - atan2f(- dotproduct(hit->b, omega_in), - dotproduct(hit->a, omega_in)));
  if(phi_rel > M_PI) phi_rel = 2.0f*M_PI - phi_rel;
  return M_PI*f->s->readout(theta_in, theta_out, phi_rel, hit->lambda_in, hit->lambda)*ENERGY_SCALE*PREBAKED_COS(z);
}

extern float brdf(float *omega_in, rayhit_t *hit, float *omega_out, void *data)
{
  // multispecBRDF *s = (multispecBRDF *)data;
  fluor_t *f = (fluor_t *)data;
  const float dot_in = - dotproduct(hit->normal, omega_in);
  const float dot_out = dotproduct(hit->normal, omega_out);
  if(dot_in < 0.001f || dot_out < 0.001f) return 0.0f;
  const float theta_in  = acosf(fminf(.999f, dot_in));
  const float theta_out = acosf(fminf(.999f, dot_out));
  float phi_rel = fabsf(atan2f(dotproduct(hit->b, omega_out), dotproduct(hit->a, omega_out)) - atan2f(- dotproduct(hit->b, omega_in), - dotproduct(hit->a, omega_in)));
  if(phi_rel > M_PI) phi_rel = 2.0f*M_PI - phi_rel;
  if(hit->adjoint) return f->s->readout(theta_out, theta_in, phi_rel, hit->lambda, hit->lambda_in)*ENERGY_SCALE*PREBAKED_COS(dot_out);
  else             return f->s->readout(theta_in, theta_out, phi_rel, hit->lambda_in, hit->lambda)*ENERGY_SCALE*PREBAKED_COS(dot_out);
}

extern int init(FILE *f, void **data)
{
  fluor_t *s = (fluor_t *)malloc(sizeof(fluor_t));
  *data = s;
  if(fscanf(f, "%s", s->filename) != 1)
  {
    printf("[fluor::init] ERROR: could not parse arguments!\n");
    printf(" expected: <data filename>\n");
    return 1;
  }
  // multispecBRDF *s = new multispecBRDF(string(filename));

#ifdef REDUCE
  s->reduce();
  s->reduce();
#endif

#if 0
  for(float lambda_out=400.0f; lambda_out < 700.0f; lambda_out+=10.0f)
  {
    for(float lambda_in=400.0f; lambda_in < 700.0f; lambda_in+=10.0f)
    {
      printf("%.2f ", s->readout(1.0, 1.0, M_PI, lambda_in, lambda_out));
    }
    printf("\n");
  }
#endif

#if 0
  float max = 0.0f;
  for(int li=0;li<15;li++)
  {
    for(int oi=0;oi<10;oi++)
    {
      float sum = 0.0f;
      float z = sqrtf(1.0 - drand48());
      const float theta_out = acosf(z);
      for(int oo=0;oo<100;oo++)
      {
        // ~ cos
        // z = sqrtf(1.0 - drand48());
        // ~ const
        z = 1.0 - drand48();
        const float theta_in = acosf(z);
        const float phi_rel = M_PI*drand48();
        for(int lo=0;lo<15;lo++)
          sum += (700.0f - 400.0f)*2.0f*M_PI*s->readout(theta_in, theta_out, phi_rel, 400.0f + li*20.0f, 400.0f + lo*20.0f)/(100.0*15.0);//*ENERGY_SCALE;//(32.f*pdf);
      }
      max = fmaxf(sum, max);
      //printf("[fluor] lo, oo = [%d %d], sum = %f\n", lo, oo, sum);
    }
  }
  printf("[fluor] max transport operator norm: %f\n", max);
#endif

  // *data = s;
  fscanf(f, "%*[^\n]\n");
  return 0;
}

extern void cleanup(void *data)
{
  multispecBRDF *s = (multispecBRDF *)data;
  delete s;
}

int uvn_init(const char *fname, unsigned int tri, shader_so_t *self)
{
  fluor_t *f = (fluor_t *)self->data;
  shader_data_t *data = shader_data(f->filename, SDAT_RAW);
  if(!data)
  {
    fprintf(stderr, "[fluor]: could not open file `%s' !\n", f->filename);
    return 1;
  }

  f->s = new multispecBRDF((char *)data->data);
  f->s->smoothSpecularHighlight();

#if 1 // print one rerad matrix
  // const float z1 = 1.0;
  const float theta_out = 0.0;//acosf(z1);
  const float z = 0.7;
  const float theta_in = acosf(z);
  const float phi_rel = M_PI;
  for(int lo=0;lo<15;lo++)
  {
    for(int li=0;li<15;li++)
    {
      const double fi = f->s->readout(theta_in, theta_out, phi_rel, 400.0f + li*20.0f, 400.0f + lo*20.0f)*PREBAKED_COS(z);
      printf("%f, ", fi);
    }
    printf("\n");
  }
#endif
#if 0 // test reciprocity
  double sum = 0.0f, sumfi = 0.0, sumfo = 0.0;
  float max = - INFINITY;
  float min = INFINITY;
  int num = 0;
  for(int li=0;li<15;li++) for(int lo=0;lo<15;lo++)
  // const int li = 7, lo = 7;
  {
    if(li == lo) continue;
    for(int oi=0;oi<100;oi++)
    {
      // float z = sqrtf(1.0 - drand48());
      float cos_out = 1.0 - drand48(); // ~ const
      cos_out = fmaxf(0.1, cos_out);
      const float theta_out = acosf(cos_out);
      for(int oo=0;oo<100;oo++)
      {
        // z = sqrtf(1.0 - drand48()); // ~ cos
        float z = 1.0 - drand48(); // ~ const
        z = fmaxf(0.1, z);
        z = fminf(0.99, cos_out + drand48()*0.1);
        const float theta_in = acosf(z);
        const float phi_rel = M_PI*drand48();
        const double fi = f->s->readout(theta_in, theta_out, phi_rel, 400.0f + li*20.0f, 400.0f + lo*20.0f)*PREBAKED_COS(z);
        const double fo = f->s->readout(theta_out, theta_in, phi_rel, 400.0f + li*20.0f, 400.0f + lo*20.0f)*PREBAKED_COS(cos_out);
        num ++;
        const double tmax = fmax(fi, fo);
        const double dev = fabs(fi-fo)/tmax;
        // const double dev = (fi-fo)/tmax;
        if(dev < min) min = dev;
        if(dev > max) max = dev;
        if(dev > 0) sum += dev;
        sumfi += (fi-fo)*(fi-fo);
        sumfo += .5f*(fo+fi);
        // printf("dev: %f %f %f\n", dev, fi, fo);
      }
    }
  }
  printf("recip within %f%% [%f -- %f]\n", 100*sum/num, min, max);
  printf("recp2: %f\n", sqrtf(sumfi/sumfo));
  exit(1);
#endif
  return 0;
}


}

