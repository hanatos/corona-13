#include "corona_common.h"
#include "shader.h"
#include "spectrum_common.h"
#include <stdlib.h>
double drand48();

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


#define ENERGY_SCALE (1./(300.0*20.0))  // no need to compensate for pdf(lambda_out), the render kernel will do this.
#define PREBAKED_COS(X) (1.f/(X))
//#define PREBAKED_COS(X) 1.f


// srgb dreggn
const int rgb_lambda_num = 89;
const float rgb_lambda_min = 380.0f;
const float rgb_lambda_step = 5.0f;

const float rgb_red[] = {
1.5000E-03, 3.8000E-03, 8.9000E-03, 1.8800E-02, 3.5000E-02, 5.3100E-02, 7.0200E-02, 7.6300E-02, 7.4500E-02, 5.6100E-02, 3.2300E-02, 4.4000E-03, 4.7800E-02, 9.7000E-02, 1.5860E-01, 2.2350E-01, 2.8480E-01, 3.3460E-01, 3.7760E-01, 4.1360E-01, 4.3170E-01, 4.4520E-01, 4.3500E-01, 4.1400E-01, 3.6730E-01, 2.8450E-01, 1.8550E-01, 4.3500E-02, 1.2700E-01, 3.1290E-01, 5.3620E-01, 7.7220E-01, 1.0059E+00, 1.2710E+00, 1.5574E+00, 1.8465E+00, 2.1511E+00, 2.4250E+00, 2.6574E+00, 2.9151E+00, 3.0779E+00, 3.1613E+00, 3.1673E+00, 3.1048E+00, 2.9462E+00, 2.7194E+00, 2.4526E+00, 2.1700E+00, 1.8358E+00, 1.5179E+00, 1.2428E+00, 1.0070E+00, 7.8270E-01, 5.9340E-01, 4.4420E-01, 3.2830E-01, 2.3940E-01, 1.7220E-01, 1.2210E-01, 8.5300E-02, 5.8600E-02, 4.0800E-02, 2.8400E-02, 1.9700E-02, 1.3500E-02, 9.2400E-03, 6.3800E-03, 4.4100E-03, 3.0700E-03, 2.1400E-03, 1.4900E-03, 1.0500E-03, 7.3900E-04, 5.2300E-04, 3.7200E-04, 2.6500E-04, 1.9000E-04, 1.3600E-04, 9.8400E-05, 7.1300E-05, 5.1800E-05, 3.7700E-05, 2.7600E-05, 2.0300E-05, 1.4900E-05, 1.1000E-05, 8.1800E-06, 6.0900E-06, 4.5500E-06};
const float rgb_green[] = {
-4.0000E-04, -1.0000E-03, -2.5000E-03, -5.9000E-03, -1.1900E-02, -2.0100E-02, -2.8900E-02, -3.3800E-02, -3.4900E-02, -2.7600E-02, -1.6900E-02, 2.4000E-03, 2.8300E-02, 6.3600E-02, 1.0820E-01, 1.6170E-01, 2.2010E-01, 2.7960E-01, 3.4280E-01, 4.0860E-01, 4.7160E-01, 5.4910E-01, 6.2600E-01, 7.0970E-01, 7.9350E-01, 8.7150E-01, 9.4770E-01, 9.9450E-01, 1.0203E+00, 1.0375E+00, 1.0517E+00, 1.0390E+00, 1.0029E+00, 9.6980E-01, 9.1620E-01, 8.5710E-01, 7.8230E-01, 6.9530E-01, 5.9660E-01, 5.0630E-01, 4.2030E-01, 3.3600E-01, 2.5910E-01, 1.9170E-01, 1.3670E-01, 9.3800E-02, 6.1100E-02, 3.7100E-02, 2.1500E-02, 1.1200E-02, 4.4000E-03, 7.8000E-05, -1.3680E-03, -1.9880E-03, -2.1680E-03, -2.0060E-03, -1.6420E-03, -1.2720E-03, -9.4700E-04, -6.8300E-04, -4.7800E-04, -3.3700E-04, -2.3500E-04, -1.6300E-04, -1.1100E-04, -7.4800E-05, -5.0800E-05, -3.4400E-05, -2.3400E-05, -1.5900E-05, -1.0700E-05, -7.2300E-06, -4.8700E-06, -3.2900E-06, -2.2200E-06, -1.5000E-06, -1.0200E-06, -6.8800E-07, -4.6500E-07, -3.1200E-07, -2.0800E-07, -1.3700E-07, -8.8000E-08, -5.5300E-08, -3.3600E-08, -1.9600E-08, -1.0900E-08, -5.7000E-09, -2.7700E-09};
const float rgb_blue[] = {
6.2000E-03 , 1.6100E-02 , 4.0000E-02 , 9.0600E-02 , 1.8020E-01 , 3.0880E-01 , 4.6700E-01 , 6.1520E-01 , 7.6380E-01 , 8.7780E-01 , 9.7550E-01 , 1.0019E+00 , 9.9960E-01 , 9.1390E-01 , 8.2970E-01 , 7.4170E-01 , 6.1340E-01 , 4.7200E-01 , 3.4950E-01 , 2.5640E-01 , 1.8190E-01 , 1.3070E-01 , 9.1000E-02 , 5.8000E-02 , 3.5700E-02 , 2.0000E-02 , 9.5000E-03 , 7.0000E-04 , 4.3000E-03 , 6.4000E-03 , 8.2000E-03 , 9.4000E-03 , 9.7000E-03 , 9.7000E-03 , 9.3000E-03 , 8.7000E-03 , 8.0000E-03 , 7.3000E-03 , 6.3000E-03 , 5.3700E-03 , 4.4500E-03 , 3.5700E-03 , 2.7700E-03 , 2.0800E-03 , 1.5000E-03 , 1.0300E-03 , 6.8000E-04 , 4.4200E-04 , 2.7200E-04 , 1.4100E-04 , 5.4900E-05 , 2.2000E-06 , 2.3700E-05 , 2.8600E-05 , 2.6100E-05 , 2.2500E-05 , 1.8200E-05 , 1.3900E-05 , 1.0300E-05 , 7.3800E-06 , 5.2200E-06 , 3.6700E-06 , 2.5600E-06 , 1.7600E-06 , 1.2000E-06 , 8.1700E-07 , 5.5500E-07 , 3.7500E-07 , 2.5400E-07 , 1.7100E-07 , 1.1600E-07 , 7.8500E-08 , 5.3100E-08 , 3.6000E-08 , 2.4400E-08 , 1.6500E-08 , 1.1200E-08 , 7.5300E-09 , 5.0700E-09 , 3.4000E-09 , 2.2700E-09 , 1.5000E-09 , 9.8600E-10 , 6.3900E-10 , 4.0700E-10 , 2.5300E-10 , 1.5200E-10 , 8.6400E-11 , 4.4200E-11};

static inline float get_spectrum(const float *s, float lambda)
{
  int l = (lambda - rgb_lambda_min)/rgb_lambda_step;
  if(l < 0 || l >= rgb_lambda_num) return 0.0f;
  return s[l];
}

// box basis rgb
void get_rgb(const float lambda, float *rgb)
{
  rgb[0] = rgb[1] = rgb[2] = 0.0f;
  const int i = (int)fmaxf(0.0f, fminf(2.99f, (lambda - 400.0f)/100.0f));
  rgb[i] = 1.;///100.;
}
float eval_rgb(const float lambda, const float *rgb)
{
  const int i = (int)fmaxf(0.0f, fminf(2.99f, (lambda - 400.0f)/100.0f));
  return rgb[i];
}

void get_rgbxrgb(const float theta_in, const float theta_out, const float phi_rel, multispecBRDF *s, float *rgb, int adjoint)
{
  // FIXME:
  adjoint = 1;
  float rgbin[3], rgbout[3];
  // const float rgbint[9] = //{175.208603, 27.509781, 0.287169 };
  // {
  //   // 136.435837, 21.422018, 0.223619,
  //   // 21.422018, 3.363504, 0.035111,
  //   // 0.223619, 0.035111, 0.000367
  //   0.893508, 0.326783, 0.180538, 
  //   0.235999, 0.086781, 0.048497, 
  //   0.172893, 0.062554, 0.036107
  // };
  for(int k=0;k<9;k++) rgb[k] = 0.0f;
  int N = 15;
  int M = 15;
  // sample lambda_out
  for(int i=0;i<M;i++)
  {
    // and get rgb weighting factors
    // const float lambda_out = 400.0f + 20.0f*i;
    const float lambda_out = 400.0f + (700.0-400.0)*((.5f+i)/(float)M);
    // const float lambda_out = 400.0f + (700.0-400.0)*(drand48());
    // spectrum_p_to_rgb(lambda_out, 1.0f, rgbout);
    // rgbout[0] = get_spectrum(rgb_red, lambda_out);
    // rgbout[1] = get_spectrum(rgb_green, lambda_out);
    // rgbout[2] = get_spectrum(rgb_blue, lambda_out);
    get_rgb(lambda_out, rgbout);
    // printf("lambda_out %f, rgb %f %f %f\n", lambda_out, rgbout[0], rgbout[1], rgbout[2]);
    // for(int k=0;k<=i;k++)
    for(int k=0;k<N;k++)
    {
      // const float lambda_in = 400.0f + 20.0f*k;
      const float lambda_in = 400.0f + (700.0-400.0)*((.5f + k)/(float)N);
      // const float lambda_in = 400.0f + (700.0-400.0)*(drand48());
      // spectrum_p_to_rgb(lambda_in, 1.0f, rgbin);
      // rgbin[0] = get_spectrum(rgb_red, lambda_in);
      // rgbin[1] = get_spectrum(rgb_green, lambda_in);
      // rgbin[2] = get_spectrum(rgb_blue, lambda_in);
      get_rgb(lambda_in, rgbin);
      float rd;
      if(adjoint) rd = s->readout(theta_out, theta_in,  phi_rel, lambda_in, lambda_out)*PREBAKED_COS(cosf(theta_out));
      else        rd = s->readout(theta_in,  theta_out, phi_rel, lambda_in, lambda_out)*PREBAKED_COS(cosf(theta_in));
      const float f = 100.0f*100.0f*rd*ENERGY_SCALE/((float)M*N);
      // printf("lambda_in, p(%f) = %f\n", lambda_in, f);
      for(int n=0;n<3;n++) for(int m=0;m<3;m++) rgb[3*n+m] += rgbin[n]*rgbout[m]*f;
    }
  }
  // for(int k=0;k<9;k++) rgb[k] *= 1.f/rgbint[k];
  // for(int n=0;n<3;n++) for(int m=0;m<3;m++) if(n < m) rgb[3*n+m] = 0.0;
  // for(int n=0;n<3;n++) for(int m=0;m<3;m++) rgb[3*n+m] *= 1.0f/(rgbint[m]*rgbint[n]);
  // printf("M:  lambda_out"); for(int n=0;n<3;n++) {printf("\n"); for(int m=0;m<3;m++) printf("%f ", rgb[3*n+m]); }
}

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

extern float sample(const float *omega_in, rayhit_t *hit, float *omega_out, const float x1, const float x2, const float rr, void *data)
{
  float x, y, z;
  fluor_t *f = (fluor_t *)data;
  // multispecBRDF *s = (multispecBRDF *)data;
  const float su = sqrtf(x1);
  x = su*cosf(2.f*M_PI*x2);
  y = su*sinf(2.f*M_PI*x2);
  z = sqrtf(1.0 - x1);
  for(int k=0;k<3;k++) omega_out[k] = hit->a[k]*x + hit->b[k]*y + hit->normal[k]*z;
  const float theta_in  = acosf(fminf(.999f, fmaxf(0.0f, - dotproduct(hit->normal, omega_in))));
  const float theta_out = acosf(z);
  float phi_rel = fabsf(2.f*M_PI*x2 - atan2f(- dotproduct(hit->b, omega_in), - dotproduct(hit->a, omega_in)));
  if(phi_rel > M_PI) phi_rel = 2.0f*M_PI - phi_rel;

  float M[9], rgbin[3], rgbout[3];
  get_rgbxrgb(theta_in, theta_out, phi_rel, f->s, M, hit->adjoint);
  // spectrum_p_to_rgb(hit->lambda_in, 1.0f, rgbin);
  get_rgb(hit->lambda_in, rgbin);
  //     rgbin[0] = get_spectrum(rgb_red,   hit->lambda_in);
  //     rgbin[1] = get_spectrum(rgb_green, hit->lambda_in);
  //     rgbin[2] = get_spectrum(rgb_blue,  hit->lambda_in);
  // for(int k=0;k<3;k++) rgbin[k] = 1.0f;
  for(int n=0;n<3;n++) {rgbout[n] = .0f; for(int m=0;m<3;m++) rgbout[n] += M[3*n+m]*rgbin[m];}
  return M_PI*eval_rgb(hit->lambda, rgbout)/3.0f;
  // return fmaxf(0.0f, spectrum_rgb_to_p(hit->lambda, rgbout))/3.0f;
}

extern float brdf(float *omega_in, rayhit_t *hit, float *omega_out, void *data)
{
  fluor_t *f = (fluor_t *)data;
  // multispecBRDF *s = (multispecBRDF *)data;
  const float dot_in = - dotproduct(hit->normal, omega_in);
  const float dot_out = dotproduct(hit->normal, omega_out);
  if(dot_in < 0.001f || dot_out < 0.001f) return 0.0f;
  const float theta_in  = acosf(fminf(.999f, dot_in));
  const float theta_out = acosf(fminf(.999f, dot_out));
  float phi_rel = fabsf(atan2f(dotproduct(hit->b, omega_out), dotproduct(hit->a, omega_out)) - atan2f(- dotproduct(hit->b, omega_in), - dotproduct(hit->a, omega_in)));
  if(phi_rel > M_PI) phi_rel = 2.0f*M_PI - phi_rel;

  float M[9], rgbin[3], rgbout[3];
  get_rgbxrgb(theta_in, theta_out, phi_rel, f->s, M, hit->adjoint);
  // spectrum_p_to_rgb(hit->lambda_in, 1.0f, rgbin);
  get_rgb(hit->lambda_in, rgbin);
  //     rgbin[0] = get_spectrum(rgb_red,   hit->lambda_in);
  //     rgbin[1] = get_spectrum(rgb_green, hit->lambda_in);
  //     rgbin[2] = get_spectrum(rgb_blue,  hit->lambda_in);
  // for(int k=0;k<3;k++) rgbin[k] = 1.0f;
  for(int n=0;n<3;n++) {rgbout[n] = 0.0f; for(int m=0;m<3;m++) rgbout[n] += M[3*n+m]*rgbin[m]; }
  return eval_rgb(hit->lambda, rgbout)/(3.0f);
  // return fmaxf(0.0f, spectrum_rgb_to_p(hit->lambda, rgbout)/M_PI)/3.0f;
}

extern int init(FILE *f, void **data)
{
  fluor_t *s = (fluor_t *)malloc(sizeof(fluor_t));
  *data = s;
  if(fscanf(f, "%s", s->filename) != 1)
  {
    printf("[fluor_rgbxrgb::init] ERROR: could not parse arguments!\n");
    printf(" expected: <data filename>\n");
    return 1;
  }
  // multispecBRDF *s = new multispecBRDF(string(filename));

#if 0
  float M[9];
  get_rgbxrgb(0.5f, 0.3f, 2.0, s->s, M, 1);
  printf("M:  "); for(int n=0;n<3;n++) {printf("\n"); for(int m=0;m<3;m++) printf("%f ", M[3*n+m]); }
  printf("\n");
  //exit(0);
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

int uvn_init(const char *fname, unsigned int tri, shader_so_t *self)
{
  fluor_t *f = (fluor_t *)self->data;
  shader_data_t *data = shader_data(f->filename, SDAT_RAW);
  if(!data)
  {
    fprintf(stderr, "[fluor_rgbxrgb]: could not open file `%s' !\n", f->filename);
    return 1;
  }

  f->s = new multispecBRDF((char *)data->data);
  f->s->smoothSpecularHighlight();
  return 0;
}

extern void cleanup(void *data)
{
}

}

