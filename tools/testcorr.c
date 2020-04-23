// gcc -std=c11 -D_XOPEN_SOURCE -Wall testcorr.c -lm -o testcorr
//
// gnuplot> binwidth=0.01
// gnuplot> bin(x,width)=width*floor(x/width)
// gnuplot> plot 'testcorr.data' using (bin($1,binwidth)):(1.0) smooth freq with boxes title 'p(x)', '' using (bin($2,binwidth)):(1.0) smooth freq with boxes title 'p(x)*p(y|x)'
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// domain is x in [0, 1]
float eval(float x)
{
  // return x;   // integral: 0.5f
  return x*x; // integral: 1/3
  // return x*x*x - 0.1*x; 
  // return expf(-x*x);
}

float sample(float r, float *pdf)
{
  // integrate f by importance sampling:
  // f(x) = x
  // p(x) = f(x)*2 (normalised)
  // x = sqrt(r)
  const float x = sqrtf(r);
  *pdf = x*2.0f;
  return x;
}

int main(int argc, char *arg[])
{
  srand48(666);

  float E1 = 0.0f;
  float E2 = 0.0f;
  const int n = 10000;
  for(int i=0;i<n;i++)
  {
    float px;
    const float x = sample(drand48(), &px);
    const float fx = eval(x);
    E1 += fx/px;
    // and with correlated sampling
    const float r2 = drand48();
    const float size = 0.5f*x;
    float y = x + size * (r2-.5f);
    if(y < 0) y++;
    if(y > 1) y--;
    const float pyx = 1.0f/size;
    // FIXME: this is wrong.
    const float py = px * pyx;
    const float fy = eval(y);
    E2 += fy / py;
    // plot sample histograms for gnuplot script above:
    fprintf(stderr, "%f %f\n", x, y);
  }
  E1 /= n;
  E2 /= n;
  // fprintf(stderr, "integral f(x) dx approx %f %f\n", E1, E2);
  exit(0);
}
