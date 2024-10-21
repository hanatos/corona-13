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
// #include "texture.h"
// #include "framebuffer.h"

#include <math.h>
#include <assert.h>

#if 1  ////===========================================================================
// random shader toy steal: https://www.shadertoy.com/view/4lB3zz
float noise(int x,int y)
{
    float fx = (float)x;
    float fy = (float)y;


    float f = sin((fx*12.9898+ fy*78.233)) * 43758.5453;
    return 2.0 * (f-(int)f) - 1.0;
}

float smoothNoise(int x,int y)
{
  return noise(x,y);
    // return noise(x,y)/4.0+(noise(x+1,y)+noise(x-1,y)+noise(x,y+1)+noise(x,y-1))/8.0+(noise(x+1,y+1)+noise(x+1,y-1)+noise(x-1,y+1)+noise(x-1,y-1))/16.0;
}

float COSInterpolation(float x,float y,float n)
{
    float r = n*3.1415926;
    float f = (1.0-cos(r))*0.5;
    return x*(1.0-f)+y*f;

}

float InterpolationNoise(float x, float y)
{
    int ix = (int)x;
    int iy = (int)y;
    float fracx = x-(int)x;
    float fracy = y-(int)y;

    float v1 = smoothNoise(ix,iy);
    float v2 = smoothNoise(ix+1,iy);
    float v3 = smoothNoise(ix,iy+1);
    float v4 = smoothNoise(ix+1,iy+1);

   	float i1 = COSInterpolation(v1,v2,fracx);
    float i2 = COSInterpolation(v3,v4,fracx);

    return COSInterpolation(i1,i2,fracy);

}

float PerlinNoise2D(float x,float y)
{
    float sum = 0.0;
    float frequency = 0.0;
    float amplitude = 0.0;
    // for(int i=firstOctave;i<octaves + firstOctave;i++)
    for(int i=0;i<2;i++)
    {
        frequency = pow(2.0,i);
        amplitude = pow(1.678,i);
        sum = sum + InterpolationNoise(x*frequency,y*frequency)*amplitude;
    }

    return sum;
}
#endif  ////===========================================================================


typedef struct tex_t
{
#if 0
  int num;
  tex_slot_t slot;
  framebuffer_t fb;
  float mul;
#endif
}
tex_t;

int init(FILE *s, void **data)
{
  int dreggn = fscanf(s, "%*[^\n]\n");
  return dreggn == -1;
#if 0
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
#endif
}

void cleanup(void *data)
{
#if 0
  tex_t *t = data;
  fb_cleanup(&t->fb);
#endif
}

float prepare(path_t *p, int v, void *data)
{
#if 0 // hack: conspire with brdf evaluation later, and store intermediates with vertex normals:
  memcpy(p->v[v].diffgeo.dpdu, p->v[v].hit.n, sizeof(float)*3);
#endif

#if 0
  // use geo normal!
  for(int k=0;k<3;k++)
    p->v[v].hit.n[k] = p->v[v].hit.gn[k];
#endif

#if 1
  // overwrite p->v[v].hit.n by bump map!
  // use simple perlin noise
  // use dpdu/dpdv as frame
  // use hit.{s,t} as texture coordinates
  float s = 13.0*(p->v[v].hit.x[2]-p->v[v].hit.x[0]); // .s;
  float t = 13.0*(p->v[v].hit.x[1]-p->v[v].hit.x[2]); // .t;
  float du = PerlinNoise2D(s, t);
  float dv = PerlinNoise2D(t, s);

  // has texture seams on sphere:
  // float *dpdu = p->v[v].diffgeo.dpdu;
  // float *dpdv = p->v[v].diffgeo.dpdv;

  // use some camera dependent bs for the test only:
  float dpdu[3], dpdv[3];
  crossproduct(p->e[v].omega, p->v[v].hit.n, dpdu);
  crossproduct(p->v[v].hit.n, dpdu, dpdv);

  // fprintf(stderr, "du dv %g %g\n", du, dv);
  for(int k=0;k<3;k++)
    p->v[v].hit.n[k] += 0.06 * (dpdu[k] * du + dpdv[k] * dv);
  normalise(p->v[v].hit.n);
#endif

#if 1 // conty's bump terminator fix doesn't work here because omega -> light not yet known.
  // Return alpha ^2 parameter from normal divergence
// float bump_alpha2 ( float3 N , float3 Nbump ) {
float cos_d = fminf(fabsf(dotproduct(p->v[v].hit.gn , p->v[v].hit.n)), 1.0f);
float tan2_d = (1 - cos_d * cos_d ) / ( cos_d * cos_d );
float alpha2 = fminf(fmaxf(0.125f * tan2_d, 0.0f), 1.0f);
// }
// Shadowing factor
// float b u m p _ s h a d o w i n g _ f u n c t i o n ( float3 N , float3 Ld , float alpha2 ) {
float cos_i = fmaxf(fabsf(dotproduct(p->v[v].hit.gn, p->e[v].omega)) , 1e-6f);
float tan2_i = (1 - cos_i * cos_i ) / ( cos_i * cos_i ) ;
p->v[v].shading.rd *= 2.0f / (1 + sqrtf (1 + alpha2 * tan2_i ) ) ;
// }
#endif

  return 1.0f;
#if 0
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
#endif
}

