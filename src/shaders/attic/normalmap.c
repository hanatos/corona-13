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
#include "prims.h"
#include "shader.h"

typedef struct
{
  int index;
  int num;
  char filename[256];
  const rt_t *rt;
  int host;
  int width;
  int height;
  float *uv;
  unsigned char *texture;
}
TextureShader;

int init(FILE *s, void **data)
{
  TextureShader *t = (TextureShader *)malloc(sizeof(TextureShader));
  *data = t;
  t->texture = NULL;
  t->uv = NULL;
  // read sample shader number
  t->filename[0] = '\0';
  if(fscanf(s, "%d %s", &(t->host), t->filename) < 2)
  {
    fprintf(stderr, "normalmap shader could not read all parameters! expecting: <num> [filename]\n");
    return 1;
  }
  int dreggn = 0;
  dreggn = fscanf(s, "%*[^\n]\n");
  if(dreggn != dreggn) return 1;
  //FILE *test = fopen(t->filename, "rb");
  //if(!test) t->filename[0] = '\0';
  //else fclose(test);
  if(t->filename[0] == '#') t->filename[0] = '\0';
  t->rt = NULL;
  return 0;
}

int shape_init(uint32_t shapeid, shader_so_t *self)
{
  TextureShader *s = (TextureShader*) self->data;
  if(s->rt)
    fprintf(stderr, "[normalmap] WARNING: this shader is applied to two .ra2 portions, which is likely to Segfault: %s\n", fname);
  // read normals
  char filename[255];
  strncpy(filename, fname, 255);
  char *c;
  for(c=filename + strlen(filename);*c!='.'&&c!=filename;c--);
  *(c+1) = 'u';
  *(c+2) = 'v';
  *(c+3) = '\0';
  shader_data_t *data = shader_data(filename, SDAT_RAW);
  if(!data)
  {
    fprintf(stderr, "could not open uv coords from `%s' !\n", filename);
    self->prepare = NULL;
    return 1;
  }

  s->rt = rt;
  s->uv  = (float *)data->data;
  s->num = data->filesize/(sizeof(float)*3*2);

  //TODO: print useful error message if filesize doesn't fit .ra2.

  data = shader_data(s->filename, SDAT_IMG);
  if(!data)
  {
    fprintf(stderr, "normalmap: can't open texture file `%s'!\n", s->filename);
    self->prepare = NULL;
    return 1;
  }
  s->texture = (unsigned char *)data->data;
  s->width = data->width;
  s->height = data->height;

  // map tri to 0
  s->index = tri;
  return 0;
}

float prepare(path_t *p, int vt, void *data)
{
  //  FIXME: generates NaN normals :(
  TextureShader *s = (TextureShader *)data;
  int index = 6*(p->v[vt].hit.prim - s->index);
  const float hu = p->v[vt].hit.u;
  const float hv = p->v[vt].hit.v;

  const float t1 = (1-hu-hv)*s->uv[index] + hv*s->uv[index+2] + hu*s->uv[index+4];
  const float t2 = (1-hu-hv)*s->uv[index+1] + hv*s->uv[index+2+1] + hu*s->uv[index+4+1];
  const float fu = (t1 - floorf(t1))*(s->width-1);
  const float fv = (t2 - floorf(t2))*(s->height-1);
  const int u = (int)fu;
  const int v = (int)fv;
  const int up = u+1 == s->width  ? u : u + 1;
  const int vp = v+1 == s->height ? v : v + 1;
  const float wu = fu - u;
  const float wv = fv - v;

  // trunc:
  //const unsigned char *rgb = s->texture + 3*s->width*v + 3*u;
  // bilinear filter
  const unsigned char *rgb00 = s->texture + 3*s->width*v + 3*u;
  const unsigned char *rgb01 = s->texture + 3*s->width*v + 3*up;
  const unsigned char *rgb10 = s->texture + 3*s->width*vp + 3*u;
  const unsigned char *rgb11 = s->texture + 3*s->width*vp + 3*up;
  float rgb[3];
  for(int k=0;k<3;k++) rgb[k] = rgb00[k]*(1-wu)*(1-wv) + rgb01[k]*wu*(1-wv) + rgb10[k]*(1-wu)*wv + rgb11[k]*wu*wv;
  float n[3];
  n[0] = 1.0f-(rgb[0]*(2.0f/255.0f));
  n[1] = (rgb[1]*(2.0f/255.0f))-1.0f;
  n[2] = (rgb[2]*(2.0f/255.0f))-1.0f;

#if 1
  // fast global axis aligned tangent system
  // float a[3], b[3];
  // getPerpendicular(hit->normal, a);
  // crossproduct(hit->normal, a, b);
#else
  float invm[9];

  // find tangent vectors a and b in world space which correspond to (0,1) and (1,0) in texture space:
  const float m[9] = {s->uv[index  ], s->uv[index+2], s->uv[index+4],
                      s->uv[index+1], s->uv[index+3], s->uv[index+5],
                      1.0f,           1.0f,           1.0f};
#define A(y, x) m[(y - 1) * 3 + (x - 1)]
#define B(y, x) invm[(y - 1) * 3 + (x - 1)]
    const float det =
      A(1, 1) * (A(3, 3) * A(2, 2) - A(3, 2) * A(2, 3)) -
      A(2, 1) * (A(3, 3) * A(1, 2) - A(3, 2) * A(1, 3)) +
      A(3, 1) * (A(2, 3) * A(1, 2) - A(2, 2) * A(1, 3));
    if(det == 0.0f) return 1.0f;

    const float invDet = 1.f / det;
    B(1, 1) =  invDet * (A(3, 3) * A(2, 2) - A(3, 2) * A(2, 3));
    B(1, 2) = -invDet * (A(3, 3) * A(1, 2) - A(3, 2) * A(1, 3));
    B(1, 3) =  invDet * (A(2, 3) * A(1, 2) - A(2, 2) * A(1, 3));

    B(2, 1) = -invDet * (A(3, 3) * A(2, 1) - A(3, 1) * A(2, 3));
    B(2, 2) =  invDet * (A(3, 3) * A(1, 1) - A(3, 1) * A(1, 3));
    B(2, 3) = -invDet * (A(2, 3) * A(1, 1) - A(2, 1) * A(1, 3));

    B(3, 1) =  invDet * (A(3, 2) * A(2, 1) - A(3, 1) * A(2, 2));
    B(3, 2) = -invDet * (A(3, 2) * A(1, 1) - A(3, 1) * A(1, 2));
    B(3, 3) =  invDet * (A(2, 2) * A(1, 1) - A(2, 1) * A(1, 2));
#undef A
#undef B
  const float hitUv[3] = {u+4.0f, v, 1.0f};
  float ba[3], a[3], b[3];
  for(int k=0;k<3;k++)
  {
    ba[k] = 0.0f;
    for(int l=0;l<3;l++) ba[k] += invm[3*k + l]*hitUv[l];
  }
  const prim_t *t = rt->prims->prims + hit->prim;
  // FIXME: only works for triangles (and beziers..?)
  for(int k=0;k<3;k++) a[k] = ba[0]*t->v[0][k] + ba[1]*t->v[1][k] + ba[2]*t->v[2][k] - hit->hit[k];
  crossproduct(a, hit->normal, b);
  normalise(b);
  crossproduct(b, hit->normal, a);

  float matrix[9];
  for(int k=0;k<3;k++)
  {
    matrix[k]   = a[k];
    matrix[3+k] = b[k];
    matrix[6+k] = hit->normal[k];
  }
  hit->normal[0] = matrix[0]*n[0] + matrix[3]*n[1] + matrix[6]*n[2];
  hit->normal[1] = matrix[1]*n[0] + matrix[4]*n[1] + matrix[7]*n[2];
  hit->normal[2] = matrix[2]*n[0] + matrix[5]*n[1] + matrix[8]*n[2];
#endif
  // assert(finite(hit->normal[0]));
  // assert(finite(hit->normal[1]));
  // assert(finite(hit->normal[2]));
  return 1.0f;
}

