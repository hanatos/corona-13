#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

#include "prims.h"
#include "geo.h"


int main (int argc, char *arg[])
{
  if(argc < 2)
  {
    fprintf(stderr, "usage: %s [sphere,cylinder] [v0: x y z] [v1: x y z] [r0 r1]\n", arg[0]);
    exit(1);
  }

  float v0[3]={-3.0,0.0,0.0}, v1[3]={-3.0, 2.0, 0.0}, r0=1.0, r1=1.0;

  uint32_t num_verts = 0, num_vtxidx = 0;

  prims_shape_t shape;
  primid_t *prims = 0;

  if(!strcmp(arg[1], "sphere"))
  {
    if(argc >= 6)
    {
      for(int k=0;k<3;k++) v0[k] = atof(arg[2+k]);
      r0 = atof(arg[5]);
    }
    num_vtxidx = 1;
    num_verts = 2;
    shape.num_prims = 1;
    shape.vtxidx = malloc(sizeof(prims_vtxidx_t) * num_vtxidx);
    shape.vtx = malloc(sizeof(prims_vtx_t)*num_verts);
    prims = malloc(sizeof(primid_t)*shape.num_prims);

    primid_t pid;
    pid.shapeid = 0;
    pid.vi = 0;
    pid.mb = 1;  // motion blur
    pid.vcnt = 1; // sphere
    prims[0] = pid;

    shape.vtxidx[0].v = 0;
    shape.vtxidx[0].uv = geo_encode_uv(0, 0);
    shape.vtx[0].n = touint(r0);  // spheres store radius instead of normal index
    shape.vtx[0].v[0] = v0[0];
    shape.vtx[0].v[1] = v0[1];
    shape.vtx[0].v[2] = v0[2];
    shape.vtx[1].v[0] = v0[0]+2.0f;
    shape.vtx[1].v[1] = v0[1];
    shape.vtx[1].v[2] = v0[2];
  }
  else if(!strcmp(arg[1], "cylinder") || !strcmp(arg[1], "cone"))
  {
    if(argc >= 10)
    {
      for(int k=0;k<3;k++) v0[k] = atof(arg[2+k]);
      for(int k=0;k<3;k++) v1[k] = atof(arg[5+k]);
      r0 = atof(arg[8]);
      r1 = atof(arg[9]);
      fprintf(stderr, "radii = %g %g\n", r0, r1);
    }
    // only difference between cylinder and cone:
    if(!strcmp(arg[1], "cone")) r1 = 0.f;

    num_vtxidx = 2;
    num_verts = 2;
    shape.num_prims = 1;
    shape.vtxidx = (prims_vtxidx_t *)malloc(sizeof(prims_vtxidx_t) * num_vtxidx);
    shape.vtx = (prims_vtx_t *)malloc(sizeof(prims_vtx_t)*num_verts);
    prims = (primid_t *)malloc(sizeof(primid_t)*shape.num_prims);

    primid_t pid;
    pid.shapeid = 0;
    pid.vi = 0;
    pid.mb = 0;  // no motion blur
    pid.vcnt = 2; // cylinder
    prims[0] = pid;

    shape.vtxidx[0].v = 0;
    shape.vtxidx[0].uv = geo_encode_uv(0, 0);
    shape.vtx[0].n = touint(r0);  // radius instead of normal index
    shape.vtx[0].v[0] = v0[0];
    shape.vtx[0].v[1] = v0[1];
    shape.vtx[0].v[2] = v0[2];

    shape.vtxidx[0].v = 1;
    shape.vtxidx[0].uv = geo_encode_uv(1, 1);
    shape.vtx[1].n = touint(r1);  // radius instead of normal index
    shape.vtx[1].v[0] = v1[0];
    shape.vtx[1].v[1] = v1[1];
    shape.vtx[1].v[2] = v1[2];
  }
  else exit(2);

  char geoname[1024];
  snprintf(geoname, 1024, "%s.geo", arg[1]);
  FILE *o = fopen(geoname, "wb");
  if(o)
  {
    prims_header_t header;
    header.magic = GEO_MAGIC;
    header.version = GEO_VERSION;
    header.num_prims = shape.num_prims;
    uint64_t cnt = sizeof(prims_header_t) + sizeof(primid_t)*shape.num_prims;
    header.vtxidx_offset = cnt;
    cnt += sizeof(prims_vtxidx_t)*num_vtxidx;
    header.vertex_offset = (cnt + 0xf) & ~0xf; // round up to next multiple of 16 bytes for sse alignment
    assert(header.vertex_offset >= cnt);
    assert((header.vertex_offset & 0xf) == 0);
    fwrite(&header, sizeof(prims_header_t), 1, o);
    fwrite(prims, sizeof(primid_t), shape.num_prims, o);
    fwrite(shape.vtxidx, sizeof(prims_vtxidx_t), num_vtxidx, o);
    if(header.vertex_offset > cnt) // pad to sse alignment
      fwrite(shape.vtx, header.vertex_offset - cnt, 1, o);
    fwrite(shape.vtx, sizeof(prims_vtx_t), num_verts, o);
    fclose(o);
  }

  exit(0);
}
