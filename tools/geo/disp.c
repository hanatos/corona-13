#include "prims.h"
#include "geo.h"
#include "texture.h"
#include "vdata.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>


int main(int argc, char *argv[])
{
  if(argc < 4)
  {
    fprintf(stderr, "[dispgeo] usage: dispgeo input.geo displacement-texture.pfm scale [--midlevel m] [--multiscale]\n");
    fprintf(stderr, "[dispgeo]        will overwrite input.geo with displaced vertices.\n");
    fprintf(stderr, "[dispgeo]        --midlevel m       : subtract this from displacement texture as an offset.\n");
    fprintf(stderr, "[dispgeo]        --multiscale       : create accompanying vdata with leadr data.\n");
    exit(1);
  }
  const float scale = atof(argv[3]);
  float midlevel = 0.0f;
  int multiscale = 0;
  for(int k=4;k<argc;k++)
  {
    if(!strcmp(argv[k], "--midlevel") && argc > k+1) midlevel = atof(argv[++k]);
    else if(!strcmp(argv[k], "--multiscale")) multiscale = 1;
  }

  prims_t prims;
  prims_init(&prims);
  prims_allocate(&prims, 1);
  if(prims_load_with_flags(&prims, argv[1], "none", 0, 'w', 0))
  {
    fprintf(stderr, "[geo-disp] could not load geo file `%s'\n", argv[1]);
    exit(2);
  }

  // open displacement texture
  texture_t tex;
  if(texture_load(&tex, argv[2], 1))
  {
    fprintf(stderr, "[geo-disp] could not open displacement texture `%s'!\n", argv[2]);
    exit(2); // will close f if still open
  }
  // open elliptical uv footprints
  vdata_t vdata = {0};
  char ewafile[1024];
  snprintf(ewafile, 1024, "%s.vdata", argv[1]);
  if(vdata_map(&vdata, ewafile))
  {
    fprintf(stderr, "[geo-disp] could not open uv footprints `%s'!\n", ewafile);
    exit(2);
  }

  const int mb = prims.shape[0].primid[0].mb+1;
  uint64_t num_verts = (prims.shape[0].data_size - ((uint8_t *)prims.shape[0].vtx - (uint8_t *)prims.shape[0].data))/
    (sizeof(prims_vtx_t) * mb);

  // this is stupid:
  vdata.num_verts = num_verts;
  vdata.num_slots = 5;

  if(multiscale)
  { // compute 5 leadr mapping moments for multiscale filtering:
    assert(tex.channels == 1);

    // compute leadr map in texture space, i.e. Sec 4.2:
    // E(u_n) E(v_n) E(u2_n) E(v2_n) and E(u_n v_n)
    // we'll compute the tangent frame statistics later per vertex using
    // Eq. (12) in that same section.
    texture_t leadr;
    texture_alloc(&leadr, tex.width, tex.height, 5);
    memset(leadr.tex, 0, sizeof(float)*leadr.channels*leadr.width*leadr.height);
    for(int j=1;j<leadr.height-1;j++)
    {
      for(int i=1;i<leadr.width-1;i++)
      { // locally compute dh/du and dh/dv for this texel:
        const float dh_du = .5f*scale*(tex.tex[tex.width*j + i+1] - tex.tex[tex.width*j + i-1]);
        const float dh_dv = .5f*scale*(tex.tex[tex.width*(j+1) + i] - tex.tex[tex.width*(j-1) + i]);
        // store in leadr map, finest level
        leadr.tex[leadr.channels*(leadr.width*j+i)+0] = dh_du;
        leadr.tex[leadr.channels*(leadr.width*j+i)+1] = dh_dv;
        leadr.tex[leadr.channels*(leadr.width*j+i)+2] = dh_du*dh_du;
        leadr.tex[leadr.channels*(leadr.width*j+i)+3] = dh_dv*dh_dv;
        leadr.tex[leadr.channels*(leadr.width*j+i)+4] = dh_du*dh_dv;
      }
    }
    char output[512];
    snprintf(output, sizeof(output), "%s-leadr.tex", argv[1]);
    texture_write(&leadr, output);
    texture_cleanup(&leadr);

    vdata_t vnormals;
    vdata_alloc(&vnormals, mb*num_verts, 3);

    for(uint64_t v=0;v<num_verts;v++)
    { // dump per-vertex base mesh normal
      float *vd = vdata_get(&vnormals, mb*v);
      for(int m=0;m<mb;m++)
        geo_decode_normal(prims.shape[0].vtx[mb*v+m].n, vd + 3*m);
    }
    snprintf(output, sizeof(output), "%s-leadr.vdata", argv[1]);
    vdata_write(&vnormals, output);
    vdata_cleanup(&vnormals);
  }

  // go through all vertices and displace them.
#pragma omp parallel for schedule(static) default(shared)
  for(uint64_t v=0;v<num_verts;v++)
  { // displace in normal direction by texture at uv coords
    float fu, fv, vn[3], tmp[3];
    const float *vd = vdata_get(&vdata, v);
    const float *d0 = vd + 1, *d1 = vd+3;
    geo_decode_uv(((uint32_t*)vd)[0], &fu, &fv);
    texture_lookup_ewa(&tex, fu, fv, d0, d1, tmp);
    // tmp[0] = sinf(M_PI*.5f*fu)*cosf(M_PI*6.0f*fv);
    const float texv = scale * (tmp[0]/12.0f - midlevel);
    for(int m=0;m<mb;m++)
    { // ignore potentially different texture filter for shutter open and close (shade once)
      geo_decode_normal(prims.shape[0].vtx[mb*v+m].n, vn);
      for(int i=0;i<3;i++) prims.shape[0].vtx[mb*v+m].v[i] += vn[i] * texv;
    }
  }

  geo_recompute_normals(prims.shape);

  // write back on close
  prims_cleanup(&prims);
  texture_cleanup(&tex);
  vdata_cleanup(&vdata);

  exit(0);
}
