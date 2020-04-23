#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

#include "prims.h"
#include "geo.h"

static inline int to_zero_base_idx(const int idx, uint64_t num_indices)
{
  return (idx <  0) ? idx + num_indices : idx - 1;
}

static inline void write_shape(
  const prims_shape_t *shape,
  const primid_t *prims,
  const int shape_num_vtxidx,
  const int shape_num_verts,
  const char *shapename,
  const int stride,
  const int recompute_normals)
{
  fprintf(stderr, "writing shape %s with %ld prims and %d vi and %d v\n", shapename, shape->num_prims, shape_num_vtxidx, shape_num_verts);
  if(shape->num_prims)
  {
    char geoname[1024];
    snprintf(geoname, 1024, "%s.geo", shapename);
    FILE *o = fopen(geoname, "wb");
    if(o)
    {
      prims_header_t header;
      header.magic = GEO_MAGIC;
      header.version = GEO_VERSION;
      header.num_prims = shape->num_prims;
      uint64_t cnt = sizeof(prims_header_t) + sizeof(primid_t)*shape->num_prims;
      header.vtxidx_offset = cnt;
      cnt += sizeof(prims_vtxidx_t)*shape_num_vtxidx;
      header.vertex_offset = (cnt + 0xf) & ~0xf; // round up to next multiple of 16 bytes for sse alignment
      assert(header.vertex_offset >= cnt);
      assert((header.vertex_offset & 0xf) == 0);
      fwrite(&header, sizeof(prims_header_t), 1, o);
      fwrite(prims, sizeof(primid_t), shape->num_prims, o);
      fwrite(shape->vtxidx, sizeof(prims_vtxidx_t), shape_num_vtxidx, o);
      if(header.vertex_offset > cnt) // pad to sse alignment
        fwrite(shape->vtx, header.vertex_offset - cnt, 1, o);
      fwrite(shape->vtx, sizeof(prims_vtx_t)*stride, shape_num_verts, o);
      fclose(o);

      if(recompute_normals)
      {
        // create vertex normals on geo by loading it back in
        prims_t prims;
        prims_init(&prims);
        prims_allocate(&prims, 1);
        char basename[1024];
        strncpy(basename, geoname, 1024);
        basename[strlen(geoname)-4] = '\0';
        if(prims_load_with_flags(&prims, basename, "none", 0, 'w', 0))
        {
          fprintf(stderr, "[write_shape] failed to read back `%s' to compute normals!\n", geoname);
          return;
        }
        geo_recompute_normals(prims.shape);
        prims_cleanup(&prims);
      }
    }

#if 0 // integrity check:
    int v_broken = 0;
    for(int v=0;v<shape_num_vtxidx;v++)
    {
      if(!v_broken && shape_num_verts <= shape->vtxidx[v].v)
      {
        v_broken = 1;
        fprintf(stderr, "vertex index range broken: %d / %d\n", shape->vtxidx[v].v, shape_num_verts);
      }
    }
#endif
  }
}

static inline void load_obj_lists(
    FILE *f,
    float *v,
    float *n,
    float *uv,
    uint64_t num_verts,
    uint64_t num_normals,
    uint64_t num_uvs)
{
  char line[2048];
  int lineno = 1;
  uint64_t ni = 0, vi = 0, ti = 0;
  while(fscanf(f, "%[^\n]", line) != EOF)
  {
    if(!strncmp(line, "vn ", 3))
    {
      if(ni >= num_normals) fprintf(stderr, "obj has fewer normals than expected (%lu)\n", num_normals);
      assert(ni < num_normals);
      const int cnt = sscanf(line, "vn %f %f %f", n + 3*ni, n+3*ni+1, n+3*ni+2);
      if(cnt == 3) ni++;
      else fprintf(stderr, "line %d: weird normal: `%s'\n", lineno, line);
    }
    else if(uv && !strncmp(line, "vt ", 3))
    {
      if(ti >= num_uvs) fprintf(stderr, "obj has fewer texture coordinates than expected (%lu)\n", num_uvs);
      assert(ti < num_uvs);
      const int cnt = sscanf(line, "vt %f %f", uv + 2*ti, uv+2*ti+1);
      if(cnt == 2) ti++;
      else fprintf(stderr, "line %d: weird uvs: `%s'\n", lineno, line);
    }
    else if(!strncmp(line, "v ", 2))
    {
      if(vi >= num_verts) fprintf(stderr, "obj has fewer vertices than expected (%lu)\n", num_verts);
      assert(vi < num_verts);
      const int cnt = sscanf(line, "v %f %f %f", v + 3*vi, v+3*vi+1, v+3*vi+2);
      if(cnt == 3) vi++;
      else fprintf(stderr, "line %d: weird vertex: `%s'\n", lineno, line);
    }
    (void)fgetc(f); // munch newline
    lineno++;
    // invalidate line in case the next line will be a blank:
    line[0] = '\0';
  }
}

int main (int argc, char *arg[])
{
  if(argc < 2)
  {
    fprintf(stderr, "usage: %s input.obj [input_shutter_close.obj] [particle radius (def 0.001)]\n", arg[0]);
    exit(1);
  }
  FILE *f = fopen(arg[1], "rb");
  if(!f)
  {
    fprintf(stderr, "could not open `%s'\n", arg[1]);
    exit(1);
  }

  // have two obj and thus motion blur?
  int motion = 0;
  FILE *f2 = 0;
  if(argc > 2)
  {
    f2 = fopen(arg[2], "rb");
    if(f2)
    {
      fprintf(stderr, "given two obj files, interpreting as motion blurred geo!\n");
      motion = 1;
    }
  }
  const int stride = motion?2:1;


  // count vertices only in first file. on mismatch we just die.
  uint64_t num_verts = 0, num_normals = 0, num_uvs = 0, num_faces = 0;
  char line[2048];
  // first pass: count vertices, normals, and uvs:
  while(fscanf(f, "%[^\n]", line) != EOF)
  {
    if(!strncmp(line, "v ", 2)) num_verts++;
    (void)fgetc(f); // munch newline
    line[0] = '\0'; // invalidate to support blank lines
  }
  num_normals = num_uvs = num_faces = num_verts;
  // rewind
  fseek(f, 0, SEEK_SET);

  fprintf(stderr, "num verts %lu num normals %lu num texture coords %lu\n", num_verts, num_normals, num_uvs);
  // alloc
  // uint32_t *vmap  = (uint32_t *)malloc(sizeof(uint32_t)*num_verts);
  // uint32_t *nmap  = (uint32_t *)malloc(sizeof(uint32_t)*num_normals);
  // uint32_t *uvmap = (uint32_t *)malloc(sizeof(uint32_t)*num_uvs);
  float *v2 = 0, *n2 = 0;
  if(motion)
  { // allocate shutter close positions and normals:
    v2 = (float *)malloc(sizeof(float)*3*num_verts);
    n2 = (float *)malloc(sizeof(float)*3*num_normals);
  }
  float *v  = (float *)malloc(sizeof(float)*3*num_verts);
  float *n  = (float *)malloc(sizeof(float)*3*num_normals);
  float *uv = (float *)malloc(sizeof(float)*2*num_uvs);

  prims_shape_t shape;
  shape.num_prims = 0;
  shape.vtxidx = (prims_vtxidx_t *)malloc(sizeof(prims_vtxidx_t) * num_faces);
  shape.vtx = (prims_vtx_t *)malloc(sizeof(prims_vtx_t)*stride*num_verts);
  uint32_t shape_num_verts = 0;
  uint32_t shape_num_vtxidx = 0;
  primid_t *prims = (primid_t *)malloc(sizeof(primid_t)*num_faces);

  const float radius = argc > (2+motion) ? atof(arg[argc-1]) : 0.001f;
  fprintf(stderr, "radius = %g\n", radius);
  const int radius_in_int = *(int *)&radius;
  char shapename[1024];
  snprintf(shapename, sizeof(shapename), "noname");

  // uint32_t last_hair_vert = -1;
  // uint32_t first_hair_vert = -1;
  // uint32_t hair_index = 0;

  int recompute_normals = 0;

  // again, load vertices this time
  load_obj_lists(f, v, n, uv, num_verts, num_normals, num_uvs);

  // in case of motion, open 2nd obj. assume vertex order lines up with first one.
  if(motion)
  {
    load_obj_lists(f2, v2, n2, 0, num_verts, num_normals, num_uvs);
    fclose(f2);
  }

  // now init shape with particles for every vertex:
  for(int vi=0;vi<num_verts;vi++)
  {
    primid_t pid;
    pid.shapeid = 0;
    pid.vi = vi;
    pid.mb = motion;
    pid.vcnt = 1; // sphere
    prims[vi] = pid;

    shape.vtxidx[vi].v = vi;
    shape.vtxidx[vi].uv = geo_encode_uv(0, 0);
    int vi2 = motion ? 2*vi : vi;
    shape.vtx[vi2].n = radius_in_int;
    shape.vtx[vi2].v[0] = v[3*vi+0];
    shape.vtx[vi2].v[1] = v[3*vi+1];
    shape.vtx[vi2].v[2] = v[3*vi+2];
    if(motion)
    {
      shape.vtx[vi2+1].n = radius_in_int;
      shape.vtx[vi2+1].v[0] = v2[3*vi+0];
      shape.vtx[vi2+1].v[1] = v2[3*vi+1];
      shape.vtx[vi2+1].v[2] = v2[3*vi+2];
    }
  }
  shape.num_prims = num_verts;
  shape_num_vtxidx = num_verts;
  shape_num_verts = motion? 2*num_verts:num_verts;

#if 0
  // third pass: load faces and write everything to shape struct 
  int lineno = 1;
  fseek(f, 0, SEEK_SET);
  int mtl = 0; // count materials
  while(fscanf(f, "%[^\n]", line) != EOF)
  {
    // if(!strncmp(line, "g ", 2) || !strncmp(line, "usemtl ", 7))
    if(!strncmp(line, "o ", 2) || !strncmp(line, "usemtl ", 7))
    // if(!strncmp(line, "usemtl ", 7))
    {
      // new file
      write_shape(
          &shape,
          prims,
          shape_num_vtxidx,
          shape_num_verts,
          shapename,
          stride,
          recompute_normals);
      recompute_normals = 0;
      char *beg = line;
      while(*beg != ' ') beg++; // skip to first whitespace
      beg++;
      // if(!strcmp(beg, "None")) // blender does this if you didn't assign a material
        // sprintf(beg, "mat_%04d", mtl);
      strcpy(shapename, beg);
      beg = shapename;
      // while(*beg != '.' && (beg < shapename + strlen(shapename))) beg++; // skip to first point
      // if(*beg == '.') *beg = 0;
      fprintf(stderr, "object `%s'\n", shapename);

      shape.num_prims = 0;
      shape_num_vtxidx = 0;
      shape_num_verts = 0;
      memset(vmap, 0xff, sizeof(uint32_t)*num_verts);
      memset(nmap, 0xff, sizeof(uint32_t)*num_normals);
      memset(uvmap, 0xff, sizeof(uint32_t)*num_uvs);
      mtl ++;
    }
    else if(!strncmp(line, "f ", 2) || !strncmp(line, "l ", 2))
    {
      // face. parse vert normal texture indices
      int vert[4] = {0,0,0,0}, norm[4] = {0,0,0,0}, uvco[4] = {0,0,0,0};

      // find number of '/' to guess format:
      char *c = line;
      int slashes = 0, doubleslashes = 0;
      for(;*c!=0;c++) if(*c == '/') { slashes++; if(c[1] == '/') doubleslashes = 1;}
      int cnt = 0;
      if(slashes != 8 && slashes != 6 && slashes != 4 && slashes != 3 && slashes != 0)
        fprintf(stderr, "line %d: wrong number of '/' in face definition! `%s'\n", lineno, line);
      // v/t/n and v//n
      if(slashes >= 6)
      {
        if(doubleslashes)
        { // v//n
          cnt = sscanf(line+1, " %d//%d %d//%d %d//%d %d//%d",
              vert+0, norm+0,
              vert+1, norm+1,
              vert+2, norm+2,
              vert+3, norm+3);
        }
        else
        { // v/t/n
          cnt = sscanf(line+1, " %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d",
              vert+0, uvco+0, norm+0,
              vert+1, uvco+1, norm+1,
              vert+2, uvco+2, norm+2,
              vert+3, uvco+3, norm+3);
        }
      }
      // v/t
      else if(slashes >= 3)
        cnt = sscanf(line+1, " %d/%d %d/%d %d/%d %d/%d",
            vert+0, uvco+0,
            vert+1, uvco+1,
            vert+2, uvco+2,
            vert+3, uvco+3);
      // v
      else
        cnt = sscanf(line+1, " %d %d %d %d",
            vert+0, vert+1, vert+2, vert+3);

      // obj indices are 1 based
      for(int k=0;k<4;k++)
      {
        vert[k] = to_zero_base_idx(vert[k], num_verts);
        uvco[k] = to_zero_base_idx(uvco[k], num_uvs);
        norm[k] = to_zero_base_idx(norm[k], num_normals);
      }

      primid_t pid;
      pid.shapeid = 0;
      pid.vi = shape_num_vtxidx;
      pid.mb = motion;
      if(cnt == 12 || cnt == 8 || cnt == 4) // quad
        pid.vcnt = 4;
      else if(cnt == 9 || cnt == 6 || cnt == 3) // tri
        pid.vcnt = 3;
      else if(cnt == 2) // uv and normal-less lines
        pid.vcnt = 2;
      else
        fprintf(stderr, "only lines, tris and quads supported so far\n");

      prims[shape.num_prims++] = pid;
      for(int k=0;k<pid.vcnt;k++)
      {
        prims_vtxidx_t vid;
        uint32_t n0 = 0, n1 = 0;
        if(pid.vcnt <= 2)
        { // those have radii, not normals
          n0 = radius_in_int;
          if(motion) n1 = radius_in_int;
        }
        else if(norm[k] >= 0)
        {
          n0 = geo_encode_normal(n + 3*norm[k]);
          if(motion) n1 = geo_encode_normal(n2 + 3*norm[k]);
        }
        else
        {
          // set normal to broken encoding, set flag to recompute on output.
          recompute_normals = 1;
          n0 = 0x8000;
          if(motion) n1 = 0x8000;
        }

        // duplicate vertex in case of inconsistent normal
        int reuse = 0;
        if(vmap[vert[k]] != (uint32_t)-1)
        {
          if(pid.vcnt <= 2 || (n0 == shape.vtx[stride*vmap[vert[k]]].n && (!motion || n1 == shape.vtx[stride*vmap[vert[k]]+1].n)))
          {
            reuse = 1;
            vid.v = vmap[vert[k]];
          }
        }

        if(!reuse)
        {
          vid.v = shape_num_verts++;
          assert(vid.v < 4*num_verts); // TODO: better bound
          for(int i=0;i<3;i++)
            shape.vtx[stride*vid.v].v[i] = v[3*vert[k] + i];
          // write second vertex in case of motion blur
          if(motion)
            for(int i=0;i<3;i++)
              shape.vtx[stride*vid.v+1].v[i] = v2[3*vert[k] + i];
          shape.vtx[stride*vid.v].n = n0;
          if(motion) shape.vtx[stride*vid.v+1].n = n1;
          vmap[vert[k]] = vid.v;
        }

        if(uvco[k] >= 0)
        {
          vid.uv = geo_encode_uv(uv[2*uvco[k]], uv[2*uvco[k]+1]);
        }
        else vid.uv = 0;
        shape.vtxidx[shape_num_vtxidx++] = vid;
      }

      // now fix some uvw into hair strands:
      if(pid.vcnt == 2)
      {
        if(first_hair_vert == -1 || vert[0] != last_hair_vert)
        { // first strand or new strand
          if(first_hair_vert != -1)
          { // finished strand, init uvw:
            for(uint32_t k=first_hair_vert;k<shape_num_vtxidx;k++)
              shape.vtxidx[k].uv = geo_encode_uvw(hair_index/10000.0f, 0, (k-first_hair_vert)/((float)shape_num_vtxidx-first_hair_vert-1.0f));
            hair_index ++;
          }
          first_hair_vert = shape_num_vtxidx-2;
          last_hair_vert = vert[1];
        }
        else if(vert[0] == last_hair_vert)
        { // continue same strand
          last_hair_vert = vert[1];
        }
      }
      else if(first_hair_vert != -1)
      { // strand finished, other faces now. init uvw:
        for(uint32_t k=first_hair_vert;k<shape_num_vtxidx-pid.vcnt;k++)
          shape.vtxidx[k].uv = geo_encode_uvw(hair_index/10000.0f, 0, (k-first_hair_vert)/((float)shape_num_vtxidx-pid.vcnt-first_hair_vert-1.0f));
        first_hair_vert = last_hair_vert = -1;
        hair_index ++;
      }
    }
    (void)fgetc(f); // munch newline
    lineno++;
    line[0] = '\0'; // support blanks
  }
#endif

  // close shape
  write_shape(
      &shape,
      prims,
      shape_num_vtxidx,
      shape_num_verts,
      shapename,
      stride,
      recompute_normals);

  fclose(f);
  exit(0);
}
