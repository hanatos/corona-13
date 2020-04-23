// use dbor buffers to denoise.
// works with anisotropic diffusion where blur radii and diffusion tensors
// are derived from the dbor input: radius is determined to cover a certain
// minimum amount of samples, diffusion tensor comes from edges in the
// image with more confidence.
#include "corona_common.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <float.h>

#include <sys/stat.h>
#include <sys/mman.h>

#define MAX_DBORS 20

// parameters:
// number of nb to determine support of blur kernel
#define DBOR_NUM_NB 3
// pull push overwrite threshold (keep glints)
// increase this to overwrite more glints
#define DBOR_GLINT_THRS 0.05f
// max blur radius clamp
#define DBOR_MAX_RAD 200
// bias for edge detection
// make the bias larger : more blur/less sharpness
// make the bias smaller: smaller blotches, more edges
#define DBOR_EDGE_BIAS 0.0005f;

typedef struct dbor_buffer_s
{
  int fd;              // file descriptor
  void *data;          // mapped buffer
  size_t data_size;    // size of file
  
  int width;
  int height;
  float* pixel;
}
dbor_buffer_t;


int main(int argc, char *arg[])
{
  if(argc < 2)
  {
    fprintf(stderr, "[dbord] usage: input_file [radius].\n");
    fprintf(stderr, "        samples: blur radius determined for this #samples, default DBOR_NUM_NB\n");
    exit(1);
  }

  const float sample_count = argc > 2 ? atof(arg[2]) : DBOR_NUM_NB;
  // we'll need these for diffusion and don't want to
  // cross their declaration by jumping to the fail: label.
  float *tensor = 0, *buf0 = 0, *buf1 = 0;
  int mip_max_level = -1;
  float **mip_buf = 0;
    
  // ============================================
  //   load dbor cascade
  // ============================================
  int num_dbors = 0;
  dbor_buffer_t dbor_cascade[MAX_DBORS] = {{ 0 }};
  while (1)
  {
    if(num_dbors >= MAX_DBORS)
    {
      fprintf(stderr, "[dbor] MAX_DBORS of %d reached.\n", MAX_DBORS);
      break;
    }

    char filename[1024];
    snprintf(filename, sizeof(char)*1024, "%s_dbor%02d.pfm", arg[1], num_dbors);
    if(access(filename, F_OK))
    {
      break;
    }
    else
    {
      dbor_buffer_t* db = &dbor_cascade[num_dbors];

      db->fd = open(filename, O_RDONLY);
      if (db->fd == -1) {
        fprintf(stderr, "[dbor] could not open dbor cascade buffer `%s'\n", filename);
        continue;
      }

      // mmap data
      db->data_size = lseek(db->fd, 0, SEEK_END);
      lseek(db->fd, 0, SEEK_SET);
      db->data = mmap(0, db->data_size, PROT_READ, MAP_SHARED | MAP_NORESERVE, db->fd, 0);
      
      // get pfm header for faster grabbing later on.
      // no error handling is done, no comments supported.
      char *endptr;
      db->width = strtol(db->data+3, &endptr, 10);
      db->height = strtol(endptr, &endptr, 10);
      endptr++; // remove newline
      while(endptr < (char *)db->data + 100 && *endptr != '\n') endptr++;
      db->pixel = (float*)(++endptr); // remove second newline

      num_dbors++;
    }
  }

  // create file for merged dbor
  dbor_buffer_t dm = { 0 };
  char filename[1024];
  snprintf(filename, sizeof(filename), "%s_dbord.pfm", arg[1]);
  dm.fd = open(filename, O_RDWR | O_CREAT | O_TRUNC, (mode_t)0644);
  if (dm.fd == -1) {
    fprintf(stderr, "[dbor] could not open dbor merged buffer `%s'\n", filename);
    goto fail;
  }
  dm.data_size = dbor_cascade[0].data_size;
  if (lseek(dm.fd, dm.data_size-1, SEEK_SET) == -1) {
    fprintf(stderr, "[dbor] could not reserve memory for dbor merged buffer `%s'\n", filename);
    goto fail;
  }
  if (write(dm.fd, "", 1) == -1) {
    fprintf(stderr, "[dbor] could write memory for dbor merged buffer `%s'\n", filename);
    goto fail;
  }
  lseek(dm.fd, 0, SEEK_SET);
  dm.data = mmap(0, dm.data_size, PROT_WRITE | PROT_READ, MAP_SHARED | MAP_NORESERVE, dm.fd, 0);
  memcpy(dm.data, dbor_cascade[1].data, dm.data_size);
  fprintf(stdout, "[dbor] writing `%s'\n", filename);
  
  // get pfm header for checking memcpy
  char *endptr;
  dm.width = strtol(dm.data+3, &endptr, 10);
  dm.height = strtol(endptr, &endptr, 10);
  endptr++; // remove newline
  while(endptr < (char *)dm.data + 100 && *endptr != '\n') endptr++;
  dm.pixel = (float*)(++endptr); // remove second newline
  assert(dm.width == dbor_cascade[0].width);
  assert(dm.height == dbor_cascade[0].height);

  // first fill buffer with lowest cascade
  for (int h = 0; h < dm.height; ++h)
  for (int w = 0; w < dm.width; ++w)
  {
    const int offset = 3*(h*dm.width+w);
    
    for (int c = 0; c < 3; ++c) 
      // XXX DEBUG: switched off so we can see only our blurred extra levels on the output
      // dm.pixel[offset+c] = 0.0f;// dbor_cascade[0].pixel[offset+c];
      dm.pixel[offset+c] = dbor_cascade[0].pixel[offset+c];
  }

  // ============================================
  //   alloc temp memory used for diffusion
  // ============================================

  // diffusion tensor:
  tensor = malloc(sizeof(float)*10*dm.width*dm.height);
  // alloc ping pong buffer for anisotropic diffusion
  buf0 = malloc(sizeof(float)*3*dm.width*dm.height);
  buf1 = malloc(sizeof(float)*3*dm.width*dm.height);

  // alloc temp memory for mip mapping/hole filling:
  int size = fminf(dm.width, dm.height);
  for(mip_max_level=0;mip_max_level<20&&size;mip_max_level++)
    size = (size-1)/2+1;
  mip_buf = malloc(sizeof(float*)*(mip_max_level+1));
  int wd = dm.width, ht = dm.height;
  for(int l=0;l<=mip_max_level;l++)
  {
    mip_buf[l] = malloc(sizeof(float)*3*wd*ht);
    memset(mip_buf[l], 0, sizeof(float)*3*wd*ht);
    wd = (wd-1)/2+1;
    ht = (ht-1)/2+1;
  }

  for (int i = 1; i < num_dbors; ++i)
  {
    fprintf(stderr, "\r|");
    for(int p=0;  p<=i;       p++) fprintf(stderr, "■");
    for(int p=i+1;p<num_dbors;p++) fprintf(stderr, "‒");
    fprintf(stderr, "|");
    // ============================================
    //  count samples for every pixel in
    //  current cascade level
    // ============================================
    {
    const float* const p0 = dbor_cascade[i-1].pixel;
    const float* const p1 = dbor_cascade[i  ].pixel;
    const float* const p2 = dbor_cascade[i+1].pixel; // will only be accessed if valid
    memset(tensor, 0, sizeof(float)*10*dm.width*dm.height);
#pragma omp parallel for default(none) collapse(2) schedule(static) shared(dm,i,tensor,num_dbors)
    for (int h = 1; h < dm.height-1; ++h)
    for (int w = 1; w < dm.width -1; ++w)
    {
      const int ind = h*dm.width+w;
      tensor[10*ind] = 0.0f;
      float cnt = 0.0f;
      float val = 1.f/3.f * (p1[3*ind]+p1[3*ind+1]+p1[3*ind+2]) / (1<<i);
      if(!(val > 0.f)) continue; // no value, no radius
      for(int jj=-1;jj<=1;jj++)
      for(int ii=-1;ii<=1;ii++)
      {
        const int offset = 3*(ind + jj*dm.width + ii);
        if(i == num_dbors-1)
        { // TODO: replace by sample counter
          cnt += 1.f/3.f*(p0[offset+0]+p0[offset+1]+p0[offset+2]) / (1 << (i-1))
              +  1.f/3.f*(p1[offset+0]+p1[offset+1]+p1[offset+2]) / (1 << i);
        }
        else
        {
          cnt += 1.f/3.f*(p0[offset+0]+p0[offset+1]+p0[offset+2]) / (1 << (i-1))
              +  1.f/3.f*(p1[offset+0]+p1[offset+1]+p1[offset+2]) / (1 << i)
              +  1.f/3.f*(p2[offset+0]+p2[offset+1]+p2[offset+2]) / (1 << (i+1));
        }
      }
      // * 4x4 pixel filter support
      tensor[10*ind] = cnt * 9; // we'll convert to radius after mip mapping
      // TODO: remember max radius
    }
    } // end scope for sample counting

    // ============================================
    //   hole filling and sub-1spp estimation:
    // ============================================
#if 0
    // TODO: DEBUG: set 0 count to something silly:
#pragma omp parallel for default(none) collapse(2) schedule(static) shared(dm,tensor)
    for (int h = 1; h < dm.height-1; ++h)
    for (int w = 1; w < dm.width -1; ++w)
    {
      // TODO: do we need bilateral blur here...?
      const int ind = h*dm.width+w;
      if(tensor[10*ind] == 0.0f) tensor[10*ind] = 0.5f;
    }
#endif
    // fill mip map, only propagate non-zero sample counts upwards
    int wd = dm.width, ht = dm.height;
#pragma omp parallel for default(none) collapse(2) schedule(static) shared(dm,tensor,mip_buf)
    for (int h = 0; h < dm.height; ++h)
    for (int w = 0; w < dm.width ; ++w)
    {
      const int ind = h*dm.width+w;
      mip_buf[0][3*ind+0] = tensor[10*ind];
    }
    for(int l=1;l<=mip_max_level;l++)
    {
      const int wdc = (wd-1)/2+1;
      const int htc = (ht-1)/2+1;
#pragma omp parallel for default(none) collapse(2) schedule(static) shared(dm,tensor,mip_buf,l,wd,ht)
      for(int j=0;j<htc;j++)
      {
        for(int i=0;i<wdc;i++)
        {
          for(int c=0;c<3;c++)
          {
            int cnt = 0;
            mip_buf[l][3*(wdc*j+i)+c] = 0.0f; // clear
            if(mip_buf[l-1][3*(wd*2*j+2*i)+c] > 0.0f)
            {
              mip_buf[l][3*(wdc*j+i)+c] += mip_buf[l-1][3*(wd*2*j+2*i)+c];
              cnt++;
            }
            if((2*i+1 < wd) && (mip_buf[l-1][3*(wd*2*j+2*i+1)+c] > 0.0f))
            {
              mip_buf[l][3*(wdc*j+i)+c] += mip_buf[l-1][3*(wd*2*j+2*i+1)+c];
              cnt++;
            }
            if((2*j+1 < ht) && (mip_buf[l-1][3*(wd*(2*j+1)+2*i)+c] > 0.0f))
            {
              mip_buf[l][3*(wdc*j+i)+c] += mip_buf[l-1][3*(wd*(2*j+1)+2*i)+c];
              cnt++;
            }
            if((2*j+1 < ht) && (2*i+1 < wd) && 
               (mip_buf[l-1][3*(wd*(2*j+1)+2*i+1)+c] > 0.0f))
            {
              mip_buf[l][3*(wdc*j+i)+c] = mip_buf[l-1][3*(wd*(2*j+1)+2*i+1)+c];
              cnt++;
            }
            if(cnt) mip_buf[l][3*(wdc*j+i)+c] /= cnt;
          }
        }
      }
      wd = wdc;
      ht = htc;
    }
    // now the push phase where we update the detailed pixels from the coarse mip map:
    // TODO: go through fine levels
    // TODO: if these are zero, fill with avg sample count from above
    // TODO: if these are below some useful threshold? too different from avg? replace, too
    for(int l=mip_max_level;l>0;l--)
    {
      // width and height for level l-1
      wd = dm.width; ht = dm.height;
      for(int ll=0;ll<l-1;ll++)
      {
        wd = (wd-1)/2+1;
        ht = (ht-1)/2+1;
      }
      const int wdc = (wd-1)/2+1; // this corresponds to level l
      // const int htc = (ht-1)/2+1;
#pragma omp parallel for default(none) collapse(2) schedule(static) shared(mip_buf,l,wd,ht)
      for(int j=0;j<ht;j++)
      {
        for(int i=0;i<wd;i++)
        {
          const float *coarse = mip_buf[l] + 3*(wdc*(j/2)+i/2);
          float *fine = mip_buf[l-1] + 3*(wd*j+i);
          // TODO: now the hard part
          for(int c=0;c<3;c++)
          {
            // TODO: need to also increase radius if fine > coarse?
            if(fine[c] <= DBOR_GLINT_THRS) fine[c] = coarse[c];
          }
        }
      }
    }

    // ============================================
    //   precompute diffusion tensor
    // ============================================
    // precompute tensors used for anisotropic diffusion: every pixel stores:
    // 1 float: diffusion radius converted to iteration count
    // 9 float: normalised weights to all neighbours
    // we'll derive these weights from:
    //  image space differences in the currently trusted buffer
    //  half vector derivatives we get from the render (need to dump buffer and average per pixel!)

    // TODO: clamp to actually needed radius
    const int max_rad = DBOR_MAX_RAD;
    float max_radius = 0.0f;

    // "boundary handling" by not handling boundaries here (same later during diffusion):
#pragma omp parallel for default(none) collapse(2) schedule(static) shared(dm,tensor,max_rad,mip_buf)
    for (int h = 1; h < dm.height-1; ++h)
    for (int w = 1; w < dm.width -1; ++w)
    {
      const uint64_t ind = h*dm.width+w;
      float (*kernel)[3] = (float (*)[3])(tensor + 10*ind+5);
#if 0
      // XXX DEBUG set to constant box filter (will result in gaussian blur)
      for(int y=-1;y<=1;y++)
        for(int x=-1;x<=1;x++)
          kernel[y][x] = 1.f/9.f;
#else
      // pixel difference in base buffer
      float sum = 0.0f;
      // XXX use dm.pixel instead once we're done debugging the blur!
      // const float *base = dbor_cascade[i-1].pixel + 3*ind;
      const float *base = dm.pixel + 3*ind;
      for(int y=-1;y<=1;y++) for(int x=-1;x<=1;x++)
      { // determine diffusion tensor by looking at already done base buffer:
        float pdiff = 0.0f;
        for(int c=0;c<3;c++)
        {
          // TODO: may need more robust criterion if base is noisy.
          // TODO: we arguably want to avoid a noisy base though.
          const float d = (base[0] - base[3*(x+dm.width*y)]);
          pdiff += d*d;
        }
        // this bias here has quite an impact on diffusion speed:
        const float bias = DBOR_EDGE_BIAS;
        kernel[y][x] = 1.0f/(bias + sqrtf(pdiff));
        sum += kernel[y][x];
      }
      // normalise:
      for(int y=-1;y<=1;y++) for(int x=-1;x<=1;x++) kernel[y][x] /= sum;
      // check for inf nan 0
      int broken = 0;
      for(int y=-1;y<=1;y++) for(int x=-1;x<=1;x++) if(!(kernel[y][x] == kernel[y][x]) || kernel[y][x] < 0.0f || kernel[y][x] > FLT_MAX) broken = 1;
      if(broken)
      { // this usually means sum close to 0 which means no difference at all so we blur it to hell:
        sum = 1.0f;
        for(int y=-1;y<=1;y++) for(int x=-1;x<=1;x++) kernel[y][x] = 1.f/9.f;
      }
#endif
      // convert count to radius (count is always > 0 because of the prepass above)
      // (rad/1px)^2 = K/cnt
      const float cnt = mip_buf[0][3*ind];// tensor[10*ind];
      tensor[10*ind] = fminf(max_rad, sqrtf(sample_count/cnt));
#pragma omp critical
      {
        max_radius = fmaxf(max_radius, tensor[10*ind]);
      }
    }

#if 1 // DEBUG output radius map
    {
      char filename[1024] = {0};
      snprintf(filename, sizeof(filename), "radius_%02d.pfm", i);
      FILE *f = fopen(filename, "wb");
      fprintf(f, "PF\n%d %d\n-1.0\n", dm.width, dm.height);
      for (int h = 0; h < dm.height; ++h) for (int w = 0; w < dm.width; ++w)
        for(int c=0;c<3;c++)
          fwrite(tensor + 10*(h*dm.width+w), sizeof(float), 1, f);
      fclose(f);
    }
#endif

    // ============================================
    //   actual anisotropic diffusion
    // ============================================

    // clear ping pong buffers
    memset(buf0, 0, sizeof(*buf0)*dm.width*dm.height*3);
    memset(buf1, 0, sizeof(*buf1)*dm.width*dm.height*3);
    const float *src = buf0;
    float *dst = buf1;
    const float *const emit = dbor_cascade[i].pixel;
    // go through all blur radii, start with largest. we'll
    // only insert the energy for diffusion if the source's radius
    // is large enough for the current iteration.
    for(int ii=max_radius;ii>0;ii--)
    { // iteration for anisotropic diffusion
#pragma omp parallel for default(none) collapse(2) schedule(static) shared(dm,tensor,ii,src,dst)
      for (int h = 1; h < dm.height-1; ++h)
      for (int w = 1; w < dm.width-1; ++w)
      {
        const uint64_t ind = h*dm.width+w;

        // kernel computed in prepass:
        const float (*kernel)[3] = (const float (*)[3])(tensor + 10*ind+5);
        // need to clamp or else it will not be picked up as emissive below.
        const int radius = fmaxf(1.0f, fminf(max_rad, tensor[10*ind]+.5f));
        for(int c=0;c<3;c++)
          dst[3*ind+c] = 0;
        for(int y=-1;y<=1;y++)
        for(int x=-1;x<=1;x++)
        { // grab power from neighbours
          for(int c=0;c<3;c++)
            dst[3*ind+c] += src[3*(ind+y*dm.width+x)+c]*kernel[y][x];
        }
        // also insert sources if their radius is large enough so they may
        // emit in this iteration already
        // TODO: gaussian via b-spline would probably insert once blur often 
        if(radius == ii) // TODO: >= and then 1/ii or == and 1 or what?
          for(int c=0;c<3;c++) // note: only the == case is energy conserving
          {
            assert(emit[3*ind+c] == emit[3*ind+c]);
            assert(emit[3*ind+c] < FLT_MAX);
            assert(emit[3*ind+c] >= 0.0f);
            dst[3*ind+c] += emit[3*ind+c];// / ii;
          }
      }
      // swap buffers
      if(src == buf1) { src = buf0; dst = buf1; }
      else            { dst = buf0; src = buf1; }
    }
#if 1 // DEBUG output blurred image buffer
    {
      char filename[1024] = {0};
      snprintf(filename, sizeof(filename), "blurred_%02d.pfm", i);
      FILE *f = fopen(filename, "wb");
      fprintf(f, "PF\n%d %d\n-1.0\n", dm.width, dm.height);
      for (int h = 0; h < dm.height; ++h) for (int w = 0; w < dm.width; ++w)
        fwrite(src + 3*(h*dm.width+w), sizeof(float), 3, f);
      fclose(f);
    }
#endif
    for (int h = 1; h < dm.height-1; ++h)
    for (int w = 1; w < dm.width-1; ++w)
    {
      const uint64_t ind = h*dm.width+w;
      for(int c=0;c<3;c++) // we swapped, so output will be in src
        dm.pixel[3*ind+c] += src[3*ind+c];
    }
    // XXX DEBUG:
    // this will jump to the end, cleanup and munmap, effectively writing the image out.
    // goto fail;
  }
  fprintf(stderr, "\n");

fail:
  free(tensor);
  free(buf0);
  free(buf1);
  for(int l=0;l<=mip_max_level;l++)
    free(mip_buf[l]);
  free(mip_buf);

  for (int i = 0; i < num_dbors; ++i)
  {
    dbor_buffer_t* db = &dbor_cascade[i];
    munmap(db->data, db->data_size);
    close(db->fd);
  }
  munmap(dm.data, dm.data_size);
  close(dm.fd);
}
