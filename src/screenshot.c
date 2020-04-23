#include "screenshot.h"
#include "corona_common.h"
#include <unistd.h>
#include <stdio.h>
#include <string.h>

void screenshot_write(
    const char *basename,
    const float *const pixel,
    const int bpp,
    const int offs,
    const int channels,
    uint64_t width,
    uint64_t height,
    float scale)
{
  char filename[512];
  snprintf(filename, sizeof(filename), "%s.pfm", basename);
  FILE* f = fopen(filename, "wb");
  if(f)
  {
    // align pfm header to sse, assuming the file will
    // be mmapped to page boundaries.
    char header[1024];
    snprintf(header, 1024, "PF\n%lu %lu\n-1.0", width, height);
    size_t len = strlen(header);
    fprintf(f, "PF\n%lu %lu\n-1.0", width, height);
    ssize_t off = 0;
    while((len + 1 + off) & 0xf) off++;
    while(off-- > 0) fprintf(f, "0");
    fprintf(f, "\n");
    for(int64_t j=0;j<height;j++)
    {
      for(int64_t i=0;i<width;i++)
      {
        int64_t p = bpp*((int)i+width*(int)j) + offs;
        float val[3] = {
          MAX(0.f, pixel[p + (0%channels)] * scale),
          MAX(0.f, pixel[p + (1%channels)] * scale),
          MAX(0.f, pixel[p + (2%channels)] * scale)
        };
        fwrite(val, sizeof(float), 3, f);
      }
    }
    fclose(f);
  }
}
