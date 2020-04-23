#include "subdiv.h"
#include <stdio.h>
#include <stdlib.h>


int main(int argc, char *argv[])
{
  if(argc < 3)
  {
    fprintf(stderr, "[geo-subdiv] usage: geo-subdiv shape levels   -- where shape.geo is your file\n");
    exit(1);
  }
  const int levels = atol(argv[2]);

  // mesh ping poing buffers
  sd_mesh_t ping, pong;
  sd_mesh_t *in = &ping, *out = &pong;
  if(sd_mesh_init_from_geo(in, argv[1])) exit(2);

  for(int l=0;l<levels;l++)
  {
    sd_mesh_subdiv(in, out);

    // swap meshes for next iteration
    sd_mesh_t *tmp = out;
    out = in;
    in = tmp;
    sd_mesh_cleanup(out);
  }

  fprintf(stderr, "[geo-subdiv] done. writing out geo\n");
  sd_mesh_create_writeout_order(in);
  // write mesh to geo format
  char output[1024];
  snprintf(output, 1024, "%s-subdiv.geo", argv[1]);
  if(sd_mesh_write_geo(in, output))
    fprintf(stderr, "[geo-subdiv] failed to write output `%s'!\n", output);
  fprintf(stderr, "[geo-subdiv] done. writing out vdata\n");
  snprintf(output, 1024, "%s-subdiv.vdata", argv[1]);
  if(sd_mesh_write_vdata(in, output))
    fprintf(stderr, "[geo-subdiv] failed to write output `%s'!\n", output);
  sd_mesh_cleanup(in);
  exit(0);
}
