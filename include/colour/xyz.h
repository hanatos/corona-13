// xyz colour space. only identity transforms
static inline void colour_xyz_to_xyz(const float *const xyz, float *out)
{
  for(int k=0;k<3;k++) out[k] = xyz[k];
}

static inline void colour_xyz_print_info(FILE *f)
{
  fprintf(f, "CIE XYZ\n");
}
