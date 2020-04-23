// this implements the academy colourspace aces,
// with xy primaries 0.73470 0.26530, 0.0 1.0, 0.00010 -0.07700,
// neutral at equal rgb coordinates or x=0.32168, y=0.33767.

// convert input color to xyz
static inline void colour_aces_to_xyz(const float *const input, float *xyz)
{
  const float RGBtoXYZ[] =
  {
    0.95255239, 0.00000000,  0.00009367,
    0.34396644, 0.72816609, -0.07213254,
    0.00000000, 0.00000000,  1.00882518,
  };
  mat3_mulv(RGBtoXYZ, input, xyz);
}

// convert xyz to input color
static inline void colour_xyz_to_aces(const float *const xyz, float *rgb)
{
  const float XYZtoRGB[] =
  {
    1.04981101, 0.00000000, -0.00009748,
   -0.49590302, 1.37331304,  0.09824003,
    0.00000000, 0.00000000,  0.99125201,
  };
  mat3_mulv(XYZtoRGB, xyz, rgb);
}

static inline void colour_aces_print_info(FILE *f)
{
  fprintf(f, "academy color encoding specification (aces)\n");
}

