// this implements the colourspace used by brian smits, i.e the rec709 (sRGB) primaries
// adapted to illuminant E (equal energy white, good idea for reflectances).
// note that these are whitepoint adapted by simple scaling, bradford adaptation messes
// with the primaries and was not the base of brian's work and will thus give wrong results.

// convert ergb to xyz
static inline void colour_ergb_to_xyz(const float *const ergb, float *xyz)
{
  const float RGBtoXYZ[] =
  {
    0.496859,  0.339094,  0.164047,
    0.256193,  0.678188,  0.065619,
    0.023290,  0.113031,  0.863978,
  };
  mat3_mulv(RGBtoXYZ, ergb, xyz);
}

// convert xyz to ergb
static inline void colour_xyz_to_ergb(const float *const xyz, float *ergb)
{
  const float XYZtoRGB[] =
  {
    2.689989, -1.276020, -0.413844,
    -1.022095,  1.978261,  0.043821,
    0.061203, -0.224411,  1.162859,
  };
  mat3_mulv(XYZtoRGB, xyz, ergb);
}

static inline void colour_ergb_print_info(FILE *f)
{
  fprintf(f, "linear rec709 adapted to illuminant E\n");
}
