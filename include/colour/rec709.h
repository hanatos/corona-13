static inline void colour_xyz_to_rec709(const float *const xyz, float *rgb)
{
  // see http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
  const float XYZtoRGB[] =
  {
    3.2404542, -1.5371385, -0.4985314,
   -0.9692660,  1.8760108,  0.0415560,
    0.0556434, -0.2040259,  1.0572252,
  };
  mat3_mulv(XYZtoRGB, xyz, rgb);
}

static inline void colour_rec709_to_xyz(const float *const rgb, float *xyz)
{
  const float RGBtoXYZ[] =
  {
    0.4124564, 0.3575761, 0.1804375,
    0.2126729, 0.7151522, 0.0721750,
    0.0193339, 0.1191920, 0.9503041,
  };
  mat3_mulv(RGBtoXYZ, rgb, xyz);
}

static inline void colour_rec709_print_info(FILE *f)
{
  fprintf(f, "linear rec709 D65\n");
}

