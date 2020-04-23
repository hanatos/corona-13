static inline void colour_xyz_to_srgb(const float *const xyz, float *rgb)
{
  // see http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
  const float XYZtoRGB[] =
  {
    // sRGB D65
    3.2404542, -1.5371385, -0.4985314,
   -0.9692660,  1.8760108,  0.0415560,
    0.0556434, -0.2040259,  1.0572252,
  };

  memset(rgb, 0, sizeof(float)*3);
  for(int k=0;k<3;k++)
    for(int i=0;i<3;i++) rgb[k] += xyz[i]*XYZtoRGB[i+3*k];

  // add srgb gamma with linear toe slope:
  for(int k=0;k<3;k++)
    rgb[k] = rgb[k] <= 0.0031308f ? 12.92f * rgb[k] : 1.055f * powf(rgb[k], 1.f/2.4f) - 0.055f;
}

static inline void colour_srgb_to_xyz(const float *const rgb, float *xyz)
{
  const float RGBtoXYZ[] =
  {
    0.4124564, 0.3575761, 0.1804375,
    0.2126729, 0.7151522, 0.0721750,
    0.0193339, 0.1191920, 0.9503041,
  };

  float linear[3];
  // undo tonecurve
  for(int k=0;k<3;k++)
    linear[k] = (rgb[k] <= 0.04045) ? rgb[k] / 12.92f : powf((rgb[k] + 0.055f)/1.055f, 2.4f);

  xyz[0] = xyz[1] = xyz[2] = 0.0f;
  for(int k=0;k<3;k++)
    for(int i=0;i<3;i++) xyz[k] += linear[i]*RGBtoXYZ[i+3*k];
}

static inline void colour_srgb_print_info(FILE *f)
{
  fprintf(f, "sRGB\n");
}
