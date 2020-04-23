#include <math.h>
static inline void colour_xyz_to_adobergb(const float *const xyz, float *rgb)
{
  // see http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
  const float XYZtoRGB[] =
  {
    // adobe rgb 1998
    2.0413690, -0.5649464, -0.3446944,
   -0.9692660,  1.8760108,  0.0415560,
    0.0134474, -0.1183897,  1.0154096,
  };

  rgb[0] = rgb[1] = rgb[2] = 0.0f;
  for(int k=0;k<3;k++)
    for(int i=0;i<3;i++) rgb[k] += xyz[i]*XYZtoRGB[i+3*k];
  
  // apply tonecurve. adobe rgb does not have a linear toe slope, but gamma of:
  const float g = 1.f/2.19921875f;
  for(int k=0;k<3;k++) rgb[k] = powf(rgb[k], g);
}

static inline void colour_adobergb_to_xyz(const float *const rgb, float *xyz)
{
  const float RGBtoXYZ[] =
  {
    0.5767309, 0.1855540, 0.1881852,
    0.2973769, 0.6273491, 0.0752741,
    0.0270343, 0.0706872, 0.9911085,
  };

  float linear[3];
  // undo tonecurve. adobe rgb does not have a linear toe slope, but gamma of:
  const float g = 2.19921875f;
  for(int k=0;k<3;k++) linear[k] = powf(rgb[k], g);

  xyz[0] = xyz[1] = xyz[2] = 0.0f;
  for(int k=0;k<3;k++)
    for(int i=0;i<3;i++) xyz[k] += linear[i]*RGBtoXYZ[i+3*k];
}

static inline void colour_adobergb_print_info(FILE *f)
{
  fprintf(f, "adobergb 1998\n");
}
