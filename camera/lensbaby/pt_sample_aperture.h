float pred_x;
float pred_y;
float pred_dx;
float pred_dy;
float sqr_err = FLT_MAX;
for(int k=0;k<5&&sqr_err > 1e-4f;k++)
{
  const float begin_x = x + dist * dx;
  const float begin_y = y + dist * dy;
  const float begin_dx = dx;
  const float begin_dy = dy;
  __attribute__((unused)) const float begin_lambda = lambda;
  pred_x =  + 2.27408e-09  + 9.40161e-08 *begin_dy + 65 *begin_dx + 2.09594e-10 *begin_y + 1 *begin_x + -4.0681e-09 *begin_y*begin_dy + 1.15011e-08 *begin_y*begin_dx + -5.29119e-11 *lens_ipow(begin_y, 2) + -2.65577e-09 *begin_x*begin_dy + -7.65603e-09 *begin_x*begin_dx + 1.33895e-10 *begin_x*begin_y + 5.85864e-11 *lens_ipow(begin_x, 2) + -1.35159e-05 *lens_ipow(begin_dx, 2)*begin_dy + -8.29487e-08 *begin_y*begin_dx*begin_dy + 1.11873e-08 *begin_y*lens_ipow(begin_dx, 2) + -2.59591e-10 *lens_ipow(begin_y, 2)*begin_dy + -1.06425e-09 *lens_ipow(begin_y, 2)*begin_dx + -2.70649e-08 *begin_x*lens_ipow(begin_dy, 2) + -1.47485e-09 *lens_ipow(begin_x, 2)*begin_dy + -2.15172e-11 *lens_ipow(begin_x, 2)*begin_y + -2.45393e-06 *begin_dx*begin_dy*lens_ipow(begin_lambda, 2) + -3.22205e-09 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2) + -9.74692e-07 *begin_x*lens_ipow(begin_dx, 2)*begin_dy + -1.35658e-12 *lens_ipow(begin_x, 3)*begin_y + 0.000234535 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 3) + -3.42176e-08 *lens_ipow(begin_x, 2)*begin_dx*lens_ipow(begin_dy, 2) + -1.07268e-09 *lens_ipow(begin_x, 3)*begin_dx*begin_dy + 2.23548e-12 *lens_ipow(begin_x, 4)*begin_dx + 2.24018e-06 *begin_x*lens_ipow(begin_dy, 5) + 8.64462e-08 *begin_x*begin_dx*lens_ipow(begin_lambda, 4) + 1.79726e-08 *begin_x*begin_y*lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 2) + -3.66991e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 4) + 3.77071e-09 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dy, 4) + -3.34941e-11 *lens_ipow(begin_x, 3)*begin_y*begin_dx*lens_ipow(begin_lambda, 2) + -0.114902 *lens_ipow(begin_dx, 7)*lens_ipow(begin_dy, 3) + 4.82923e-16 *begin_x*lens_ipow(begin_y, 7)*begin_dx*begin_dy + -4.46802e-14 *lens_ipow(begin_x, 6)*begin_y*begin_dx*lens_ipow(begin_dy, 2) + -3.07233e-05 *begin_y*lens_ipow(begin_dx, 4)*lens_ipow(begin_lambda, 6) + 0.000638466 *begin_y*lens_ipow(begin_dx, 10) + -3.68773e-09 *lens_ipow(begin_x, 4)*lens_ipow(begin_dx, 3)*begin_dy*lens_ipow(begin_lambda, 3);
  pred_y =  + -1.6903e-09  + 65 *begin_dy + -2.09602e-08 *begin_dx + 1 *begin_y + -1.16464e-09 *begin_x + 5.45112e-07 *begin_dx*begin_dy + 5.12608e-09 *begin_y*begin_dy + 7.68168e-11 *lens_ipow(begin_y, 2) + -1.31964e-08 *begin_x*begin_dy + -2.16388e-09 *begin_x*begin_dx + -6.60343e-11 *begin_x*begin_y + 3.38607e-08 *begin_y*lens_ipow(begin_dx, 2) + 1.38741e-07 *begin_x*lens_ipow(begin_dy, 2) + -2.99743e-08 *begin_x*begin_dx*begin_dy + -7.11455e-08 *begin_x*lens_ipow(begin_dx, 2) + 1.22323e-09 *begin_x*begin_y*begin_dy + -1.52558e-09 *lens_ipow(begin_x, 2)*begin_dy + 9.62841e-12 *lens_ipow(begin_x, 3) + 2.15305e-08 *begin_y*begin_dx*lens_ipow(begin_lambda, 2) + 3.72342e-07 *begin_x*lens_ipow(begin_dy, 3) + 2.53751e-09 *begin_x*begin_y*lens_ipow(begin_dx, 2) + 5.08928e-09 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 2) + -6.15124e-13 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2) + 1.48272e-05 *lens_ipow(begin_dy, 5) + 8.65769e-05 *begin_dx*lens_ipow(begin_dy, 4) + 1.45594e-06 *begin_y*lens_ipow(begin_dx, 3)*begin_dy + -4.58885e-12 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*begin_dx + -4.88226e-14 *lens_ipow(begin_x, 4)*begin_y + -6.54416e-06 *begin_x*lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 2) + -1.98809e-08 *begin_x*begin_y*begin_dx*begin_dy*lens_ipow(begin_lambda, 2) + -2.86538e-13 *begin_x*lens_ipow(begin_y, 4)*begin_dy + -0.000279357 *lens_ipow(begin_dx, 7) + 2.24792e-08 *lens_ipow(begin_x, 2)*begin_y*begin_dx*lens_ipow(begin_dy, 3) + -0.000310699 *lens_ipow(begin_dy, 6)*lens_ipow(begin_lambda, 2) + 6.35914e-10 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 3) + -1.24787e-14 *lens_ipow(begin_x, 7)*lens_ipow(begin_dy, 2) + -0.124796 *lens_ipow(begin_dx, 7)*lens_ipow(begin_dy, 3) + -3.89281e-09 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 5) + -0.000459617 *begin_dx*lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 7) + -4.09566e-06 *lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 5)*lens_ipow(begin_dy, 3);
  pred_dx =  + 6.00192e-11  + 2.50376e-10 *begin_dy + 1 *begin_dx + 1.03502e-11 *begin_y + -2.40826e-12 *begin_x + -5.73958e-09 *begin_dx*begin_dy + -3.62406e-09 *lens_ipow(begin_dx, 2) + 5.78082e-11 *begin_y*begin_dx + -1.36359e-14 *lens_ipow(begin_y, 2) + -3.65037e-13 *lens_ipow(begin_x, 2) + -9.73885e-09 *begin_dx*lens_ipow(begin_dy, 2) + 8.15663e-11 *begin_x*lens_ipow(begin_dx, 2) + -4.87636e-12 *begin_x*begin_y*begin_dy + 2.36671e-11 *begin_x*begin_y*begin_dx + 4.40988e-07 *lens_ipow(begin_dx, 3)*begin_dy + -4.21174e-11 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2) + 1.45922e-11 *begin_x*begin_y*lens_ipow(begin_dx, 2) + 4.03979e-11 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 2) + 1.26099e-13 *lens_ipow(begin_x, 2)*begin_y*begin_dy + -1.15912e-07 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 3) + 1.23141e-10 *begin_x*begin_y*lens_ipow(begin_dx, 2)*begin_dy + -6.83272e-14 *begin_x*lens_ipow(begin_y, 3)*begin_dx + -1.1099e-10 *begin_y*lens_ipow(begin_lambda, 5) + 4.66621e-08 *begin_x*lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 2) + -4.7903e-15 *lens_ipow(begin_x, 5)*lens_ipow(begin_dy, 2) + 2.07788e-11 *begin_x*lens_ipow(begin_y, 2)*begin_dx*begin_dy*lens_ipow(begin_lambda, 3) + -1.32945e-08 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 6)*begin_lambda + -4.1643e-17 *lens_ipow(begin_y, 7)*begin_dx*begin_dy + -1.65258e-06 *begin_x*begin_dx*lens_ipow(begin_dy, 7) + -2.94633e-05 *begin_y*lens_ipow(begin_dx, 6)*lens_ipow(begin_dy, 3) + 3.4978e-18 *lens_ipow(begin_y, 8)*lens_ipow(begin_dx, 2) + 3.37778e-18 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 6)*begin_dx*begin_dy + 9.27307e-10 *lens_ipow(begin_x, 3)*lens_ipow(begin_dx, 4)*begin_dy*lens_ipow(begin_lambda, 2) + -0.000611784 *lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 6)*lens_ipow(begin_lambda, 2) + -0.00220415 *lens_ipow(begin_dx, 7)*lens_ipow(begin_dy, 3)*begin_lambda + 6.00972e-05 *lens_ipow(begin_dx, 9)*lens_ipow(begin_lambda, 2) + 1.08132e-06 *begin_y*lens_ipow(begin_dx, 6)*lens_ipow(begin_lambda, 4) + -3.23839e-08 *lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 6)*lens_ipow(begin_dy, 2) + -2.32767e-10 *begin_x*lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 6)*begin_dy + -1.37479e-17 *lens_ipow(begin_x, 7)*begin_y*begin_dx*lens_ipow(begin_dy, 2);
  pred_dy =  + 7.45533e-11  + 1 *begin_dy + -1.14801e-10 *begin_dx + -1.71562e-12 *begin_y + -1.38706e-12 *begin_x + -5.2179e-09 *lens_ipow(begin_dy, 2) + -4.05273e-09 *begin_dx*begin_dy + -2.86197e-09 *lens_ipow(begin_dx, 2) + -3.96147e-13 *lens_ipow(begin_y, 2) + -2.66841e-11 *begin_x*begin_dy + -6.40023e-13 *lens_ipow(begin_x, 2) + -4.58721e-09 *lens_ipow(begin_dx, 2)*begin_dy + 6.93122e-10 *begin_x*lens_ipow(begin_dy, 2) + -8.46011e-11 *begin_x*begin_dx*begin_lambda + 4.39597e-12 *begin_x*begin_y*begin_dx + 4.56749e-14 *lens_ipow(begin_x, 2)*begin_y + 1.59617e-07 *lens_ipow(begin_dx, 3)*begin_dy + -1.27506e-09 *begin_y*begin_dx*begin_dy*begin_lambda + -1.23964e-11 *begin_x*begin_y*begin_dx*begin_dy + 7.03863e-11 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 2) + -2.23676e-15 *lens_ipow(begin_x, 3)*begin_y + 1.0605e-10 *lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 3) + 3.45635e-11 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*begin_lambda + 2.04907e-10 *begin_x*begin_y*lens_ipow(begin_dy, 3) + -2.90979e-10 *begin_x*begin_y*begin_dx*lens_ipow(begin_dy, 2) + -4.21009e-10 *begin_x*begin_y*lens_ipow(begin_dx, 2)*begin_dy + 9.5847e-07 *lens_ipow(begin_dy, 6) + 1.54819e-09 *begin_y*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 3) + 1.47964e-08 *begin_y*lens_ipow(begin_dx, 5) + 1.23202e-07 *begin_x*begin_dx*lens_ipow(begin_dy, 4) + -9.24124e-15 *lens_ipow(begin_x, 5)*lens_ipow(begin_dy, 2) + -2.12014e-06 *begin_y*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 5) + 3.15581e-08 *begin_x*lens_ipow(begin_dy, 5)*lens_ipow(begin_lambda, 2) + 4.93641e-10 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 2) + -8.48157e-18 *lens_ipow(begin_x, 6)*begin_y*begin_dx + -8.28261e-09 *lens_ipow(begin_x, 2)*begin_dx*lens_ipow(begin_dy, 4)*lens_ipow(begin_lambda, 2) + 5.01729e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 5)*lens_ipow(begin_dy, 2) + -6.22663e-09 *lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 4) + 2.83951e-13 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 3) + -1.81622e-19 *lens_ipow(begin_x, 6)*lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 2);
  float dx1_domega0[2][2];
  dx1_domega0[0][0] =  + 65  + 1.15011e-08 *begin_y + -7.65603e-09 *begin_x + -2.70318e-05 *begin_dx*begin_dy + -8.29487e-08 *begin_y*begin_dy + 2.23746e-08 *begin_y*begin_dx + -1.06425e-09 *lens_ipow(begin_y, 2) + -2.45393e-06 *begin_dy*lens_ipow(begin_lambda, 2) + -6.44409e-09 *lens_ipow(begin_y, 2)*begin_dx + -1.94938e-06 *begin_x*begin_dx*begin_dy + 0.000469069 *begin_dx*lens_ipow(begin_dy, 3) + -3.42176e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 2) + -1.07268e-09 *lens_ipow(begin_x, 3)*begin_dy + 2.23548e-12 *lens_ipow(begin_x, 4) + 8.64462e-08 *begin_x*lens_ipow(begin_lambda, 4) + 3.59452e-08 *begin_x*begin_y*begin_dx*lens_ipow(begin_lambda, 2) + -3.34941e-11 *lens_ipow(begin_x, 3)*begin_y*lens_ipow(begin_lambda, 2) + -0.804311 *lens_ipow(begin_dx, 6)*lens_ipow(begin_dy, 3) + 4.82923e-16 *begin_x*lens_ipow(begin_y, 7)*begin_dy + -4.46802e-14 *lens_ipow(begin_x, 6)*begin_y*lens_ipow(begin_dy, 2) + -0.000122893 *begin_y*lens_ipow(begin_dx, 3)*lens_ipow(begin_lambda, 6) + 0.00638466 *begin_y*lens_ipow(begin_dx, 9) + -1.10632e-08 *lens_ipow(begin_x, 4)*lens_ipow(begin_dx, 2)*begin_dy*lens_ipow(begin_lambda, 3)+0.0f;
  dx1_domega0[0][1] =  + 9.40161e-08  + -4.0681e-09 *begin_y + -2.65577e-09 *begin_x + -1.35159e-05 *lens_ipow(begin_dx, 2) + -8.29487e-08 *begin_y*begin_dx + -2.59591e-10 *lens_ipow(begin_y, 2) + -5.41298e-08 *begin_x*begin_dy + -1.47485e-09 *lens_ipow(begin_x, 2) + -2.45393e-06 *begin_dx*lens_ipow(begin_lambda, 2) + -9.74692e-07 *begin_x*lens_ipow(begin_dx, 2) + 0.000703604 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + -6.84352e-08 *lens_ipow(begin_x, 2)*begin_dx*begin_dy + -1.07268e-09 *lens_ipow(begin_x, 3)*begin_dx + 1.12009e-05 *begin_x*lens_ipow(begin_dy, 4) + -1.46796e-07 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 3) + 1.50828e-08 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dy, 3) + -0.344705 *lens_ipow(begin_dx, 7)*lens_ipow(begin_dy, 2) + 4.82923e-16 *begin_x*lens_ipow(begin_y, 7)*begin_dx + -8.93603e-14 *lens_ipow(begin_x, 6)*begin_y*begin_dx*begin_dy + -3.68773e-09 *lens_ipow(begin_x, 4)*lens_ipow(begin_dx, 3)*lens_ipow(begin_lambda, 3)+0.0f;
  dx1_domega0[1][0] =  + -2.09602e-08  + 5.45112e-07 *begin_dy + -2.16388e-09 *begin_x + 6.77214e-08 *begin_y*begin_dx + -2.99743e-08 *begin_x*begin_dy + -1.42291e-07 *begin_x*begin_dx + 2.15305e-08 *begin_y*lens_ipow(begin_lambda, 2) + 5.07502e-09 *begin_x*begin_y*begin_dx + 8.65769e-05 *lens_ipow(begin_dy, 4) + 4.36781e-06 *begin_y*lens_ipow(begin_dx, 2)*begin_dy + -4.58885e-12 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2) + -1.96325e-05 *begin_x*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + -1.98809e-08 *begin_x*begin_y*begin_dy*lens_ipow(begin_lambda, 2) + -0.0019555 *lens_ipow(begin_dx, 6) + 2.24792e-08 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dy, 3) + 6.35914e-10 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 3) + -0.873574 *lens_ipow(begin_dx, 6)*lens_ipow(begin_dy, 3) + -0.000459617 *lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 7) + -2.04783e-05 *lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 3)+0.0f;
  dx1_domega0[1][1] =  + 65  + 5.45112e-07 *begin_dx + 5.12608e-09 *begin_y + -1.31964e-08 *begin_x + 2.77483e-07 *begin_x*begin_dy + -2.99743e-08 *begin_x*begin_dx + 1.22323e-09 *begin_x*begin_y + -1.52558e-09 *lens_ipow(begin_x, 2) + 1.11703e-06 *begin_x*lens_ipow(begin_dy, 2) + 1.01786e-08 *lens_ipow(begin_x, 2)*begin_dy + 7.41358e-05 *lens_ipow(begin_dy, 4) + 0.000346308 *begin_dx*lens_ipow(begin_dy, 3) + 1.45594e-06 *begin_y*lens_ipow(begin_dx, 3) + -1.30883e-05 *begin_x*lens_ipow(begin_dx, 3)*begin_dy + -1.98809e-08 *begin_x*begin_y*begin_dx*lens_ipow(begin_lambda, 2) + -2.86538e-13 *begin_x*lens_ipow(begin_y, 4) + 6.74377e-08 *lens_ipow(begin_x, 2)*begin_y*begin_dx*lens_ipow(begin_dy, 2) + -0.0018642 *lens_ipow(begin_dy, 5)*lens_ipow(begin_lambda, 2) + 1.90774e-09 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 2) + -2.49574e-14 *lens_ipow(begin_x, 7)*begin_dy + -0.374389 *lens_ipow(begin_dx, 7)*lens_ipow(begin_dy, 2) + -7.78562e-09 *begin_x*lens_ipow(begin_y, 2)*begin_dy*lens_ipow(begin_lambda, 5) + -0.00137885 *begin_dx*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 7) + -1.2287e-05 *lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 5)*lens_ipow(begin_dy, 2)+0.0f;
  float invJ[2][2];
  const float invdet = 1.0f/(dx1_domega0[0][0]*dx1_domega0[1][1] - dx1_domega0[0][1]*dx1_domega0[1][0]);
  invJ[0][0] =  dx1_domega0[1][1]*invdet;
  invJ[1][1] =  dx1_domega0[0][0]*invdet;
  invJ[0][1] = -dx1_domega0[0][1]*invdet;
  invJ[1][0] = -dx1_domega0[1][0]*invdet;
  const float dx1[2] = {out_x - pred_x, out_y - pred_y};
  for(int i=0;i<2;i++)
  {
    dx += invJ[0][i]*dx1[i];
    dy += invJ[1][i]*dx1[i];
  }
  sqr_err = dx1[0]*dx1[0] + dx1[1]*dx1[1];
}
out_dx = pred_dx;
out_dy = pred_dy;