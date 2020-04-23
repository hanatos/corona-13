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
  pred_x =  + -1.45977e-05  + 33.0198 *begin_dx + 4.34706e-06 *begin_y + 0.590865 *begin_x + -0.00612442 *lens_ipow(begin_dy, 2) + 0.000107017 *begin_y*begin_dx + 3.60094e-06 *lens_ipow(begin_y, 2) + 6.94731e-06 *begin_x*begin_y + 2.39929e-06 *lens_ipow(begin_x, 2) + 28.7241 *begin_dx*lens_ipow(begin_dy, 2) + 28.7695 *lens_ipow(begin_dx, 3) + 2.72925 *begin_y*begin_dx*begin_dy + 0.0443946 *lens_ipow(begin_y, 2)*begin_dx + 1.03501 *begin_x*lens_ipow(begin_dy, 2) + -0.000629516 *begin_x*begin_dx*begin_dy + 3.78078 *begin_x*lens_ipow(begin_dx, 2) + 0.0773786 *begin_x*begin_y*begin_dy + 0.00113911 *begin_x*lens_ipow(begin_y, 2) + 0.121896 *lens_ipow(begin_x, 2)*begin_dx + 5.35929e-07 *lens_ipow(begin_x, 2)*begin_y + 0.00112678 *lens_ipow(begin_x, 3) + -0.0121689 *begin_y*lens_ipow(begin_dx, 2)*begin_dy + 1.03302e-05 *begin_x*lens_ipow(begin_y, 2)*begin_dx + -0.075106 *begin_x*lens_ipow(begin_lambda, 4) + 0.112328 *begin_x*lens_ipow(begin_dy, 4) + -1.89372e-05 *lens_ipow(begin_x, 2)*begin_y*begin_dx*begin_dy + 0.000288058 *lens_ipow(begin_x, 3)*lens_ipow(begin_dy, 2) + 1.28042e-07 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 2) + 2.15719e-07 *lens_ipow(begin_x, 5) + -2.21332 *begin_dx*lens_ipow(begin_lambda, 5) + -9.00309e-07 *lens_ipow(begin_x, 3)*begin_y*lens_ipow(begin_lambda, 2) + 8.0143e-08 *lens_ipow(begin_x, 6)*begin_dx + 4.365e-07 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 3)*begin_dy*begin_lambda + 4.29905e-12 *begin_x*lens_ipow(begin_y, 8) + -9.78558e-08 *lens_ipow(begin_x, 6)*begin_y*begin_dx*begin_dy + 3.86969e-13 *lens_ipow(begin_x, 9)*begin_y + 8.87442 *begin_dx*lens_ipow(begin_lambda, 10) + 0.359032 *begin_x*lens_ipow(begin_lambda, 10) + 7.74365e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 6)*begin_dx*lens_ipow(begin_dy, 2) + 6.52176e-13 *lens_ipow(begin_x, 5)*lens_ipow(begin_y, 6);
  pred_y =  + -1.27715e-05  + 33.028 *begin_dy + 0.592062 *begin_y + 5.99218e-06 *begin_x + 0.00200405 *lens_ipow(begin_dx, 2) + 0.000153303 *begin_y*begin_dx + 3.71186e-06 *begin_x*begin_y + 28.6036 *lens_ipow(begin_dy, 3) + 28.8268 *lens_ipow(begin_dx, 2)*begin_dy + 3.76492 *begin_y*lens_ipow(begin_dy, 2) + 1.04157 *begin_y*lens_ipow(begin_dx, 2) + 0.121886 *lens_ipow(begin_y, 2)*begin_dy + 0.00113931 *lens_ipow(begin_y, 3) + 2.73214 *begin_x*begin_dx*begin_dy + 2.9143e-05 *begin_x*begin_y*begin_dy + 0.0776132 *begin_x*begin_y*begin_dx + 9.04475e-07 *begin_x*lens_ipow(begin_y, 2) + 0.0445754 *lens_ipow(begin_x, 2)*begin_dy + 0.00114321 *lens_ipow(begin_x, 2)*begin_y + -0.0476118 *begin_y*lens_ipow(begin_lambda, 3) + -1.26068 *begin_dy*lens_ipow(begin_lambda, 4) + 0.0873411 *begin_y*lens_ipow(begin_dx, 4) + 0.000142096 *lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 2) + 2.19267e-08 *lens_ipow(begin_y, 5) + -0.000251576 *begin_x*lens_ipow(begin_y, 2)*begin_dx*begin_dy + -5.49403e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 3) + 7.58404e-06 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 2) + -0.000390339 *lens_ipow(begin_x, 3)*begin_y*begin_dx*lens_ipow(begin_dy, 2) + -6.75449e-05 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*begin_dy*lens_ipow(begin_lambda, 4) + 0.000799558 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + 2.39644e-07 *lens_ipow(begin_x, 5)*begin_y*begin_dx*lens_ipow(begin_lambda, 2) + 1.0824e-08 *lens_ipow(begin_x, 7)*begin_dx*begin_dy + 5.69878e-12 *lens_ipow(begin_x, 8)*begin_y + 1.03377e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 6)*begin_dy*begin_lambda + 5.96074 *begin_dy*lens_ipow(begin_lambda, 10) + 0.246937 *begin_y*lens_ipow(begin_lambda, 10) + 6.1849e-13 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 9) + -6.52784e-09 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 6)*begin_dx*begin_dy + 0.000742958 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 3) + 7.37959e-13 *lens_ipow(begin_x, 6)*lens_ipow(begin_y, 5);
  pred_dx =  + 2.14133e-06  + -0.0430707 *begin_dx + -3.91902e-07 *begin_y + -0.0310928 *begin_x + -3.21064 *begin_dx*lens_ipow(begin_dy, 2) + -3.2119 *lens_ipow(begin_dx, 3) + -0.0544335 *begin_y*begin_dx*begin_dy + 0.000181522 *lens_ipow(begin_y, 2)*begin_dx + -0.0223645 *begin_x*lens_ipow(begin_dy, 2) + -0.076129 *begin_x*lens_ipow(begin_dx, 2) + 0.00119072 *begin_x*begin_y*begin_dy + 2.75323e-05 *begin_x*lens_ipow(begin_y, 2) + 0.00140174 *lens_ipow(begin_x, 2)*begin_dx + 6.68032e-08 *lens_ipow(begin_x, 2)*begin_y + 2.77685e-05 *lens_ipow(begin_x, 3) + -0.182264 *begin_dx*lens_ipow(begin_lambda, 3) + -0.00643498 *begin_x*lens_ipow(begin_lambda, 4) + 0.000988239 *begin_x*begin_y*lens_ipow(begin_dx, 2)*begin_dy + -0.0217416 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 5) + 0.314136 *begin_x*lens_ipow(begin_dy, 6) + 7.31978e-09 *begin_x*lens_ipow(begin_y, 5)*begin_dy + 4.03418e-08 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 2)*begin_dx + -0.00915566 *lens_ipow(begin_x, 2)*begin_dx*lens_ipow(begin_dy, 4)*begin_lambda + 8.97312e-08 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 3)*begin_dy*begin_lambda + 4.56841e-08 *lens_ipow(begin_x, 5)*lens_ipow(begin_dy, 2)*begin_lambda + 73.9822 *lens_ipow(begin_dx, 9) + 0.0244948 *lens_ipow(begin_x, 3)*lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 2) + 1.37942e-11 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 6) + -0.000397675 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 3) + -0.00036476 *lens_ipow(begin_x, 4)*lens_ipow(begin_dx, 5)*begin_lambda + 0.975594 *begin_dx*lens_ipow(begin_lambda, 10) + 1901 *lens_ipow(begin_dx, 7)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 2) + -3.14508 *begin_y*begin_dx*lens_ipow(begin_dy, 5)*lens_ipow(begin_lambda, 4) + 0.0317773 *begin_x*lens_ipow(begin_lambda, 10) + 3.73242e-12 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 8)*begin_dx + -3.36884e-10 *lens_ipow(begin_x, 7)*lens_ipow(begin_lambda, 4) + -1.09166e-10 *lens_ipow(begin_x, 8)*begin_y*begin_dx*begin_dy + -6.2297e-11 *lens_ipow(begin_x, 9)*lens_ipow(begin_dx, 2) + 6.41592e-14 *lens_ipow(begin_x, 9)*lens_ipow(begin_y, 2) + 8.05435e-15 *lens_ipow(begin_x, 11);
  pred_dy =  + 4.66018e-07  + -0.0438014 *begin_dy + -0.031091 *begin_y + 3.38804e-07 *begin_x + 6.15983e-06 *begin_y*begin_dx + -3.18856 *lens_ipow(begin_dy, 3) + -3.19362 *lens_ipow(begin_dx, 2)*begin_dy + -0.0762984 *begin_y*lens_ipow(begin_dy, 2) + -0.0224616 *begin_y*lens_ipow(begin_dx, 2) + 0.00138872 *lens_ipow(begin_y, 2)*begin_dy + 2.74592e-05 *lens_ipow(begin_y, 3) + -0.0539869 *begin_x*begin_dx*begin_dy + 0.00120299 *begin_x*begin_y*begin_dx + 0.000197944 *lens_ipow(begin_x, 2)*begin_dy + 2.76321e-05 *lens_ipow(begin_x, 2)*begin_y + -0.181453 *begin_dy*lens_ipow(begin_lambda, 3) + -0.00634061 *begin_y*lens_ipow(begin_lambda, 4) + 0.00101406 *begin_x*begin_y*begin_dx*lens_ipow(begin_dy, 2) + 9.78031e-10 *begin_x*lens_ipow(begin_y, 4) + 0.441651 *begin_y*lens_ipow(begin_dx, 6) + -0.0751138 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 3) + -0.0197984 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 5) + -4.66264e-07 *lens_ipow(begin_x, 4)*begin_y*lens_ipow(begin_dx, 2) + -5.034e-05 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 3)*begin_dy*begin_lambda + 5.80044e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 4)*begin_dy*begin_lambda + -0.000199781 *lens_ipow(begin_y, 4)*lens_ipow(begin_dy, 5) + -6.17192e-09 *lens_ipow(begin_y, 7)*lens_ipow(begin_dy, 2) + -5.60323e-09 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 5)*lens_ipow(begin_dx, 2) + -0.000442067 *lens_ipow(begin_x, 4)*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 3) + -3.10694e-08 *lens_ipow(begin_x, 5)*lens_ipow(begin_y, 2)*begin_dx*begin_dy + 8.9558e-12 *lens_ipow(begin_x, 6)*lens_ipow(begin_y, 3) + 1.00562e-07 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 3)*begin_dx*lens_ipow(begin_lambda, 3) + 0.966882 *begin_dy*lens_ipow(begin_lambda, 10) + 12347.6 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 9) + 2639.39 *lens_ipow(begin_dx, 6)*lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 2) + 0.0306743 *begin_y*lens_ipow(begin_lambda, 10) + 5.37018e-15 *lens_ipow(begin_y, 11) + -1.09944e-10 *begin_x*lens_ipow(begin_y, 8)*begin_dx*begin_dy + -0.000674305 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*begin_dy*lens_ipow(begin_lambda, 4) + 6.14979e-14 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 9);
  float dx1_domega0[2][2];
  dx1_domega0[0][0] =  + 33.0198  + 0.000107017 *begin_y + 28.7241 *lens_ipow(begin_dy, 2) + 86.3086 *lens_ipow(begin_dx, 2) + 2.72925 *begin_y*begin_dy + 0.0443946 *lens_ipow(begin_y, 2) + -0.000629516 *begin_x*begin_dy + 7.56156 *begin_x*begin_dx + 0.121896 *lens_ipow(begin_x, 2) + -0.0243377 *begin_y*begin_dx*begin_dy + 1.03302e-05 *begin_x*lens_ipow(begin_y, 2) + -1.89372e-05 *lens_ipow(begin_x, 2)*begin_y*begin_dy + -2.21332 *lens_ipow(begin_lambda, 5) + 8.0143e-08 *lens_ipow(begin_x, 6) + -9.78558e-08 *lens_ipow(begin_x, 6)*begin_y*begin_dy + 8.87442 *lens_ipow(begin_lambda, 10) + 7.74365e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 6)*lens_ipow(begin_dy, 2)+0.0f;
  dx1_domega0[0][1] =  + -0.0122488 *begin_dy + 57.4482 *begin_dx*begin_dy + 2.72925 *begin_y*begin_dx + 2.07003 *begin_x*begin_dy + -0.000629516 *begin_x*begin_dx + 0.0773786 *begin_x*begin_y + -0.0121689 *begin_y*lens_ipow(begin_dx, 2) + 0.449312 *begin_x*lens_ipow(begin_dy, 3) + -1.89372e-05 *lens_ipow(begin_x, 2)*begin_y*begin_dx + 0.000576116 *lens_ipow(begin_x, 3)*begin_dy + 4.365e-07 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 3)*begin_lambda + -9.78558e-08 *lens_ipow(begin_x, 6)*begin_y*begin_dx + 1.54873e-07 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 6)*begin_dx*begin_dy+0.0f;
  dx1_domega0[1][0] =  + 0.00400811 *begin_dx + 0.000153303 *begin_y + 57.6535 *begin_dx*begin_dy + 2.08315 *begin_y*begin_dx + 2.73214 *begin_x*begin_dy + 0.0776132 *begin_x*begin_y + 0.349364 *begin_y*lens_ipow(begin_dx, 3) + 0.000284192 *lens_ipow(begin_y, 3)*begin_dx + -0.000251576 *begin_x*lens_ipow(begin_y, 2)*begin_dy + -0.000390339 *lens_ipow(begin_x, 3)*begin_y*lens_ipow(begin_dy, 2) + 0.00159912 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 3)*begin_dx*lens_ipow(begin_dy, 2) + 2.39644e-07 *lens_ipow(begin_x, 5)*begin_y*lens_ipow(begin_lambda, 2) + 1.0824e-08 *lens_ipow(begin_x, 7)*begin_dy + -6.52784e-09 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 6)*begin_dy + 0.00148592 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 3)+0.0f;
  dx1_domega0[1][1] =  + 33.028  + 85.8108 *lens_ipow(begin_dy, 2) + 28.8268 *lens_ipow(begin_dx, 2) + 7.52983 *begin_y*begin_dy + 0.121886 *lens_ipow(begin_y, 2) + 2.73214 *begin_x*begin_dx + 2.9143e-05 *begin_x*begin_y + 0.0445754 *lens_ipow(begin_x, 2) + -1.26068 *lens_ipow(begin_lambda, 4) + -0.000251576 *begin_x*lens_ipow(begin_y, 2)*begin_dx + 1.51681e-05 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 3)*begin_dy + -0.000780678 *lens_ipow(begin_x, 3)*begin_y*begin_dx*begin_dy + -6.75449e-05 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*lens_ipow(begin_lambda, 4) + 0.00159912 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 2)*begin_dy + 1.0824e-08 *lens_ipow(begin_x, 7)*begin_dx + 1.03377e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 6)*begin_lambda + 5.96074 *lens_ipow(begin_lambda, 10) + -6.52784e-09 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 6)*begin_dx + 0.00222887 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2)+0.0f;
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