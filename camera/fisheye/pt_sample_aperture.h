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
  pred_x =  + 0.000155008  + 44.0953 *begin_dx + 6.81398e-06 *begin_y + 0.656767 *begin_x + -0.00611466 *lens_ipow(begin_dy, 2) + -3.59568e-06 *begin_x*begin_y + 18.506 *begin_dx*lens_ipow(begin_dy, 2) + 18.7993 *lens_ipow(begin_dx, 3) + 1.45232 *begin_y*begin_dx*begin_dy + 0.0173492 *lens_ipow(begin_y, 2)*begin_dx + 0.274624 *begin_x*lens_ipow(begin_dy, 2) + 1.71924 *begin_x*lens_ipow(begin_dx, 2) + 0.0214982 *begin_x*begin_y*begin_dy + 0.000216627 *begin_x*lens_ipow(begin_y, 2) + 0.0384358 *lens_ipow(begin_x, 2)*begin_dx + 0.000207915 *lens_ipow(begin_x, 3) + -3.78161 *begin_dx*lens_ipow(begin_lambda, 3) + 4.90937 *lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 2) + -0.157309 *begin_x*lens_ipow(begin_lambda, 4) + 8.60202e-06 *begin_x*lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 3) + -3.30109e-06 *begin_x*lens_ipow(begin_y, 4)*lens_ipow(begin_dx, 2) + 0.0155756 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dx, 3)*begin_dy + 2472.91 *lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 6) + -0.00688976 *begin_x*lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 3) + -0.00124823 *lens_ipow(begin_x, 4)*lens_ipow(begin_dx, 5) + -4.61264e-08 *lens_ipow(begin_x, 5)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2) + 7.16037e-11 *lens_ipow(begin_x, 5)*lens_ipow(begin_y, 4) + -6.70113e-08 *lens_ipow(begin_x, 6)*begin_y*begin_dx*begin_dy + 1.43112e-09 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 4)*lens_ipow(begin_dy, 2) + 20.6195 *begin_dx*lens_ipow(begin_lambda, 10) + 26507.4 *lens_ipow(begin_dx, 7)*lens_ipow(begin_dy, 4) + 12409.3 *lens_ipow(begin_dx, 9)*lens_ipow(begin_dy, 2) + 1837.33 *lens_ipow(begin_dx, 11) + 0.770513 *begin_x*lens_ipow(begin_lambda, 10) + -0.0523361 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 5)*lens_ipow(begin_dy, 2) + -0.0148458 *lens_ipow(begin_x, 4)*lens_ipow(begin_dx, 5)*lens_ipow(begin_dy, 2) + -3.1209e-09 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 5)*begin_dx*begin_dy + 0.000196274 *lens_ipow(begin_x, 5)*lens_ipow(begin_dy, 6) + -2.18867e-10 *lens_ipow(begin_x, 9)*lens_ipow(begin_dx, 2) + 6.17016e-14 *lens_ipow(begin_x, 11);
  pred_y =  + -1.81429e-06  + 44.0985 *begin_dy + 0.656259 *begin_y + -0.00037931 *begin_x*begin_dy + 4.1033e-06 *begin_x*begin_y + 18.9055 *lens_ipow(begin_dy, 3) + 0.0455439 *begin_dx*lens_ipow(begin_dy, 2) + 18.414 *lens_ipow(begin_dx, 2)*begin_dy + 1.73807 *begin_y*lens_ipow(begin_dy, 2) + 0.262874 *begin_y*lens_ipow(begin_dx, 2) + 0.0390663 *lens_ipow(begin_y, 2)*begin_dy + -3.20286e-05 *lens_ipow(begin_y, 2)*begin_dx + 0.000229581 *lens_ipow(begin_y, 3) + 1.44886 *begin_x*begin_dx*begin_dy + 0.0211923 *begin_x*begin_y*begin_dx + 0.0173283 *lens_ipow(begin_x, 2)*begin_dy + 0.000224903 *lens_ipow(begin_x, 2)*begin_y + -3.8073 *begin_dy*lens_ipow(begin_lambda, 3) + 0.00848631 *begin_x*lens_ipow(begin_dy, 3) + -0.160619 *begin_y*lens_ipow(begin_lambda, 4) + 0.00237679 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + 0.00652786 *begin_x*lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 3) + 0.000642894 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dx, 4) + 4.72613e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 4)*begin_dy + 2572.86 *lens_ipow(begin_dx, 6)*lens_ipow(begin_dy, 3) + 0.0124217 *lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 6) + 2.60945e-05 *lens_ipow(begin_y, 5)*lens_ipow(begin_dy, 4) + 4.43839e-10 *lens_ipow(begin_y, 8)*begin_dy + -0.0045496 *lens_ipow(begin_x, 3)*begin_y*lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 2) + 1.88779e-09 *lens_ipow(begin_x, 5)*lens_ipow(begin_y, 3)*begin_dx + 7.33568e-10 *lens_ipow(begin_x, 6)*lens_ipow(begin_y, 2)*begin_dy + 20.9501 *begin_dy*lens_ipow(begin_lambda, 10) + 1167.6 *lens_ipow(begin_dy, 11) + 16011.4 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 9) + 0.789592 *begin_y*lens_ipow(begin_lambda, 10) + 2.20174 *begin_y*lens_ipow(begin_dy, 6)*lens_ipow(begin_lambda, 4) + 1.2535e-10 *lens_ipow(begin_y, 9)*lens_ipow(begin_dx, 2) + 1.04095 *begin_x*lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 7) + 2.56328e-08 *begin_x*lens_ipow(begin_y, 7)*begin_dx*lens_ipow(begin_dy, 2) + 6.16022 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 9);
  pred_dx =  + -5.38199e-06  + 0.0547469 *begin_dx + -0.0220138 *begin_x + 0.000302318 *lens_ipow(begin_dx, 2) + 2.97585e-05 *begin_y*begin_dx + -3.48758e-07 *begin_x*begin_y + 1.69644 *begin_dx*lens_ipow(begin_dy, 2) + 1.6834 *lens_ipow(begin_dx, 3) + 0.0909783 *begin_y*begin_dx*begin_dy + 0.000810681 *lens_ipow(begin_y, 2)*begin_dx + 0.0460194 *begin_x*lens_ipow(begin_dy, 2) + 0.136832 *begin_x*lens_ipow(begin_dx, 2) + 0.00211819 *begin_x*begin_y*begin_dy + 1.57983e-05 *begin_x*lens_ipow(begin_y, 2) + 0.00294841 *lens_ipow(begin_x, 2)*begin_dx + 1.58523e-05 *lens_ipow(begin_x, 3) + -0.260206 *begin_dx*lens_ipow(begin_lambda, 3) + 0.34215 *begin_dx*lens_ipow(begin_dy, 4) + 1.2156 *lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 2) + 0.400805 *lens_ipow(begin_dx, 5) + 0.000670459 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 3) + -1.50629e-07 *lens_ipow(begin_y, 4)*begin_dx + -0.00913681 *begin_x*lens_ipow(begin_lambda, 4) + -0.0136069 *begin_y*lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 2) + 0.00163636 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dx, 3)*begin_dy + 8.86087e-09 *lens_ipow(begin_x, 5)*begin_y*begin_dy + 1.27346e-05 *begin_x*lens_ipow(begin_y, 4)*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + -5.77914e-13 *begin_x*lens_ipow(begin_y, 8) + -9.83637e-05 *lens_ipow(begin_x, 4)*lens_ipow(begin_dx, 5) + 1.43868 *begin_dx*lens_ipow(begin_lambda, 10) + 1662.33 *lens_ipow(begin_dx, 5)*lens_ipow(begin_dy, 6) + 570.521 *lens_ipow(begin_dx, 9)*lens_ipow(begin_dy, 2) + 132.969 *lens_ipow(begin_dx, 11) + 0.0450834 *begin_x*lens_ipow(begin_lambda, 10) + -1.31685e-10 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 5)*begin_dx*begin_dy + 1.36862e-05 *lens_ipow(begin_x, 5)*lens_ipow(begin_dy, 6) + -1.69581e-09 *lens_ipow(begin_x, 7)*lens_ipow(begin_dx, 3)*begin_dy + 1.38263e-09 *lens_ipow(begin_x, 7)*begin_y*lens_ipow(begin_dx, 2)*begin_dy + 1.12276e-11 *lens_ipow(begin_x, 9)*lens_ipow(begin_dy, 2) + -1.13654e-11 *lens_ipow(begin_x, 9)*lens_ipow(begin_dx, 2);
  pred_dy =  + 5.13233e-06  + 0.0556752 *begin_dy + -9.12943e-06 *begin_dx + -0.0219592 *begin_y + 8.0622e-07 *begin_x*begin_y + 1.70502 *lens_ipow(begin_dy, 3) + 1.72927 *lens_ipow(begin_dx, 2)*begin_dy + 0.136777 *begin_y*lens_ipow(begin_dy, 2) + 0.0453447 *begin_y*lens_ipow(begin_dx, 2) + 0.00291988 *lens_ipow(begin_y, 2)*begin_dy + 1.54337e-05 *lens_ipow(begin_y, 3) + 0.0917248 *begin_x*begin_dx*begin_dy + 0.00209898 *begin_x*begin_y*begin_dx + 0.000794969 *lens_ipow(begin_x, 2)*begin_dy + 1.49161e-05 *lens_ipow(begin_x, 2)*begin_y + -0.265342 *begin_dy*lens_ipow(begin_lambda, 3) + 0.000421919 *begin_x*lens_ipow(begin_dy, 3) + 0.333965 *lens_ipow(begin_dy, 5) + 0.0407783 *begin_dx*lens_ipow(begin_dy, 4) + 0.793703 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 3) + 0.0373036 *lens_ipow(begin_dx, 4)*begin_dy + -0.00947545 *begin_y*lens_ipow(begin_lambda, 4) + 0.000363611 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 3) + 8.39396e-08 *lens_ipow(begin_x, 4)*begin_dy + 0.000160114 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + 0.0013312 *begin_x*lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 3) + -9.58619e-05 *lens_ipow(begin_y, 4)*lens_ipow(begin_dy, 5) + -8.39532e-09 *begin_x*lens_ipow(begin_y, 6)*begin_dx*begin_dy + 1.1055e-05 *lens_ipow(begin_x, 4)*begin_y*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + 1.46349 *begin_dy*lens_ipow(begin_lambda, 10) + 117.389 *lens_ipow(begin_dy, 11) + 841.701 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 9) + 3122.88 *lens_ipow(begin_dx, 6)*lens_ipow(begin_dy, 5) + 0.0472013 *begin_y*lens_ipow(begin_lambda, 10) + 1.51103e-05 *lens_ipow(begin_y, 5)*lens_ipow(begin_dx, 6) + -1.21439e-11 *lens_ipow(begin_y, 9)*lens_ipow(begin_dy, 2) + 1.28916e-11 *lens_ipow(begin_y, 9)*lens_ipow(begin_dx, 2) + 0.137746 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 5) + 0.127837 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 9) + -1.75791e-10 *lens_ipow(begin_x, 5)*lens_ipow(begin_y, 4)*begin_dx*begin_dy;
  float dx1_domega0[2][2];
  dx1_domega0[0][0] =  + 44.0953  + 18.506 *lens_ipow(begin_dy, 2) + 56.398 *lens_ipow(begin_dx, 2) + 1.45232 *begin_y*begin_dy + 0.0173492 *lens_ipow(begin_y, 2) + 3.43848 *begin_x*begin_dx + 0.0384358 *lens_ipow(begin_x, 2) + -3.78161 *lens_ipow(begin_lambda, 3) + 14.7281 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + -6.60218e-06 *begin_x*lens_ipow(begin_y, 4)*begin_dx + 0.0467267 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dx, 2)*begin_dy + 7418.74 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 6) + -0.0137795 *begin_x*lens_ipow(begin_y, 3)*begin_dx*lens_ipow(begin_dy, 3) + -0.00624113 *lens_ipow(begin_x, 4)*lens_ipow(begin_dx, 4) + -9.22528e-08 *lens_ipow(begin_x, 5)*lens_ipow(begin_y, 2)*begin_dx + -6.70113e-08 *lens_ipow(begin_x, 6)*begin_y*begin_dy + 20.6195 *lens_ipow(begin_lambda, 10) + 185552 *lens_ipow(begin_dx, 6)*lens_ipow(begin_dy, 4) + 111683 *lens_ipow(begin_dx, 8)*lens_ipow(begin_dy, 2) + 20210.6 *lens_ipow(begin_dx, 10) + -0.261681 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 2) + -0.0742288 *lens_ipow(begin_x, 4)*lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 2) + -3.1209e-09 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 5)*begin_dy + -4.37734e-10 *lens_ipow(begin_x, 9)*begin_dx+0.0f;
  dx1_domega0[0][1] =  + -0.0122293 *begin_dy + 37.012 *begin_dx*begin_dy + 1.45232 *begin_y*begin_dx + 0.549248 *begin_x*begin_dy + 0.0214982 *begin_x*begin_y + 9.81873 *lens_ipow(begin_dx, 3)*begin_dy + 2.58061e-05 *begin_x*lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 2) + 0.0155756 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dx, 3) + 14837.5 *lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 5) + -0.0206693 *begin_x*lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + -6.70113e-08 *lens_ipow(begin_x, 6)*begin_y*begin_dx + 2.86224e-09 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 4)*begin_dy + 106029 *lens_ipow(begin_dx, 7)*lens_ipow(begin_dy, 3) + 24818.5 *lens_ipow(begin_dx, 9)*begin_dy + -0.104672 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 5)*begin_dy + -0.0296915 *lens_ipow(begin_x, 4)*lens_ipow(begin_dx, 5)*begin_dy + -3.1209e-09 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 5)*begin_dx + 0.00117765 *lens_ipow(begin_x, 5)*lens_ipow(begin_dy, 5)+0.0f;
  dx1_domega0[1][0] =  + 0.0455439 *lens_ipow(begin_dy, 2) + 36.828 *begin_dx*begin_dy + 0.525747 *begin_y*begin_dx + -3.20286e-05 *lens_ipow(begin_y, 2) + 1.44886 *begin_x*begin_dy + 0.0211923 *begin_x*begin_y + 0.00475358 *lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 2) + 0.00652786 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 3) + 0.00257158 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dx, 3) + 15437.2 *lens_ipow(begin_dx, 5)*lens_ipow(begin_dy, 3) + 0.0745304 *lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 5) + -0.0136488 *lens_ipow(begin_x, 3)*begin_y*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + 1.88779e-09 *lens_ipow(begin_x, 5)*lens_ipow(begin_y, 3) + 32022.8 *begin_dx*lens_ipow(begin_dy, 9) + 2.50699e-10 *lens_ipow(begin_y, 9)*begin_dx + 1.04095 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 7) + 2.56328e-08 *begin_x*lens_ipow(begin_y, 7)*lens_ipow(begin_dy, 2)+0.0f;
  dx1_domega0[1][1] =  + 44.0985  + -0.00037931 *begin_x + 56.7165 *lens_ipow(begin_dy, 2) + 0.0910878 *begin_dx*begin_dy + 18.414 *lens_ipow(begin_dx, 2) + 3.47614 *begin_y*begin_dy + 0.0390663 *lens_ipow(begin_y, 2) + 1.44886 *begin_x*begin_dx + 0.0173283 *lens_ipow(begin_x, 2) + -3.8073 *lens_ipow(begin_lambda, 3) + 0.0254589 *begin_x*lens_ipow(begin_dy, 2) + 0.00475358 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*begin_dy + 0.0195836 *begin_x*lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 2) + 4.72613e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 4) + 7718.59 *lens_ipow(begin_dx, 6)*lens_ipow(begin_dy, 2) + 0.000104378 *lens_ipow(begin_y, 5)*lens_ipow(begin_dy, 3) + 4.43839e-10 *lens_ipow(begin_y, 8) + -0.0090992 *lens_ipow(begin_x, 3)*begin_y*lens_ipow(begin_dx, 3)*begin_dy + 7.33568e-10 *lens_ipow(begin_x, 6)*lens_ipow(begin_y, 2) + 20.9501 *lens_ipow(begin_lambda, 10) + 12843.6 *lens_ipow(begin_dy, 10) + 144103 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 8) + 13.2104 *begin_y*lens_ipow(begin_dy, 5)*lens_ipow(begin_lambda, 4) + 7.28668 *begin_x*lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 6) + 5.12656e-08 *begin_x*lens_ipow(begin_y, 7)*begin_dx*begin_dy + 55.442 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 8)+0.0f;
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