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
  pred_x =  + -3.44136e-06  + -0.000225677 *begin_dy + 43.0804 *begin_dx + 0.790034 *begin_x + -0.00255173 *begin_dx*begin_dy + -0.00350182 *lens_ipow(begin_dx, 2) + -3.14678e-05 *begin_x*begin_dy + -6.54549e-05 *begin_x*begin_dx + -3.98629 *begin_dx*lens_ipow(begin_dy, 2) + -3.61141 *lens_ipow(begin_dx, 3) + 0.0963336 *begin_y*begin_dx*begin_dy + 0.000418481 *lens_ipow(begin_y, 2)*begin_dx + -0.304969 *begin_x*lens_ipow(begin_dy, 2) + -0.0264966 *begin_x*lens_ipow(begin_dx, 2) + -0.00942134 *begin_x*begin_y*begin_dy + -0.000135801 *begin_x*lens_ipow(begin_y, 2) + -0.00346055 *lens_ipow(begin_x, 2)*begin_dx + -0.000122032 *lens_ipow(begin_x, 3) + 1.26908 *begin_dx*lens_ipow(begin_lambda, 3) + 0.0110795 *begin_x*lens_ipow(begin_lambda, 3) + -6.06289e-09 *lens_ipow(begin_x, 3)*begin_y + -0.00609 *begin_x*begin_y*lens_ipow(begin_dx, 2)*begin_dy + -6.6452e-05 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2) + 0.000158005 *lens_ipow(begin_x, 3)*lens_ipow(begin_dx, 2) + -0.0624715 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 5) + -0.125383 *begin_x*begin_y*lens_ipow(begin_dy, 5) + 1.72986e-05 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_lambda, 4) + -3.53209e-11 *begin_x*lens_ipow(begin_y, 6) + -0.429206 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 5) + -8.62198e-05 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 2) + 2.20033e-06 *lens_ipow(begin_x, 4)*begin_y*begin_dx*begin_dy + 1.88496e-08 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 2)*begin_dx + -7.77368e-07 *lens_ipow(begin_x, 5)*lens_ipow(begin_dy, 2) + 0.00178277 *begin_x*begin_y*begin_dy*lens_ipow(begin_lambda, 6) + -1.05692e-12 *lens_ipow(begin_x, 5)*lens_ipow(begin_y, 4) + -5.88897e-13 *lens_ipow(begin_x, 9) + -1.96868 *begin_x*lens_ipow(begin_dx, 4)*lens_ipow(begin_lambda, 5) + 3.91758e-08 *lens_ipow(begin_x, 5)*lens_ipow(begin_lambda, 5) + -6.7011 *begin_dx*lens_ipow(begin_lambda, 10) + -0.0646457 *begin_x*lens_ipow(begin_lambda, 10);
  pred_y =  + 5.15165e-06  + 46.704 *begin_dy + 0.90776 *begin_y + 5.66776e-07 *begin_x + -0.00683992 *begin_dx*begin_dy + -4.22064e-05 *begin_y*begin_dx + 3.70906e-05 *begin_x*begin_dy + -2.96909e-07 *lens_ipow(begin_x, 2) + 0.742631 *lens_ipow(begin_dy, 3) + -0.363209 *lens_ipow(begin_dx, 2)*begin_dy + 0.0154544 *begin_y*lens_ipow(begin_dy, 2) + -0.144491 *begin_y*lens_ipow(begin_dx, 2) + -0.00564172 *lens_ipow(begin_y, 2)*begin_dy + -0.000106534 *lens_ipow(begin_y, 3) + -0.00092955 *begin_x*lens_ipow(begin_dy, 2) + 0.242385 *begin_x*begin_dx*begin_dy + -1.478e-05 *begin_x*begin_y*begin_dy + -0.00636909 *begin_x*begin_y*begin_dx + 0.00356562 *lens_ipow(begin_x, 2)*begin_dy + -0.000105678 *lens_ipow(begin_x, 2)*begin_y + 0.0175459 *begin_y*lens_ipow(begin_lambda, 3) + 2.37133 *begin_dy*lens_ipow(begin_lambda, 4) + -0.0116404 *begin_x*begin_y*begin_dx*lens_ipow(begin_dy, 2) + -0.00376833 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 2)*begin_dy + -0.0900307 *lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 5) + 3.28891e-05 *begin_x*lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 3) + 7.85147e-07 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 2)*begin_dx*begin_dy + 1.59512e-05 *lens_ipow(begin_x, 4)*lens_ipow(begin_dy, 3) + -1.78105e-10 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 3) + -1.43124e-11 *lens_ipow(begin_x, 6)*begin_y + -0.0773146 *begin_y*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 5) + 0.00131756 *begin_x*begin_y*begin_dx*lens_ipow(begin_lambda, 5) + 2.4764e-08 *lens_ipow(begin_y, 5)*lens_ipow(begin_lambda, 4) + -7.76224e-14 *lens_ipow(begin_y, 9) + -3.62324e-13 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 7) + 1.02381e-11 *lens_ipow(begin_x, 8)*begin_dy + 1.43774e-07 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 3)*lens_ipow(begin_lambda, 5) + 2.49957e-06 *lens_ipow(begin_x, 4)*begin_y*lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 3) + -11.4874 *begin_dy*lens_ipow(begin_lambda, 10) + -0.0958202 *begin_y*lens_ipow(begin_lambda, 10);
  pred_dx =  + 3.42879e-08  + -9.73018e-06 *begin_dy + 0.0820856 *begin_dx + -0.0217531 *begin_x + -0.000760571 *begin_dx*begin_dy + -0.000563343 *lens_ipow(begin_dx, 2) + -6.38104e-06 *begin_x*begin_dy + -1.1525e-05 *begin_x*begin_dx + -0.288547 *begin_dx*lens_ipow(begin_dy, 2) + -0.328124 *lens_ipow(begin_dx, 3) + 0.01282 *begin_y*begin_dx*begin_dy + 6.49129e-05 *lens_ipow(begin_y, 2)*begin_dx + 0.00316486 *begin_x*lens_ipow(begin_dy, 2) + 0.00996627 *begin_x*lens_ipow(begin_dx, 2) + 0.000370356 *begin_x*begin_y*begin_dy + -2.57057e-06 *begin_x*lens_ipow(begin_y, 2) + 0.000246851 *lens_ipow(begin_x, 2)*begin_dx + -5.29607e-06 *lens_ipow(begin_x, 3) + 0.00309456 *begin_x*lens_ipow(begin_lambda, 3) + -1.20209e-09 *lens_ipow(begin_x, 3)*begin_y + 0.210206 *begin_dx*lens_ipow(begin_lambda, 4) + 1.09497e-06 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_lambda, 2) + 1.17936e-05 *lens_ipow(begin_x, 2)*begin_y*begin_dx*begin_dy + 1.15206e-06 *lens_ipow(begin_x, 3)*lens_ipow(begin_lambda, 2) + 8.28538e-06 *lens_ipow(begin_x, 3)*lens_ipow(begin_dx, 2) + -4.98327e-06 *begin_x*begin_y*begin_dy*lens_ipow(begin_lambda, 3) + -0.0109514 *begin_x*begin_y*lens_ipow(begin_dy, 5) + 0.000276551 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 4) + -1.43042e-11 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 4) + -2.87827e-06 *lens_ipow(begin_x, 4)*begin_dx*lens_ipow(begin_dy, 2) + -1.46751e-10 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 2)*begin_dy + 1.1459e-09 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 2)*begin_dx + -5.29478e-12 *lens_ipow(begin_x, 7) + -2.47003e-07 *lens_ipow(begin_x, 4)*begin_y*begin_dx*lens_ipow(begin_dy, 2) + 1.17968e-12 *begin_x*lens_ipow(begin_y, 7)*begin_dy + 7.41204e-12 *lens_ipow(begin_x, 5)*lens_ipow(begin_y, 3)*begin_dy + 3.42937e-12 *lens_ipow(begin_x, 7)*begin_y*begin_dy + 1.7579e-12 *lens_ipow(begin_x, 8)*begin_dx + -1.01807 *begin_dx*lens_ipow(begin_lambda, 10) + -0.017547 *begin_x*lens_ipow(begin_lambda, 10);
  pred_dy =  + 3.20886e-06  + 0.222157 *begin_dy + -0.0171876 *begin_y + -4.5178e-07 *begin_x + -0.000631418 *begin_dx*begin_dy + 7.16646e-06 *begin_y*begin_dy + -5.52e-06 *begin_y*begin_dx + 5.97241e-06 *begin_x*begin_dy + -4.4707e-08 *lens_ipow(begin_x, 2) + -0.199668 *lens_ipow(begin_dy, 3) + -0.20003 *lens_ipow(begin_dx, 2)*begin_dy + 0.0190858 *begin_y*lens_ipow(begin_dy, 2) + -6.46025e-05 *begin_y*begin_dx*begin_dy + 0.0053591 *begin_y*lens_ipow(begin_dx, 2) + 0.000476776 *lens_ipow(begin_y, 2)*begin_dy + -1.39398e-06 *lens_ipow(begin_y, 3) + 0.00893331 *begin_x*begin_dx*begin_dy + 0.00027695 *begin_x*begin_y*begin_dx + 2.69198e-09 *begin_x*lens_ipow(begin_y, 2) + -4.55926e-06 *lens_ipow(begin_x, 2)*begin_y + 0.00369838 *begin_y*lens_ipow(begin_lambda, 3) + -0.000686672 *begin_y*lens_ipow(begin_dx, 2)*begin_dy + -1.22166e-05 *begin_x*begin_y*lens_ipow(begin_dx, 2) + 8.1788e-10 *lens_ipow(begin_x, 3)*begin_y + 0.237473 *begin_dy*lens_ipow(begin_lambda, 4) + -0.00019531 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 3) + -0.000321049 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 2)*begin_dy + -5.25444e-08 *lens_ipow(begin_x, 4)*begin_dy + -0.000820706 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*begin_dy*begin_lambda + -0.000977944 *begin_x*begin_y*lens_ipow(begin_dx, 3)*lens_ipow(begin_lambda, 2) + -1.42901e-08 *begin_x*lens_ipow(begin_y, 4)*begin_dx*begin_dy + -1.73111e-11 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 5) + 6.14492e-10 *lens_ipow(begin_x, 5)*begin_y*begin_dx + 1.98679e-09 *lens_ipow(begin_y, 5)*lens_ipow(begin_lambda, 3) + -9.62433e-15 *lens_ipow(begin_y, 9) + 1.14817e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 3)*lens_ipow(begin_lambda, 4) + -1.044e-13 *lens_ipow(begin_x, 6)*lens_ipow(begin_y, 3) + -1.15222 *begin_dy*lens_ipow(begin_lambda, 10) + -0.0202784 *begin_y*lens_ipow(begin_lambda, 10) + 2.29323e-12 *lens_ipow(begin_x, 5)*lens_ipow(begin_y, 4)*lens_ipow(begin_dx, 2);
  float dx1_domega0[2][2];
  dx1_domega0[0][0] =  + 43.0804  + -0.00255173 *begin_dy + -0.00700365 *begin_dx + -6.54549e-05 *begin_x + -3.98629 *lens_ipow(begin_dy, 2) + -10.8342 *lens_ipow(begin_dx, 2) + 0.0963336 *begin_y*begin_dy + 0.000418481 *lens_ipow(begin_y, 2) + -0.0529932 *begin_x*begin_dx + -0.00346055 *lens_ipow(begin_x, 2) + 1.26908 *lens_ipow(begin_lambda, 3) + -0.01218 *begin_x*begin_y*begin_dx*begin_dy + -0.000132904 *begin_x*lens_ipow(begin_y, 2)*begin_dx + 0.00031601 *lens_ipow(begin_x, 3)*begin_dx + -0.312357 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 4) + -2.14603 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 4) + -8.62198e-05 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 2) + 2.20033e-06 *lens_ipow(begin_x, 4)*begin_y*begin_dy + 1.88496e-08 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 2) + -7.87471 *begin_x*lens_ipow(begin_dx, 3)*lens_ipow(begin_lambda, 5) + -6.7011 *lens_ipow(begin_lambda, 10)+0.0f;
  dx1_domega0[0][1] =  + -0.000225677  + -0.00255173 *begin_dx + -3.14678e-05 *begin_x + -7.97259 *begin_dx*begin_dy + 0.0963336 *begin_y*begin_dx + -0.609938 *begin_x*begin_dy + -0.00942134 *begin_x*begin_y + -0.00609 *begin_x*begin_y*lens_ipow(begin_dx, 2) + -0.626917 *begin_x*begin_y*lens_ipow(begin_dy, 4) + -0.00017244 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*begin_dx*begin_dy + 2.20033e-06 *lens_ipow(begin_x, 4)*begin_y*begin_dx + -1.55474e-06 *lens_ipow(begin_x, 5)*begin_dy + 0.00178277 *begin_x*begin_y*lens_ipow(begin_lambda, 6)+0.0f;
  dx1_domega0[1][0] =  + -0.00683992 *begin_dy + -4.22064e-05 *begin_y + -0.726419 *begin_dx*begin_dy + -0.288983 *begin_y*begin_dx + 0.242385 *begin_x*begin_dy + -0.00636909 *begin_x*begin_y + -0.0116404 *begin_x*begin_y*lens_ipow(begin_dy, 2) + -0.00753667 *lens_ipow(begin_x, 2)*begin_dx*begin_dy + 9.86674e-05 *begin_x*lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 2) + 7.85147e-07 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 2)*begin_dy + 0.00131756 *begin_x*begin_y*lens_ipow(begin_lambda, 5) + 4.99913e-06 *lens_ipow(begin_x, 4)*begin_y*begin_dx*lens_ipow(begin_lambda, 3)+0.0f;
  dx1_domega0[1][1] =  + 46.704  + -0.00683992 *begin_dx + 3.70906e-05 *begin_x + 2.22789 *lens_ipow(begin_dy, 2) + -0.363209 *lens_ipow(begin_dx, 2) + 0.0309089 *begin_y*begin_dy + -0.00564172 *lens_ipow(begin_y, 2) + -0.0018591 *begin_x*begin_dy + 0.242385 *begin_x*begin_dx + -1.478e-05 *begin_x*begin_y + 0.00356562 *lens_ipow(begin_x, 2) + 2.37133 *lens_ipow(begin_lambda, 4) + -0.0232807 *begin_x*begin_y*begin_dx*begin_dy + -0.00376833 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 2) + -0.450153 *lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 4) + 7.85147e-07 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 2)*begin_dx + 4.78536e-05 *lens_ipow(begin_x, 4)*lens_ipow(begin_dy, 2) + -0.154629 *begin_y*begin_dy*lens_ipow(begin_lambda, 5) + 1.02381e-11 *lens_ipow(begin_x, 8) + -11.4874 *lens_ipow(begin_lambda, 10)+0.0f;
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