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
  pred_x =  + -8.12453e-05  + 0.000911091 *begin_dy + 198.216 *begin_dx + 0.829483 *begin_x + 0.000260821 *begin_y*begin_dx + 199.292 *begin_dx*lens_ipow(begin_dy, 2) + 202.748 *lens_ipow(begin_dx, 3) + 2.35854 *begin_y*begin_dx*begin_dy + 0.00602257 *lens_ipow(begin_y, 2)*begin_dx + 0.829079 *begin_x*lens_ipow(begin_dy, 2) + 3.18621 *begin_x*lens_ipow(begin_dx, 2) + 0.00934177 *begin_x*begin_y*begin_dy + 2.44959e-05 *begin_x*lens_ipow(begin_y, 2) + 0.0159097 *lens_ipow(begin_x, 2)*begin_dx + 1.24937e-07 *lens_ipow(begin_x, 2)*begin_y + 2.33169e-05 *lens_ipow(begin_x, 3) + 9.01885 *begin_dx*lens_ipow(begin_lambda, 3) + -0.0269557 *lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 2) + 0.0445591 *begin_x*lens_ipow(begin_lambda, 4) + -0.0429779 *begin_x*begin_y*lens_ipow(begin_dx, 2)*begin_dy + -0.000174347 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 2) + -0.0206592 *lens_ipow(begin_x, 2)*begin_dx*lens_ipow(begin_dy, 2) + -0.0462175 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 3) + -0.000464986 *lens_ipow(begin_x, 2)*begin_y*begin_dx*begin_dy + -0.211982 *begin_y*lens_ipow(begin_dx, 5) + -0.000230059 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*begin_lambda + 2175.58 *begin_dx*lens_ipow(begin_dy, 6) + 4983.1 *lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 4) + 6209.88 *lens_ipow(begin_dx, 5)*lens_ipow(begin_dy, 2) + -0.339065 *begin_x*lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 4) + -0.160331 *begin_x*begin_y*lens_ipow(begin_dy, 5) + -0.679798 *begin_x*begin_y*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 3) + -6.68873e-07 *lens_ipow(begin_x, 5)*lens_ipow(begin_dx, 2) + -37.6688 *begin_dx*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 5) + -47.2747 *lens_ipow(begin_dx, 3)*lens_ipow(begin_lambda, 5) + -0.000680194 *lens_ipow(begin_y, 3)*begin_dx*begin_dy*lens_ipow(begin_lambda, 3) + 18239.3 *lens_ipow(begin_dx, 9) + -41.5475 *begin_dx*lens_ipow(begin_lambda, 10) + -0.20553 *begin_x*lens_ipow(begin_lambda, 10) + 1.19791 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dx, 7)*begin_dy;
  pred_y =  + 9.81663e-05  + 198.201 *begin_dy + 0.82946 *begin_y + 8.079e-06 *begin_x + -0.000117274 *begin_x*begin_dy + 203.019 *lens_ipow(begin_dy, 3) + 199.248 *lens_ipow(begin_dx, 2)*begin_dy + 3.19676 *begin_y*lens_ipow(begin_dy, 2) + 0.829491 *begin_y*lens_ipow(begin_dx, 2) + 0.0159878 *lens_ipow(begin_y, 2)*begin_dy + 2.31989e-05 *lens_ipow(begin_y, 3) + 2.33801 *begin_x*begin_dx*begin_dy + 0.00935804 *begin_x*begin_y*begin_dx + 0.00607021 *lens_ipow(begin_x, 2)*begin_dy + 2.44656e-05 *lens_ipow(begin_x, 2)*begin_y + 9.06814 *begin_dy*lens_ipow(begin_lambda, 3) + 0.0451561 *begin_y*lens_ipow(begin_lambda, 4) + -0.0457117 *lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 3) + -0.0182267 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*begin_dy + -0.043742 *begin_x*begin_y*begin_dx*lens_ipow(begin_dy, 2) + -0.0260825 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 2)*begin_dy + -0.000144031 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dy, 2) + -0.000170775 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dx, 2) + 6101.6 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 5) + 4852.31 *lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 3) + 2156.32 *lens_ipow(begin_dx, 6)*begin_dy + -0.408 *begin_y*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 4) + -7.055e-07 *lens_ipow(begin_y, 5)*lens_ipow(begin_dy, 2) + -0.631091 *begin_x*begin_y*lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 2) + -0.163327 *begin_x*begin_y*lens_ipow(begin_dx, 5) + -0.000953331 *begin_x*lens_ipow(begin_y, 2)*begin_dx*begin_dy*lens_ipow(begin_lambda, 2) + -6.39755e-07 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 2)*begin_dx*begin_dy + -52.4451 *lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 5) + -40.5936 *lens_ipow(begin_dx, 2)*begin_dy*lens_ipow(begin_lambda, 5) + 18162.3 *lens_ipow(begin_dy, 9) + -41.2104 *begin_dy*lens_ipow(begin_lambda, 10) + -0.2093 *begin_y*lens_ipow(begin_lambda, 10) + 1.06096 *begin_x*lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 7) + -3.94619e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 4)*begin_dy*lens_ipow(begin_lambda, 4) + -0.00106098 *lens_ipow(begin_x, 3)*begin_dx*begin_dy*lens_ipow(begin_lambda, 6);
  pred_dx =  + -4.49536e-06  + -3.25995e-05 *begin_dy + -0.0517651 *begin_dx + -6.16509e-07 *begin_y + -0.0052278 *begin_x + -0.460396 *begin_dx*lens_ipow(begin_dy, 2) + -0.452101 *lens_ipow(begin_dx, 3) + 0.00938814 *begin_y*begin_dx*begin_dy + 4.14432e-07 *lens_ipow(begin_y, 2)*begin_dy + 3.60217e-05 *lens_ipow(begin_y, 2)*begin_dx + 0.00505714 *begin_x*lens_ipow(begin_dy, 2) + 0.0143396 *begin_x*lens_ipow(begin_dx, 2) + 9.39294e-05 *begin_x*begin_y*begin_dy + 5.91683e-07 *begin_x*begin_y*begin_dx + 2.2259e-07 *begin_x*lens_ipow(begin_y, 2) + 7.2069e-09 *lens_ipow(begin_x, 2)*begin_y + 2.36428e-07 *lens_ipow(begin_x, 3) + 0.281101 *begin_dx*lens_ipow(begin_lambda, 3) + 0.00027694 *lens_ipow(begin_x, 2)*begin_dx*begin_lambda + -9.91916e-05 *lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 2) + 0.00248041 *begin_x*lens_ipow(begin_lambda, 4) + -0.000494963 *begin_x*begin_y*lens_ipow(begin_dx, 2)*begin_dy + -4.38711e-06 *lens_ipow(begin_x, 2)*begin_y*begin_dx*begin_dy + -0.00174928 *begin_y*lens_ipow(begin_dx, 5) + 10.5182 *begin_dx*lens_ipow(begin_dy, 6) + 20.5008 *lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 4) + 31.3906 *lens_ipow(begin_dx, 5)*lens_ipow(begin_dy, 2) + 4.94399e-06 *lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 4) + -0.00120294 *begin_x*begin_y*lens_ipow(begin_dy, 5) + -1.12464e-05 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 4) + -0.000861726 *lens_ipow(begin_x, 2)*begin_dx*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 2) + -0.00118595 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 3)*lens_ipow(begin_lambda, 2) + -6.09736e-09 *lens_ipow(begin_x, 5)*lens_ipow(begin_dx, 2) + -0.0290002 *begin_y*lens_ipow(begin_dx, 3)*begin_dy*lens_ipow(begin_lambda, 3) + 89.7027 *lens_ipow(begin_dx, 9) + -0.00356106 *begin_x*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 6) + -1.38116 *begin_dx*lens_ipow(begin_lambda, 10) + -0.0119907 *begin_x*lens_ipow(begin_lambda, 10) + -0.0212661 *begin_x*lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 8) + -0.00120855 *lens_ipow(begin_x, 2)*begin_dx*lens_ipow(begin_lambda, 8);
  pred_dy =  + -1.14516e-06  + -0.0518182 *begin_dy + -0.00522682 *begin_y + -3.74672e-07 *begin_x + 4.25837e-08 *lens_ipow(begin_y, 2) + -4.72128e-06 *begin_x*begin_dy + 2.38432e-06 *begin_x*begin_dx + -0.451699 *lens_ipow(begin_dy, 3) + -0.462198 *lens_ipow(begin_dx, 2)*begin_dy + 0.0143555 *begin_y*lens_ipow(begin_dy, 2) + -2.86649e-05 *begin_y*begin_dx*begin_dy + 0.00495887 *begin_y*lens_ipow(begin_dx, 2) + 2.3188e-07 *lens_ipow(begin_y, 3) + 4.82065e-05 *begin_x*lens_ipow(begin_dy, 2) + 0.0092764 *begin_x*begin_dx*begin_dy + 4.55767e-07 *begin_x*begin_y*begin_dy + 9.33646e-05 *begin_x*begin_y*begin_dx + 3.78189e-05 *lens_ipow(begin_x, 2)*begin_dy + 2.22293e-07 *lens_ipow(begin_x, 2)*begin_y + 0.280513 *begin_dy*lens_ipow(begin_lambda, 3) + 0.000154186 *begin_y*begin_dx*lens_ipow(begin_dy, 2) + 0.000279613 *lens_ipow(begin_y, 2)*begin_dy*begin_lambda + 0.00249751 *begin_y*lens_ipow(begin_lambda, 4) + -0.000451272 *begin_x*begin_y*begin_dx*lens_ipow(begin_dy, 2) + -4.47914e-06 *begin_x*lens_ipow(begin_y, 2)*begin_dx*begin_dy + -0.000136698 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 2)*begin_dy + 35.7671 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 5) + 11.2207 *lens_ipow(begin_dx, 6)*begin_dy + -0.00126652 *lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 2) + -0.000897379 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*begin_dy*lens_ipow(begin_lambda, 2) + -6.96379e-09 *lens_ipow(begin_y, 5)*lens_ipow(begin_dy, 2) + -0.00168473 *begin_x*begin_y*lens_ipow(begin_dx, 5) + -5.0409e-06 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 3) + 101.422 *lens_ipow(begin_dy, 9) + 69.6049 *lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 2) + -3.00568e-10 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 4)*begin_dy*lens_ipow(begin_lambda, 3) + -1.38057 *begin_dy*lens_ipow(begin_lambda, 10) + -0.0123033 *begin_y*lens_ipow(begin_lambda, 10) + -0.0218027 *begin_y*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 8) + -0.00120058 *lens_ipow(begin_y, 2)*begin_dy*lens_ipow(begin_lambda, 8);
  float dx1_domega0[2][2];
  dx1_domega0[0][0] =  + 198.216  + 0.000260821 *begin_y + 199.292 *lens_ipow(begin_dy, 2) + 608.243 *lens_ipow(begin_dx, 2) + 2.35854 *begin_y*begin_dy + 0.00602257 *lens_ipow(begin_y, 2) + 6.37243 *begin_x*begin_dx + 0.0159097 *lens_ipow(begin_x, 2) + 9.01885 *lens_ipow(begin_lambda, 3) + -0.0269557 *lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 2) + -0.0859558 *begin_x*begin_y*begin_dx*begin_dy + -0.0206592 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 2) + -0.138653 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 2) + -0.000464986 *lens_ipow(begin_x, 2)*begin_y*begin_dy + -1.05991 *begin_y*lens_ipow(begin_dx, 4) + -0.000460118 *begin_x*lens_ipow(begin_y, 2)*begin_dx*begin_lambda + 2175.58 *lens_ipow(begin_dy, 6) + 14949.3 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 4) + 31049.4 *lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 2) + -0.67813 *begin_x*begin_dx*lens_ipow(begin_lambda, 4) + -1.3596 *begin_x*begin_y*begin_dx*lens_ipow(begin_dy, 3) + -1.33775e-06 *lens_ipow(begin_x, 5)*begin_dx + -37.6688 *lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 5) + -141.824 *lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 5) + -0.000680194 *lens_ipow(begin_y, 3)*begin_dy*lens_ipow(begin_lambda, 3) + 164154 *lens_ipow(begin_dx, 8) + -41.5475 *lens_ipow(begin_lambda, 10) + 8.38536 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dx, 6)*begin_dy+0.0f;
  dx1_domega0[0][1] =  + 0.000911091  + 398.584 *begin_dx*begin_dy + 2.35854 *begin_y*begin_dx + 1.65816 *begin_x*begin_dy + 0.00934177 *begin_x*begin_y + -0.0539115 *lens_ipow(begin_y, 2)*begin_dx*begin_dy + -0.0429779 *begin_x*begin_y*lens_ipow(begin_dx, 2) + -0.000348694 *begin_x*lens_ipow(begin_y, 2)*begin_dy + -0.0413184 *lens_ipow(begin_x, 2)*begin_dx*begin_dy + -0.000464986 *lens_ipow(begin_x, 2)*begin_y*begin_dx + 13053.5 *begin_dx*lens_ipow(begin_dy, 5) + 19932.4 *lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 3) + 12419.8 *lens_ipow(begin_dx, 5)*begin_dy + -0.801654 *begin_x*begin_y*lens_ipow(begin_dy, 4) + -2.0394 *begin_x*begin_y*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + -75.3376 *begin_dx*begin_dy*lens_ipow(begin_lambda, 5) + -0.000680194 *lens_ipow(begin_y, 3)*begin_dx*lens_ipow(begin_lambda, 3) + 1.19791 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dx, 7)+0.0f;
  dx1_domega0[1][0] =  + 398.496 *begin_dx*begin_dy + 1.65898 *begin_y*begin_dx + 2.33801 *begin_x*begin_dy + 0.00935804 *begin_x*begin_y + -0.0364535 *lens_ipow(begin_y, 2)*begin_dx*begin_dy + -0.043742 *begin_x*begin_y*lens_ipow(begin_dy, 2) + -0.0521651 *lens_ipow(begin_x, 2)*begin_dx*begin_dy + -0.000341549 *lens_ipow(begin_x, 2)*begin_y*begin_dx + 12203.2 *begin_dx*lens_ipow(begin_dy, 5) + 19409.2 *lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 3) + 12937.9 *lens_ipow(begin_dx, 5)*begin_dy + -1.89327 *begin_x*begin_y*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + -0.816636 *begin_x*begin_y*lens_ipow(begin_dx, 4) + -0.000953331 *begin_x*lens_ipow(begin_y, 2)*begin_dy*lens_ipow(begin_lambda, 2) + -6.39755e-07 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 2)*begin_dy + -81.1871 *begin_dx*begin_dy*lens_ipow(begin_lambda, 5) + 1.06096 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 7) + -0.00106098 *lens_ipow(begin_x, 3)*begin_dy*lens_ipow(begin_lambda, 6)+0.0f;
  dx1_domega0[1][1] =  + 198.201  + -0.000117274 *begin_x + 609.058 *lens_ipow(begin_dy, 2) + 199.248 *lens_ipow(begin_dx, 2) + 6.39351 *begin_y*begin_dy + 0.0159878 *lens_ipow(begin_y, 2) + 2.33801 *begin_x*begin_dx + 0.00607021 *lens_ipow(begin_x, 2) + 9.06814 *lens_ipow(begin_lambda, 3) + -0.137135 *lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 2) + -0.0182267 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2) + -0.087484 *begin_x*begin_y*begin_dx*begin_dy + -0.0260825 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 2) + -0.000288061 *lens_ipow(begin_x, 2)*begin_y*begin_dy + 30508 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 4) + 14556.9 *lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 2) + 2156.32 *lens_ipow(begin_dx, 6) + -0.816 *begin_y*begin_dy*lens_ipow(begin_lambda, 4) + -1.411e-06 *lens_ipow(begin_y, 5)*begin_dy + -1.26218 *begin_x*begin_y*lens_ipow(begin_dx, 3)*begin_dy + -0.000953331 *begin_x*lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_lambda, 2) + -6.39755e-07 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 2)*begin_dx + -157.335 *lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 5) + -40.5936 *lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 5) + 163461 *lens_ipow(begin_dy, 8) + -41.2104 *lens_ipow(begin_lambda, 10) + 7.42671 *begin_x*lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 6) + -3.94619e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 4)*lens_ipow(begin_lambda, 4) + -0.00106098 *lens_ipow(begin_x, 3)*begin_dx*lens_ipow(begin_lambda, 6)+0.0f;
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