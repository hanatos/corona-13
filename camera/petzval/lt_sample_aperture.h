//input: scene_[x,y,z] - point in scene, ap_[x,y] - point on aperture
//output: [x,y,dx,dy] point and direction on sensor
#ifndef DEBUG_LOG
#define DEBUG_LOG
#endif
float view[3] =
{
  scene_x,
  scene_y,
  scene_z + lens_outer_pupil_curvature_radius
};
normalise(view);
int error = 0;
if(1 || view[2] >= lens_field_of_view)
{
  const float eps = 1e-8;
  float sqr_err = 1e30, sqr_ap_err = 1e30;
  float prev_sqr_err = 1e32, prev_sqr_ap_err = 1e32;
  for(int k=0;k<100&&(sqr_err>eps||sqr_ap_err>eps)&&error==0;k++)
  {
    prev_sqr_err = sqr_err, prev_sqr_ap_err = sqr_ap_err;
    const float begin_x = x;
    const float begin_y = y;
    const float begin_dx = dx;
    const float begin_dy = dy;
    const float begin_lambda = lambda;
    const float pred_ap[2] = {
       + -2.39596e-05  + 0.000170721 *begin_dy + 50.3715 *begin_dx + 1.20573e-06 *begin_y + 0.792054 *begin_x + 0.003814 *lens_ipow(begin_dx, 2) + -1.02565e-06 *lens_ipow(begin_x, 2) + 0.101002 *begin_dx*lens_ipow(begin_dy, 2) + 0.240926 *lens_ipow(begin_dx, 3) + 0.534569 *begin_y*begin_dx*begin_dy + 0.00916713 *lens_ipow(begin_y, 2)*begin_dx + 0.123821 *begin_x*lens_ipow(begin_dy, 2) + 0.658175 *begin_x*lens_ipow(begin_dx, 2) + 0.013701 *begin_x*begin_y*begin_dy + 0.000166393 *begin_x*lens_ipow(begin_y, 2) + -7.99133e-06 *lens_ipow(begin_x, 2)*begin_dy + 0.022954 *lens_ipow(begin_x, 2)*begin_dx + 0.000168299 *lens_ipow(begin_x, 3) + -0.0537504 *begin_dx*lens_ipow(begin_dy, 3) + 0.0538104 *lens_ipow(begin_dx, 3)*begin_dy + 0.0553584 *begin_x*lens_ipow(begin_lambda, 3) + 3.57819 *begin_dx*lens_ipow(begin_lambda, 4) + 0.00194404 *lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 2) + 6.06162e-06 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_lambda, 3) + 0.00031874 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 4) + 0.00750986 *begin_x*begin_y*lens_ipow(begin_dy, 5) + 0.0268014 *lens_ipow(begin_x, 2)*begin_dx*lens_ipow(begin_dy, 4) + 4.90502e-05 *lens_ipow(begin_x, 3)*begin_y*lens_ipow(begin_dx, 2)*begin_dy + 1.80118e-09 *lens_ipow(begin_x, 6)*begin_dx + 0.0583032 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 6)*begin_dy + -0.144778 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 5)*lens_ipow(begin_dy, 2) + 1.88903e-10 *lens_ipow(begin_x, 5)*lens_ipow(begin_y, 2)*lens_ipow(begin_lambda, 3) + 9.66942e-11 *lens_ipow(begin_x, 7)*lens_ipow(begin_lambda, 3) + -17.3519 *begin_dx*lens_ipow(begin_lambda, 10) + -2.25406 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 9) + -0.298715 *begin_x*lens_ipow(begin_lambda, 10) + 0.691107 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 9) + 1.13005e-15 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 8) + 2.36409e-08 *lens_ipow(begin_x, 7)*lens_ipow(begin_dy, 4) + 1.25173e-15 *lens_ipow(begin_x, 9)*lens_ipow(begin_y, 2),
       + 5.36935e-05  + 50.4606 *begin_dy + -7.98426e-05 *begin_dx + 0.793472 *begin_y + -7.02206e-07 *begin_x + 2.57101e-05 *begin_y*begin_dy + -6.73559e-06 *begin_x*begin_dy + -11.1018 *lens_ipow(begin_dy, 3) + -1.28628 *lens_ipow(begin_dx, 2)*begin_dy + 0.104314 *begin_y*lens_ipow(begin_dx, 2) + 0.0185838 *lens_ipow(begin_y, 2)*begin_dy + 0.00014928 *lens_ipow(begin_y, 3) + 0.504188 *begin_x*begin_dx*begin_dy + 0.013111 *begin_x*begin_y*begin_dx + 0.00900447 *lens_ipow(begin_x, 2)*begin_dy + 0.000161239 *lens_ipow(begin_x, 2)*begin_y + 0.0490546 *begin_y*lens_ipow(begin_lambda, 3) + 1.01948 *begin_y*lens_ipow(begin_dy, 2)*begin_lambda + 2.62835 *begin_dy*lens_ipow(begin_lambda, 4) + 35.2614 *lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 2) + -0.00169852 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 3) + 3.99942 *lens_ipow(begin_dx, 2)*begin_dy*lens_ipow(begin_lambda, 3) + 0.0127049 *lens_ipow(begin_y, 2)*begin_dy*lens_ipow(begin_lambda, 3) + 0.0430536 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 3) + 0.00229183 *begin_x*begin_y*lens_ipow(begin_dx, 5) + 4.37134e-05 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*begin_dy + 0.0655604 *begin_y*lens_ipow(begin_dx, 4)*lens_ipow(begin_lambda, 4) + -0.00308816 *lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 6) + -0.00228635 *lens_ipow(begin_x, 2)*begin_dy*lens_ipow(begin_lambda, 6) + -17.2515 *lens_ipow(begin_dy, 5)*lens_ipow(begin_lambda, 5) + 0.00147988 *lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 5) + -3.60928e-14 *lens_ipow(begin_y, 9)*begin_lambda + -9.44344e-14 *begin_x*lens_ipow(begin_y, 8)*begin_dx + -12.1588 *begin_dy*lens_ipow(begin_lambda, 10) + -77.059 *lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 8) + -0.24787 *begin_y*lens_ipow(begin_lambda, 10) + 2.87159e-07 *lens_ipow(begin_y, 5)*lens_ipow(begin_lambda, 6) + 0.000274665 *begin_x*lens_ipow(begin_y, 3)*begin_dx*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 4) + 1.22456e-12 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 7)*lens_ipow(begin_lambda, 2) + 4.54444e-15 *lens_ipow(begin_x, 6)*lens_ipow(begin_y, 5)
    };
    const float delta_ap[] = {ap_x - pred_ap[0], ap_y - pred_ap[1]};
    sqr_ap_err = delta_ap[0]*delta_ap[0]+delta_ap[1]*delta_ap[1];
    float dx1_domega0[2][2];
    dx1_domega0[0][0] =  + 50.3715  + 0.007628 *begin_dx + 0.101002 *lens_ipow(begin_dy, 2) + 0.722778 *lens_ipow(begin_dx, 2) + 0.534569 *begin_y*begin_dy + 0.00916713 *lens_ipow(begin_y, 2) + 1.31635 *begin_x*begin_dx + 0.022954 *lens_ipow(begin_x, 2) + -0.0537504 *lens_ipow(begin_dy, 3) + 0.161431 *lens_ipow(begin_dx, 2)*begin_dy + 3.57819 *lens_ipow(begin_lambda, 4) + 0.00194404 *lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 2) + 0.0268014 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 4) + 9.81004e-05 *lens_ipow(begin_x, 3)*begin_y*begin_dx*begin_dy + 1.80118e-09 *lens_ipow(begin_x, 6) + 0.349819 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 5)*begin_dy + -0.723889 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 2) + -17.3519 *lens_ipow(begin_lambda, 10) + -20.2866 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 8) + 6.21996 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 8)+0.0f;
    dx1_domega0[0][1] =  + 0.000170721  + 0.202004 *begin_dx*begin_dy + 0.534569 *begin_y*begin_dx + 0.247643 *begin_x*begin_dy + 0.013701 *begin_x*begin_y + -7.99133e-06 *lens_ipow(begin_x, 2) + -0.161251 *begin_dx*lens_ipow(begin_dy, 2) + 0.0538104 *lens_ipow(begin_dx, 3) + 0.00388809 *lens_ipow(begin_y, 2)*begin_dx*begin_dy + 0.00127496 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 3) + 0.0375493 *begin_x*begin_y*lens_ipow(begin_dy, 4) + 0.107206 *lens_ipow(begin_x, 2)*begin_dx*lens_ipow(begin_dy, 3) + 4.90502e-05 *lens_ipow(begin_x, 3)*begin_y*lens_ipow(begin_dx, 2) + 0.0583032 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 6) + -0.289556 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 5)*begin_dy + 9.45637e-08 *lens_ipow(begin_x, 7)*lens_ipow(begin_dy, 3)+0.0f;
    dx1_domega0[1][0] =  + -7.98426e-05  + -2.57256 *begin_dx*begin_dy + 0.208628 *begin_y*begin_dx + 0.504188 *begin_x*begin_dy + 0.013111 *begin_x*begin_y + 7.99884 *begin_dx*begin_dy*lens_ipow(begin_lambda, 3) + 0.0861072 *lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 3) + 0.0114591 *begin_x*begin_y*lens_ipow(begin_dx, 4) + 8.74268e-05 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*begin_dx*begin_dy + 0.262242 *begin_y*lens_ipow(begin_dx, 3)*lens_ipow(begin_lambda, 4) + 0.00295976 *lens_ipow(begin_y, 3)*begin_dx*lens_ipow(begin_lambda, 5) + -9.44344e-14 *begin_x*lens_ipow(begin_y, 8) + 0.000274665 *begin_x*lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 4)+0.0f;
    dx1_domega0[1][1] =  + 50.4606  + 2.57101e-05 *begin_y + -6.73559e-06 *begin_x + -33.3055 *lens_ipow(begin_dy, 2) + -1.28628 *lens_ipow(begin_dx, 2) + 0.0185838 *lens_ipow(begin_y, 2) + 0.504188 *begin_x*begin_dx + 0.00900447 *lens_ipow(begin_x, 2) + 2.03896 *begin_y*begin_dy*begin_lambda + 2.62835 *lens_ipow(begin_lambda, 4) + 105.784 *lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 2) + -0.00509555 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 2) + 3.99942 *lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 3) + 0.0127049 *lens_ipow(begin_y, 2)*lens_ipow(begin_lambda, 3) + 0.129161 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + 4.37134e-05 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2) + -0.018529 *lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 5) + -0.00228635 *lens_ipow(begin_x, 2)*lens_ipow(begin_lambda, 6) + -86.2577 *lens_ipow(begin_dy, 4)*lens_ipow(begin_lambda, 5) + -12.1588 *lens_ipow(begin_lambda, 10) + -231.177 *lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 8) + 0.000549329 *begin_x*lens_ipow(begin_y, 3)*begin_dx*begin_dy*lens_ipow(begin_lambda, 4)+0.0f;
    float invApJ[2][2];
    const float invdetap = 1.0f/(dx1_domega0[0][0]*dx1_domega0[1][1] - dx1_domega0[0][1]*dx1_domega0[1][0]);
    invApJ[0][0] =  dx1_domega0[1][1]*invdetap;
    invApJ[1][1] =  dx1_domega0[0][0]*invdetap;
    invApJ[0][1] = -dx1_domega0[0][1]*invdetap;
    invApJ[1][0] = -dx1_domega0[1][0]*invdetap;
    for(int i=0;i<2;i++)
    {
      dx += invApJ[0][i]*delta_ap[i];
      dy += invApJ[1][i]*delta_ap[i];
    }
    out[0] =  + 4.12859e-06  + 9.45701e-05 *begin_dy + 63.788 *begin_dx + 0.580638 *begin_x + 0.00691092 *begin_dx*begin_dy + 0.00593982 *lens_ipow(begin_dx, 2) + 3.54369e-05 *begin_y*begin_dy + 6.50137e-06 *begin_y*begin_dx + -2.53243e-06 *lens_ipow(begin_x, 2) + -42.1459 *begin_dx*lens_ipow(begin_dy, 2) + -41.5653 *lens_ipow(begin_dx, 3) + 0.0101871 *lens_ipow(begin_y, 2)*begin_dx + 0.402587 *begin_x*lens_ipow(begin_dy, 2) + 0.440776 *begin_x*lens_ipow(begin_dx, 2) + 0.0278479 *begin_x*begin_y*begin_dy + 0.000310157 *begin_x*lens_ipow(begin_y, 2) + -7.89379e-06 *lens_ipow(begin_x, 2)*begin_dy + 0.038931 *lens_ipow(begin_x, 2)*begin_dx + 0.000316531 *lens_ipow(begin_x, 3) + -0.167559 *begin_dx*lens_ipow(begin_dy, 3) + 0.10932 *begin_x*lens_ipow(begin_lambda, 3) + 5.43354e-05 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 2) + 6.98846 *begin_dx*lens_ipow(begin_lambda, 4) + -13.4067 *lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 2) + -0.00034174 *lens_ipow(begin_x, 2)*begin_y*begin_dx*begin_dy + 2.87462e-05 *lens_ipow(begin_x, 3)*lens_ipow(begin_dy, 2) + 6.7319e-08 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 2) + 8.47919e-06 *lens_ipow(begin_x, 3)*lens_ipow(begin_lambda, 3) + -0.861676 *begin_x*lens_ipow(begin_dy, 6) + 0.0117883 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 5) + -3.31947e-05 *lens_ipow(begin_y, 4)*lens_ipow(begin_dx, 5) + 6.70399e-08 *begin_x*lens_ipow(begin_y, 4)*lens_ipow(begin_lambda, 4) + 3.11338e-09 *lens_ipow(begin_x, 7)*lens_ipow(begin_dy, 2)*begin_lambda + -34.0239 *begin_dx*lens_ipow(begin_lambda, 10) + 73.4031 *lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 6) + -4.85342 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 9) + -0.588629 *begin_x*lens_ipow(begin_lambda, 10) + 10.5021 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 7)*lens_ipow(begin_dy, 2) + 1.21706e-15 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 8) + 1.94912e-16 *lens_ipow(begin_x, 11);
    out[1] =  + 9.84173e-05  + 63.7892 *begin_dy + -0.000241407 *begin_dx + 0.580787 *begin_y + -8.3677e-06 *begin_x + -0.00215718 *lens_ipow(begin_dy, 2) + 4.15587e-05 *begin_x*begin_dy + -7.87069e-07 *lens_ipow(begin_x, 2) + -41.6097 *lens_ipow(begin_dy, 3) + -42.3765 *lens_ipow(begin_dx, 2)*begin_dy + 0.442704 *begin_y*lens_ipow(begin_dy, 2) + 0.406335 *begin_y*lens_ipow(begin_dx, 2) + 0.0390684 *lens_ipow(begin_y, 2)*begin_dy + 0.000315021 *lens_ipow(begin_y, 3) + 0.0280539 *begin_x*begin_y*begin_dx + 2.18427e-07 *begin_x*lens_ipow(begin_y, 2) + 0.0101966 *lens_ipow(begin_x, 2)*begin_dy + 0.000313623 *lens_ipow(begin_x, 2)*begin_y + 0.107869 *begin_y*lens_ipow(begin_lambda, 3) + 2.88198e-05 *begin_x*begin_y*lens_ipow(begin_dx, 2) + 6.99505 *begin_dy*lens_ipow(begin_lambda, 4) + 1.82732e-05 *lens_ipow(begin_y, 3)*lens_ipow(begin_lambda, 2) + -0.000291516 *begin_x*lens_ipow(begin_y, 2)*begin_dx*begin_dy + 0.000136211 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dy, 2) + -3.33244e-06 *lens_ipow(begin_x, 3)*lens_ipow(begin_dy, 2) + 0.000454962 *begin_x*begin_dx*begin_dy*lens_ipow(begin_lambda, 3) + -1.13353e-05 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 3) + -1.3032 *begin_y*lens_ipow(begin_dx, 6) + 0.00721634 *lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 5) + 0.00126756 *lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 4) + -4.01578e-06 *begin_x*lens_ipow(begin_y, 3)*begin_dx*lens_ipow(begin_lambda, 2) + 4.2306e-09 *lens_ipow(begin_x, 6)*begin_dy*begin_lambda + 5.32008e-14 *lens_ipow(begin_y, 9) + 2.58723e-09 *lens_ipow(begin_y, 7)*lens_ipow(begin_dx, 2)*begin_lambda + -34.2711 *begin_dy*lens_ipow(begin_lambda, 10) + -0.583892 *begin_y*lens_ipow(begin_lambda, 10) + 4.9731 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 2) + -0.000408997 *lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 6) + 3.73e-15 *lens_ipow(begin_x, 6)*lens_ipow(begin_y, 5) + 7.46349e-13 *lens_ipow(begin_x, 8)*begin_y*lens_ipow(begin_lambda, 2);
    out[2] =  + 3.75124e-07  + -1.62788 *begin_dx + -4.36351e-09 *begin_y + -0.0305627 *begin_x + -1.14114 *begin_dx*lens_ipow(begin_dy, 2) + 0.960382 *lens_ipow(begin_dx, 3) + -0.0477622 *begin_y*begin_dx*begin_dy + -0.000394659 *lens_ipow(begin_y, 2)*begin_dx + -0.0334408 *begin_x*lens_ipow(begin_dy, 2) + -0.0017499 *begin_x*lens_ipow(begin_dx, 2) + -0.00115196 *begin_x*begin_y*begin_dy + -6.69053e-06 *begin_x*lens_ipow(begin_y, 2) + -0.00060646 *lens_ipow(begin_x, 2)*begin_dx + -3.10061e-06 *lens_ipow(begin_x, 3) + 0.0370781 *begin_dx*lens_ipow(begin_lambda, 3) + -0.000200509 *begin_x*lens_ipow(begin_lambda, 3) + 1.72885 *lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 2) + 1.57155e-05 *lens_ipow(begin_x, 2)*begin_y*begin_dx*begin_dy + -3.79636e-09 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 2)*begin_lambda + 0.000512223 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 2) + -6.4637e-05 *lens_ipow(begin_x, 3)*lens_ipow(begin_dy, 4) + 27.7906 *lens_ipow(begin_dx, 9) + -0.000317397 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 6) + -1.35798e-11 *lens_ipow(begin_x, 7)*lens_ipow(begin_dy, 2) + -0.242612 *begin_dx*lens_ipow(begin_lambda, 10) + -4.99303 *lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 6);
    out[3] =  + -1.08219e-06  + -1.6281 *begin_dy + 4.21623e-06 *begin_dx + -0.0305648 *begin_y + 2.02282e-05 *lens_ipow(begin_dy, 2) + 0.966154 *lens_ipow(begin_dy, 3) + 3.08821 *lens_ipow(begin_dx, 2)*begin_dy + -0.00137119 *begin_y*lens_ipow(begin_dy, 2) + 0.00642377 *begin_y*lens_ipow(begin_dx, 2) + -0.00059824 *lens_ipow(begin_y, 2)*begin_dy + -3.05808e-06 *lens_ipow(begin_y, 3) + 0.0735958 *begin_x*begin_dx*begin_dy + 0.000384833 *lens_ipow(begin_x, 2)*begin_dy + 7.84626e-07 *lens_ipow(begin_x, 2)*begin_y + 0.0379537 *begin_dy*lens_ipow(begin_lambda, 3) + -0.000203596 *begin_y*lens_ipow(begin_lambda, 3) + -6.60185e-09 *lens_ipow(begin_y, 5)*lens_ipow(begin_dx, 2) + -0.104594 *begin_x*lens_ipow(begin_dx, 5)*begin_dy + -0.0755948 *lens_ipow(begin_dx, 2)*begin_dy*lens_ipow(begin_lambda, 5) + 29.2819 *lens_ipow(begin_dy, 9) + 414.364 *lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 5) + -1.80799e-12 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 4)*begin_dy + 3.66008e-15 *lens_ipow(begin_x, 8)*begin_y + -0.250029 *begin_dy*lens_ipow(begin_lambda, 10);
    float pred_out_cs[7] = {0.0f};
    lens_sphereToCs(out, out+2, pred_out_cs, pred_out_cs+3, - lens_outer_pupil_curvature_radius, lens_outer_pupil_curvature_radius);
    float view[3] =
    {
      scene_x - pred_out_cs[0],
      scene_y - pred_out_cs[1],
      scene_z - pred_out_cs[2]
    };
    normalise(view);
    float out_new[5];
    lens_csToSphere(pred_out_cs, view, out_new, out_new+2, - lens_outer_pupil_curvature_radius, lens_outer_pupil_curvature_radius);
    const float delta_out[] = {out_new[2] - out[2], out_new[3] - out[3]};
    sqr_err = delta_out[0]*delta_out[0]+delta_out[1]*delta_out[1];
    float domega2_dx0[2][2];
    domega2_dx0[0][0] =  + -0.0305627  + -0.0334408 *lens_ipow(begin_dy, 2) + -0.0017499 *lens_ipow(begin_dx, 2) + -0.00115196 *begin_y*begin_dy + -6.69053e-06 *lens_ipow(begin_y, 2) + -0.00121292 *begin_x*begin_dx + -9.30182e-06 *lens_ipow(begin_x, 2) + -0.000200509 *lens_ipow(begin_lambda, 3) + 3.1431e-05 *begin_x*begin_y*begin_dx*begin_dy + -1.13891e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*begin_lambda + -0.000193911 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 4) + -0.000317397 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 6) + -9.50584e-11 *lens_ipow(begin_x, 6)*lens_ipow(begin_dy, 2)+0.0f;
    domega2_dx0[0][1] =  + -4.36351e-09  + -0.0477622 *begin_dx*begin_dy + -0.000789318 *begin_y*begin_dx + -0.00115196 *begin_x*begin_dy + -1.33811e-05 *begin_x*begin_y + 1.57155e-05 *lens_ipow(begin_x, 2)*begin_dx*begin_dy + -7.59273e-09 *lens_ipow(begin_x, 3)*begin_y*begin_lambda + 0.00102445 *begin_y*lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 2) + -0.000634794 *begin_x*begin_y*lens_ipow(begin_dx, 6)+0.0f;
    domega2_dx0[1][0] =  + 0.0735958 *begin_dx*begin_dy + 0.000769666 *begin_x*begin_dy + 1.56925e-06 *begin_x*begin_y + -0.104594 *lens_ipow(begin_dx, 5)*begin_dy + -7.23198e-12 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 4)*begin_dy + 2.92807e-14 *lens_ipow(begin_x, 7)*begin_y+0.0f;
    domega2_dx0[1][1] =  + -0.0305648  + -0.00137119 *lens_ipow(begin_dy, 2) + 0.00642377 *lens_ipow(begin_dx, 2) + -0.00119648 *begin_y*begin_dy + -9.17423e-06 *lens_ipow(begin_y, 2) + 7.84626e-07 *lens_ipow(begin_x, 2) + -0.000203596 *lens_ipow(begin_lambda, 3) + -3.30093e-08 *lens_ipow(begin_y, 4)*lens_ipow(begin_dx, 2) + -7.23198e-12 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 3)*begin_dy + 3.66008e-15 *lens_ipow(begin_x, 8)+0.0f;
    float invJ[2][2];
    const float invdet = 1.0f/(domega2_dx0[0][0]*domega2_dx0[1][1] - domega2_dx0[0][1]*domega2_dx0[1][0]);
    invJ[0][0] =  domega2_dx0[1][1]*invdet;
    invJ[1][1] =  domega2_dx0[0][0]*invdet;
    invJ[0][1] = -domega2_dx0[0][1]*invdet;
    invJ[1][0] = -domega2_dx0[1][0]*invdet;
    for(int i=0;i<2;i++)
    {
      x += invJ[0][i]*delta_out[i];
      y += invJ[1][i]*delta_out[i];
    }
    if(sqr_err>prev_sqr_err) error |= 1;
    if(sqr_ap_err>prev_sqr_ap_err) error |= 2;
    if(out[0]!=out[0]) error |= 4;
    if(out[1]!=out[1]) error |= 8;
    DEBUG_LOG;
    // reset error code for first few iterations.
    if(k<10) error = 0;
  }
}
else
  error = 128;
if(out[0]*out[0]+out[1]*out[1] > lens_outer_pupil_radius*lens_outer_pupil_radius) error |= 16;
const float begin_x = x;
const float begin_y = y;
const float begin_dx = dx;
const float begin_dy = dy;
const float begin_lambda = lambda;
if(error==0)
  out[4] =  + 0.695607  + 0.137175 *begin_lambda + 5.30727e-06 *begin_dx + 4.9533e-07 *begin_y + -0.0143904 *lens_ipow(begin_dy, 2) + -0.0113874 *lens_ipow(begin_dx, 2) + -0.000567518 *begin_y*begin_dy + -1.74105e-05 *lens_ipow(begin_y, 2) + -5.06366e-06 *begin_x*begin_dy + -0.000574729 *begin_x*begin_dx + -1.74568e-05 *lens_ipow(begin_x, 2) + -0.115482 *lens_ipow(begin_lambda, 3) + 1.14242e-05 *begin_y*begin_dx*begin_dy + -0.000129152 *begin_y*begin_dy*lens_ipow(begin_lambda, 2) + 0.000222906 *lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 2) + 8.82067e-05 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2) + -0.000178598 *begin_x*begin_dx*lens_ipow(begin_lambda, 2) + 0.000156553 *begin_x*begin_y*begin_dx*begin_dy + 8.91984e-05 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 2) + 0.000204498 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 2) + -1.11135e-10 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 3) + -7.68996 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 4) + -2.74706 *lens_ipow(begin_dx, 6) + -4.20595e-10 *lens_ipow(begin_y, 6) + -1.02151e-05 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 2) + -7.41889e-05 *lens_ipow(begin_x, 2)*begin_dx*lens_ipow(begin_dy, 3) + -1.16018e-09 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 4) + -1.17593e-09 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 2) + -4.42264e-10 *lens_ipow(begin_x, 6) + 0.00511412 *begin_x*begin_y*lens_ipow(begin_dx, 3)*begin_dy*begin_lambda + -17.8435 *lens_ipow(begin_dy, 8) + -58.0627 *lens_ipow(begin_dx, 6)*lens_ipow(begin_dy, 2) + 2.54104e-05 *begin_x*lens_ipow(begin_y, 3)*begin_dx*lens_ipow(begin_dy, 3) + 6.79814e-06 *lens_ipow(begin_x, 4)*lens_ipow(begin_dx, 4) + -5.7828e-11 *lens_ipow(begin_x, 6)*begin_y*begin_dy + -0.985058 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 5) + -0.186187 *lens_ipow(begin_dx, 4)*lens_ipow(begin_lambda, 5) + 1.5948e-08 *lens_ipow(begin_y, 6)*lens_ipow(begin_dy, 4) + -5.97847e-06 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 2) + 0.212409 *lens_ipow(begin_lambda, 11);
else
  out[4] = 0.0f;