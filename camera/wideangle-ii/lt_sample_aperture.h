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
       + -2.75777e-05  + 35.5802 *begin_dx + 9.16576e-06 *begin_y + 0.57611 *begin_x + 0.0045315 *begin_dx*begin_dy + 8.06954e-05 *begin_y*begin_dx + 2.53463e-06 *lens_ipow(begin_y, 2) + -7.80895e-07 *lens_ipow(begin_x, 2) + 38.1235 *begin_dx*lens_ipow(begin_dy, 2) + 38.2222 *lens_ipow(begin_dx, 3) + 3.0928 *begin_y*begin_dx*begin_dy + 0.0454737 *lens_ipow(begin_y, 2)*begin_dx + 1.19486 *begin_x*lens_ipow(begin_dy, 2) + 4.31218 *begin_x*lens_ipow(begin_dx, 2) + 0.0795652 *begin_x*begin_y*begin_dy + 5.09669e-05 *begin_x*begin_y*begin_dx + 0.00107686 *begin_x*lens_ipow(begin_y, 2) + 0.125814 *lens_ipow(begin_x, 2)*begin_dx + 0.00107984 *lens_ipow(begin_x, 3) + -0.0644539 *begin_x*lens_ipow(begin_lambda, 3) + -1.96618 *begin_dx*lens_ipow(begin_lambda, 4) + 0.000801215 *lens_ipow(begin_x, 2)*begin_y*begin_dx*begin_dy + 7.95431e-06 *lens_ipow(begin_x, 3)*begin_y*begin_dy + 1.38357e-08 *lens_ipow(begin_x, 4)*begin_y + 1.66823e-08 *lens_ipow(begin_x, 5) + -0.316163 *begin_y*lens_ipow(begin_dx, 4)*begin_dy + 7.54716e-08 *lens_ipow(begin_x, 4)*begin_y*begin_dy + 1.52972e-07 *begin_x*lens_ipow(begin_y, 5)*begin_dy*lens_ipow(begin_lambda, 2) + 3.87508e-12 *begin_x*lens_ipow(begin_y, 8) + -5.32515e-07 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 3)*begin_dx*begin_dy + 4.76218e-11 *lens_ipow(begin_x, 5)*lens_ipow(begin_y, 4) + 2.15057e-08 *lens_ipow(begin_x, 7)*lens_ipow(begin_dy, 2) + 51.781 *begin_x*begin_y*lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 3)*begin_lambda + 2.53133e-14 *lens_ipow(begin_x, 10) + 9.4202 *begin_dx*lens_ipow(begin_lambda, 10) + 0.342727 *begin_x*lens_ipow(begin_lambda, 10) + -0.00100112 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 3) + 4.27468e-08 *lens_ipow(begin_x, 6)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 3) + -3.60283e-10 *lens_ipow(begin_x, 8)*begin_y*begin_dx*begin_dy + 8.46017e-14 *lens_ipow(begin_x, 9)*lens_ipow(begin_y, 2),
       + -5.95104e-06  + 35.5792 *begin_dy + 0.57582 *begin_y + 4.90396e-06 *begin_x + 0.00412335 *lens_ipow(begin_dx, 2) + 0.00020035 *begin_y*begin_dx + -3.6962e-05 *begin_x*begin_dy + 37.9917 *lens_ipow(begin_dy, 3) + 38.0517 *lens_ipow(begin_dx, 2)*begin_dy + 4.30109 *begin_y*lens_ipow(begin_dy, 2) + 1.20187 *begin_y*lens_ipow(begin_dx, 2) + 0.125568 *lens_ipow(begin_y, 2)*begin_dy + 0.00107905 *lens_ipow(begin_y, 3) + 3.10618 *begin_x*begin_dx*begin_dy + 0.0801155 *begin_x*begin_y*begin_dx + 0.0455699 *lens_ipow(begin_x, 2)*begin_dy + 0.00108298 *lens_ipow(begin_x, 2)*begin_y + -0.0622474 *begin_y*lens_ipow(begin_lambda, 3) + 1.09019e-07 *begin_x*lens_ipow(begin_y, 3) + -1.88828 *begin_dy*lens_ipow(begin_lambda, 4) + 1.87954e-08 *lens_ipow(begin_y, 5) + 0.000867341 *begin_x*lens_ipow(begin_y, 2)*begin_dx*begin_dy + 3.22142e-06 *begin_x*lens_ipow(begin_y, 3)*begin_dx + 9.51235e-06 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*begin_dy + 1.96136 *begin_x*begin_y*begin_dx*lens_ipow(begin_dy, 4)*begin_lambda + -5.8194e-08 *begin_x*lens_ipow(begin_y, 6)*begin_dx*begin_dy + 0.000534861 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + -6.20265e-07 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 4)*begin_dx*begin_dy + 7.09764e-11 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 5) + -5.33428 *begin_x*lens_ipow(begin_dx, 3)*begin_dy*lens_ipow(begin_lambda, 5) + 0.00927332 *begin_x*lens_ipow(begin_y, 2)*begin_dx*begin_dy*lens_ipow(begin_lambda, 5) + -2.48897e-05 *lens_ipow(begin_x, 4)*begin_dy*lens_ipow(begin_lambda, 5) + 8.90388 *begin_dy*lens_ipow(begin_lambda, 10) + 0.323866 *begin_y*lens_ipow(begin_lambda, 10) + 7.43515e-06 *lens_ipow(begin_y, 6)*lens_ipow(begin_dx, 2)*begin_dy*lens_ipow(begin_lambda, 2) + 2.70745e-10 *lens_ipow(begin_y, 9)*lens_ipow(begin_dx, 2) + -0.000591592 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 3)*begin_dx*lens_ipow(begin_dy, 4) + -1.82172e-08 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 3)*lens_ipow(begin_lambda, 4) + 1.1258e-08 *lens_ipow(begin_x, 8)*lens_ipow(begin_dy, 3) + 4.68681e-10 *lens_ipow(begin_x, 8)*begin_y*lens_ipow(begin_dy, 2)
    };
    const float delta_ap[] = {ap_x - pred_ap[0], ap_y - pred_ap[1]};
    sqr_ap_err = delta_ap[0]*delta_ap[0]+delta_ap[1]*delta_ap[1];
    float dx1_domega0[2][2];
    dx1_domega0[0][0] =  + 35.5802  + 0.0045315 *begin_dy + 8.06954e-05 *begin_y + 38.1235 *lens_ipow(begin_dy, 2) + 114.667 *lens_ipow(begin_dx, 2) + 3.0928 *begin_y*begin_dy + 0.0454737 *lens_ipow(begin_y, 2) + 8.62436 *begin_x*begin_dx + 5.09669e-05 *begin_x*begin_y + 0.125814 *lens_ipow(begin_x, 2) + -1.96618 *lens_ipow(begin_lambda, 4) + 0.000801215 *lens_ipow(begin_x, 2)*begin_y*begin_dy + -1.26465 *begin_y*lens_ipow(begin_dx, 3)*begin_dy + -5.32515e-07 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 3)*begin_dy + 207.124 *begin_x*begin_y*lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 3)*begin_lambda + 9.4202 *lens_ipow(begin_lambda, 10) + -0.00200224 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 3)*begin_dx*lens_ipow(begin_dy, 3) + 1.2824e-07 *lens_ipow(begin_x, 6)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2) + -3.60283e-10 *lens_ipow(begin_x, 8)*begin_y*begin_dy+0.0f;
    dx1_domega0[0][1] =  + 0.0045315 *begin_dx + 76.247 *begin_dx*begin_dy + 3.0928 *begin_y*begin_dx + 2.38972 *begin_x*begin_dy + 0.0795652 *begin_x*begin_y + 0.000801215 *lens_ipow(begin_x, 2)*begin_y*begin_dx + 7.95431e-06 *lens_ipow(begin_x, 3)*begin_y + -0.316163 *begin_y*lens_ipow(begin_dx, 4) + 7.54716e-08 *lens_ipow(begin_x, 4)*begin_y + 1.52972e-07 *begin_x*lens_ipow(begin_y, 5)*lens_ipow(begin_lambda, 2) + -5.32515e-07 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 3)*begin_dx + 4.30114e-08 *lens_ipow(begin_x, 7)*begin_dy + 155.343 *begin_x*begin_y*lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 2)*begin_lambda + -0.00300336 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + -3.60283e-10 *lens_ipow(begin_x, 8)*begin_y*begin_dx+0.0f;
    dx1_domega0[1][0] =  + 0.0082467 *begin_dx + 0.00020035 *begin_y + 76.1035 *begin_dx*begin_dy + 2.40374 *begin_y*begin_dx + 3.10618 *begin_x*begin_dy + 0.0801155 *begin_x*begin_y + 0.000867341 *begin_x*lens_ipow(begin_y, 2)*begin_dy + 3.22142e-06 *begin_x*lens_ipow(begin_y, 3) + 1.96136 *begin_x*begin_y*lens_ipow(begin_dy, 4)*begin_lambda + -5.8194e-08 *begin_x*lens_ipow(begin_y, 6)*begin_dy + 0.00106972 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 3)*begin_dx*lens_ipow(begin_dy, 2) + -6.20265e-07 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 4)*begin_dy + -16.0029 *begin_x*lens_ipow(begin_dx, 2)*begin_dy*lens_ipow(begin_lambda, 5) + 0.00927332 *begin_x*lens_ipow(begin_y, 2)*begin_dy*lens_ipow(begin_lambda, 5) + 1.48703e-05 *lens_ipow(begin_y, 6)*begin_dx*begin_dy*lens_ipow(begin_lambda, 2) + 5.4149e-10 *lens_ipow(begin_y, 9)*begin_dx + -0.000591592 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 4)+0.0f;
    dx1_domega0[1][1] =  + 35.5792  + -3.6962e-05 *begin_x + 113.975 *lens_ipow(begin_dy, 2) + 38.0517 *lens_ipow(begin_dx, 2) + 8.60218 *begin_y*begin_dy + 0.125568 *lens_ipow(begin_y, 2) + 3.10618 *begin_x*begin_dx + 0.0455699 *lens_ipow(begin_x, 2) + -1.88828 *lens_ipow(begin_lambda, 4) + 0.000867341 *begin_x*lens_ipow(begin_y, 2)*begin_dx + 9.51235e-06 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2) + 7.84544 *begin_x*begin_y*begin_dx*lens_ipow(begin_dy, 3)*begin_lambda + -5.8194e-08 *begin_x*lens_ipow(begin_y, 6)*begin_dx + 0.00106972 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 2)*begin_dy + -6.20265e-07 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 4)*begin_dx + -5.33428 *begin_x*lens_ipow(begin_dx, 3)*lens_ipow(begin_lambda, 5) + 0.00927332 *begin_x*lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_lambda, 5) + -2.48897e-05 *lens_ipow(begin_x, 4)*lens_ipow(begin_lambda, 5) + 8.90388 *lens_ipow(begin_lambda, 10) + 7.43515e-06 *lens_ipow(begin_y, 6)*lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 2) + -0.00236637 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 3)*begin_dx*lens_ipow(begin_dy, 3) + 3.37739e-08 *lens_ipow(begin_x, 8)*lens_ipow(begin_dy, 2) + 9.37362e-10 *lens_ipow(begin_x, 8)*begin_y*begin_dy+0.0f;
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
    out[0] =  + 0.000273071  + 11.3942 *begin_dx + -1.79206 *begin_x + -0.026403 *lens_ipow(begin_dy, 2) + 0.000508213 *begin_y*begin_dx + 18.3599 *begin_dx*lens_ipow(begin_dy, 2) + 1.40919 *begin_y*begin_dx*begin_dy + 0.051411 *lens_ipow(begin_y, 2)*begin_dx + 0.702177 *begin_x*lens_ipow(begin_dy, 2) + 1.96079 *begin_x*lens_ipow(begin_dx, 2) + 0.0634463 *begin_x*begin_y*begin_dy + 0.00273604 *begin_x*lens_ipow(begin_y, 2) + 0.114242 *lens_ipow(begin_x, 2)*begin_dx + 0.00274454 *lens_ipow(begin_x, 3) + -0.0747238 *begin_y*lens_ipow(begin_dx, 2)*begin_dy + -0.126033 *begin_x*lens_ipow(begin_lambda, 3) + 60.2941 *lens_ipow(begin_dx, 3)*lens_ipow(begin_lambda, 2) + 0.00172066 *lens_ipow(begin_x, 2)*begin_dx*lens_ipow(begin_lambda, 2) + -0.00806772 *lens_ipow(begin_x, 2)*begin_y*begin_dx*begin_dy + 0.00313786 *lens_ipow(begin_x, 3)*lens_ipow(begin_dy, 2) + 4.07219 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 2) + -15.3263 *begin_x*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 2) + -0.00904287 *begin_x*lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 2)*begin_dy + 0.74991 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 5) + 6.27855e-05 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2) + -12.738 *begin_y*begin_dx*lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 3) + 0.0457966 *lens_ipow(begin_y, 4)*lens_ipow(begin_dx, 5) + 6.7248e-08 *lens_ipow(begin_y, 7)*begin_dx*begin_dy + -0.0273947 *begin_x*begin_y*begin_dy*lens_ipow(begin_lambda, 6) + -4.36807e-11 *begin_x*lens_ipow(begin_y, 8) + -0.0206066 *lens_ipow(begin_x, 3)*begin_y*lens_ipow(begin_dy, 5) + -2.85077e-10 *lens_ipow(begin_x, 5)*lens_ipow(begin_y, 4) + -1.08046e-05 *lens_ipow(begin_x, 6)*lens_ipow(begin_dx, 3) + -9.21904e-11 *lens_ipow(begin_x, 7)*lens_ipow(begin_y, 2) + -2.95374e-11 *lens_ipow(begin_x, 9) + -2.35678e+06 *lens_ipow(begin_dx, 5)*lens_ipow(begin_dy, 6) + -7.89248 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 5)*lens_ipow(begin_lambda, 4) + 0.526558 *begin_x*lens_ipow(begin_lambda, 10) + 8.03647 *begin_x*lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 8) + 2.45229e-05 *lens_ipow(begin_x, 6)*begin_y*lens_ipow(begin_dx, 3)*begin_dy;
    out[1] =  + 0.000173604  + 11.0908 *begin_dy + -1.78474 *begin_y + -0.000642049 *begin_x*begin_dy + 18.7614 *lens_ipow(begin_dy, 3) + 18.3002 *lens_ipow(begin_dx, 2)*begin_dy + -0.103219 *begin_y*lens_ipow(begin_lambda, 2) + 2.16252 *begin_y*lens_ipow(begin_dy, 2) + 0.687603 *begin_y*lens_ipow(begin_dx, 2) + 0.12054 *lens_ipow(begin_y, 2)*begin_dy + 8.00862e-05 *lens_ipow(begin_y, 2)*begin_dx + 0.00277844 *lens_ipow(begin_y, 3) + 1.42244 *begin_x*begin_dx*begin_dy + 0.070137 *begin_x*begin_y*begin_dx + 0.0504552 *lens_ipow(begin_x, 2)*begin_dy + 0.00283847 *lens_ipow(begin_x, 2)*begin_y + 1.57168 *begin_dy*lens_ipow(begin_lambda, 3) + -0.0108884 *lens_ipow(begin_y, 2)*begin_dy*lens_ipow(begin_lambda, 2) + 0.00471391 *lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 2) + -0.012487 *begin_x*begin_y*begin_dx*lens_ipow(begin_lambda, 2) + -0.349981 *begin_x*begin_y*begin_dx*lens_ipow(begin_dy, 2) + -0.0642196 *begin_x*begin_y*lens_ipow(begin_dx, 3) + -0.00822533 *begin_x*lens_ipow(begin_y, 2)*begin_dx*begin_dy + 0.189502 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 3) + -19.79 *begin_y*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 2) + -10.9848 *begin_x*lens_ipow(begin_dx, 3)*begin_dy*lens_ipow(begin_lambda, 2) + 0.000115866 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 2) + -2.96425e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 5) + -2.39147e-08 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 3) + -8.56936e-06 *lens_ipow(begin_y, 6)*lens_ipow(begin_dy, 3) + -1.20812e-07 *lens_ipow(begin_y, 7)*lens_ipow(begin_dx, 2) + -2.76624e-11 *lens_ipow(begin_y, 9) + -1.40083e-08 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 5)*begin_dx + -3.22486e-06 *lens_ipow(begin_x, 4)*begin_y*lens_ipow(begin_lambda, 4) + -3.87192e-05 *lens_ipow(begin_x, 5)*begin_y*begin_dx*lens_ipow(begin_dy, 2) + -3.21587e-11 *lens_ipow(begin_x, 8)*begin_y + -101.879 *begin_y*lens_ipow(begin_dy, 6)*lens_ipow(begin_lambda, 3) + 0.582681 *begin_y*lens_ipow(begin_lambda, 10) + -2.17193e-07 *begin_x*lens_ipow(begin_y, 7)*begin_dx*lens_ipow(begin_dy, 2) + -0.0611717 *lens_ipow(begin_x, 2)*begin_dy*lens_ipow(begin_lambda, 8);
    out[2] =  + 3.1777e-06  + -0.38828 *begin_dx + -0.0269329 *begin_x + 4.96973e-08 *begin_x*begin_y + -2.32394e-07 *lens_ipow(begin_x, 2) + 0.443308 *begin_dx*lens_ipow(begin_dy, 2) + 0.40462 *lens_ipow(begin_dx, 3) + 0.00861408 *begin_y*begin_dx*begin_dy + 0.000199611 *lens_ipow(begin_y, 2)*begin_dx + -0.000182299 *begin_x*begin_y*begin_dy + -3.7398e-06 *begin_x*lens_ipow(begin_y, 2) + -0.000142886 *lens_ipow(begin_x, 2)*begin_dx + 1.71332e-05 *lens_ipow(begin_x, 3) + 0.00518719 *begin_x*lens_ipow(begin_lambda, 3) + 0.0153536 *begin_dx*lens_ipow(begin_lambda, 4) + -0.0166865 *begin_x*lens_ipow(begin_dy, 4) + -5.03297e-05 *lens_ipow(begin_x, 2)*begin_y*begin_dx*begin_dy + 3.38484e-05 *lens_ipow(begin_x, 3)*lens_ipow(begin_dy, 2) + -9.23377e-07 *lens_ipow(begin_x, 3)*begin_y*begin_dy + -6.88028e-09 *lens_ipow(begin_x, 4)*begin_y*begin_dy + 0.0165216 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 5) + -9.26389e-05 *begin_x*lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 2)*begin_dy + 3.84119e-05 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 2) + -1.35264e-05 *lens_ipow(begin_x, 3)*begin_y*lens_ipow(begin_dy, 3) + 7.3658e-07 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2) + -7.77223e-06 *lens_ipow(begin_x, 4)*lens_ipow(begin_dx, 3) + 0.000167283 *lens_ipow(begin_x, 3)*lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 3) + -1054.51 *lens_ipow(begin_dx, 5)*lens_ipow(begin_dy, 4) + 4.06104e-06 *lens_ipow(begin_y, 5)*begin_dx*lens_ipow(begin_dy, 3) + -1.62985e-05 *lens_ipow(begin_y, 5)*lens_ipow(begin_dx, 3)*begin_dy + -2.32588e-09 *lens_ipow(begin_x, 7)*lens_ipow(begin_dy, 2) + -1.30349e-12 *lens_ipow(begin_x, 7)*lens_ipow(begin_y, 2) + -3.9334e-13 *lens_ipow(begin_x, 9) + -0.593565 *begin_x*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 5) + -6.93166e-08 *lens_ipow(begin_x, 5)*lens_ipow(begin_lambda, 5) + -1778.97 *lens_ipow(begin_dx, 11) + -0.0255423 *begin_x*lens_ipow(begin_lambda, 10) + -4.11955 *begin_x*lens_ipow(begin_dx, 6)*lens_ipow(begin_lambda, 4) + 3.20863e-06 *begin_x*lens_ipow(begin_y, 3)*begin_dy*lens_ipow(begin_lambda, 6) + 4.6848e-06 *lens_ipow(begin_x, 4)*begin_y*begin_dx*begin_dy*lens_ipow(begin_lambda, 4);
    out[3] =  + 3.45882e-06  + -0.387565 *begin_dy + -1.97645e-05 *begin_dx + -0.0268871 *begin_y + -1.91681e-07 *begin_x*begin_y + 0.39321 *lens_ipow(begin_dy, 3) + 0.379793 *lens_ipow(begin_dx, 2)*begin_dy + -0.00767475 *begin_y*lens_ipow(begin_dx, 2) + -0.000185387 *lens_ipow(begin_y, 2)*begin_dy + 1.65477e-05 *lens_ipow(begin_y, 3) + -3.06053e-05 *begin_x*begin_y*begin_dx + -3.09542e-08 *begin_x*lens_ipow(begin_y, 2) + -0.000210542 *lens_ipow(begin_x, 2)*begin_dy + 3.9166e-05 *lens_ipow(begin_x, 2)*begin_y + 0.0051492 *begin_y*lens_ipow(begin_lambda, 3) + 0.0173095 *begin_dy*lens_ipow(begin_lambda, 4) + -0.00126511 *begin_x*begin_y*begin_dx*lens_ipow(begin_dy, 2) + 4.3721e-05 *lens_ipow(begin_x, 3)*begin_dx*begin_dy + -0.00265214 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 3)*begin_dy + -4.80643e-05 *begin_x*lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 3) + 1.05682e-06 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 2) + -5.01031e-10 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 5) + -5.54933e-05 *lens_ipow(begin_x, 3)*begin_y*begin_dx*lens_ipow(begin_dy, 2) + -6.29819e-10 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 3) + -3.8479e-10 *lens_ipow(begin_x, 6)*begin_y + 12.5116 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 4) + 8.60045e-07 *lens_ipow(begin_y, 5)*lens_ipow(begin_dy, 4) + -4.62414e-13 *lens_ipow(begin_y, 9) + -3.58501e-07 *begin_x*lens_ipow(begin_y, 5)*begin_dx*lens_ipow(begin_dy, 2) + 1.96721e-07 *lens_ipow(begin_x, 6)*lens_ipow(begin_dy, 3) + -0.154681 *begin_y*lens_ipow(begin_dy, 4)*lens_ipow(begin_lambda, 5) + 1.41152e-06 *lens_ipow(begin_y, 5)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 3) + 4.07165e-08 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 4)*begin_dx*begin_dy*begin_lambda + -8.0891e-08 *lens_ipow(begin_x, 4)*begin_y*lens_ipow(begin_lambda, 5) + -26752.4 *lens_ipow(begin_dx, 6)*lens_ipow(begin_dy, 5) + -0.0247506 *begin_y*lens_ipow(begin_lambda, 10) + -9.41248e-08 *lens_ipow(begin_y, 5)*lens_ipow(begin_lambda, 6) + -0.000591273 *lens_ipow(begin_x, 2)*begin_dy*lens_ipow(begin_lambda, 8) + 1.21852e-10 *lens_ipow(begin_x, 7)*lens_ipow(begin_y, 2)*begin_dx*begin_dy + 4.86351e-11 *lens_ipow(begin_x, 8)*begin_y*lens_ipow(begin_dx, 2);
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
    domega2_dx0[0][0] =  + -0.0269329  + 4.96973e-08 *begin_y + -4.64787e-07 *begin_x + -0.000182299 *begin_y*begin_dy + -3.7398e-06 *lens_ipow(begin_y, 2) + -0.000285772 *begin_x*begin_dx + 5.13997e-05 *lens_ipow(begin_x, 2) + 0.00518719 *lens_ipow(begin_lambda, 3) + -0.0166865 *lens_ipow(begin_dy, 4) + -0.000100659 *begin_x*begin_y*begin_dx*begin_dy + 0.000101545 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 2) + -2.77013e-06 *lens_ipow(begin_x, 2)*begin_y*begin_dy + -2.75211e-08 *lens_ipow(begin_x, 3)*begin_y*begin_dy + -9.26389e-05 *lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 2)*begin_dy + 7.68238e-05 *begin_x*lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 2) + -4.05792e-05 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dy, 3) + 2.20974e-06 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2) + -3.10889e-05 *lens_ipow(begin_x, 3)*lens_ipow(begin_dx, 3) + 0.00050185 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 3) + -1.62812e-08 *lens_ipow(begin_x, 6)*lens_ipow(begin_dy, 2) + -9.12446e-12 *lens_ipow(begin_x, 6)*lens_ipow(begin_y, 2) + -3.54006e-12 *lens_ipow(begin_x, 8) + -0.593565 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 5) + -3.46583e-07 *lens_ipow(begin_x, 4)*lens_ipow(begin_lambda, 5) + -0.0255423 *lens_ipow(begin_lambda, 10) + -4.11955 *lens_ipow(begin_dx, 6)*lens_ipow(begin_lambda, 4) + 3.20863e-06 *lens_ipow(begin_y, 3)*begin_dy*lens_ipow(begin_lambda, 6) + 1.87392e-05 *lens_ipow(begin_x, 3)*begin_y*begin_dx*begin_dy*lens_ipow(begin_lambda, 4)+0.0f;
    domega2_dx0[0][1] =  + 4.96973e-08 *begin_x + 0.00861408 *begin_dx*begin_dy + 0.000399222 *begin_y*begin_dx + -0.000182299 *begin_x*begin_dy + -7.47959e-06 *begin_x*begin_y + -5.03297e-05 *lens_ipow(begin_x, 2)*begin_dx*begin_dy + -9.23377e-07 *lens_ipow(begin_x, 3)*begin_dy + -6.88028e-09 *lens_ipow(begin_x, 4)*begin_dy + 0.0330432 *begin_y*lens_ipow(begin_dx, 5) + -0.000277917 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*begin_dy + 7.68238e-05 *lens_ipow(begin_x, 2)*begin_y*begin_dx*lens_ipow(begin_dy, 2) + -1.35264e-05 *lens_ipow(begin_x, 3)*lens_ipow(begin_dy, 3) + 1.47316e-06 *lens_ipow(begin_x, 3)*begin_y*lens_ipow(begin_dx, 2) + 2.03052e-05 *lens_ipow(begin_y, 4)*begin_dx*lens_ipow(begin_dy, 3) + -8.14927e-05 *lens_ipow(begin_y, 4)*lens_ipow(begin_dx, 3)*begin_dy + -2.60699e-12 *lens_ipow(begin_x, 7)*begin_y + 9.62589e-06 *begin_x*lens_ipow(begin_y, 2)*begin_dy*lens_ipow(begin_lambda, 6) + 4.6848e-06 *lens_ipow(begin_x, 4)*begin_dx*begin_dy*lens_ipow(begin_lambda, 4)+0.0f;
    domega2_dx0[1][0] =  + -1.91681e-07 *begin_y + -3.06053e-05 *begin_y*begin_dx + -3.09542e-08 *lens_ipow(begin_y, 2) + -0.000421085 *begin_x*begin_dy + 7.83319e-05 *begin_x*begin_y + -0.00126511 *begin_y*begin_dx*lens_ipow(begin_dy, 2) + 0.000131163 *lens_ipow(begin_x, 2)*begin_dx*begin_dy + -0.00265214 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 3)*begin_dy + -4.80643e-05 *lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 3) + 2.11365e-06 *begin_x*lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 2) + -1.00206e-09 *begin_x*lens_ipow(begin_y, 5) + -0.00016648 *lens_ipow(begin_x, 2)*begin_y*begin_dx*lens_ipow(begin_dy, 2) + -2.51928e-09 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 3) + -2.30874e-09 *lens_ipow(begin_x, 5)*begin_y + -3.58501e-07 *lens_ipow(begin_y, 5)*begin_dx*lens_ipow(begin_dy, 2) + 1.18033e-06 *lens_ipow(begin_x, 5)*lens_ipow(begin_dy, 3) + 1.22149e-07 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 4)*begin_dx*begin_dy*begin_lambda + -3.23564e-07 *lens_ipow(begin_x, 3)*begin_y*lens_ipow(begin_lambda, 5) + -0.00118255 *begin_x*begin_dy*lens_ipow(begin_lambda, 8) + 8.52964e-10 *lens_ipow(begin_x, 6)*lens_ipow(begin_y, 2)*begin_dx*begin_dy + 3.89081e-10 *lens_ipow(begin_x, 7)*begin_y*lens_ipow(begin_dx, 2)+0.0f;
    domega2_dx0[1][1] =  + -0.0268871  + -1.91681e-07 *begin_x + -0.00767475 *lens_ipow(begin_dx, 2) + -0.000370774 *begin_y*begin_dy + 4.9643e-05 *lens_ipow(begin_y, 2) + -3.06053e-05 *begin_x*begin_dx + -6.19085e-08 *begin_x*begin_y + 3.9166e-05 *lens_ipow(begin_x, 2) + 0.0051492 *lens_ipow(begin_lambda, 3) + -0.00126511 *begin_x*begin_dx*lens_ipow(begin_dy, 2) + -0.00530428 *begin_x*begin_y*lens_ipow(begin_dx, 3)*begin_dy + -0.000144193 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 3) + 3.17047e-06 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 2) + -2.50516e-09 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 4) + -5.54933e-05 *lens_ipow(begin_x, 3)*begin_dx*lens_ipow(begin_dy, 2) + -1.88946e-09 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 2) + -3.8479e-10 *lens_ipow(begin_x, 6) + 4.30022e-06 *lens_ipow(begin_y, 4)*lens_ipow(begin_dy, 4) + -4.16173e-12 *lens_ipow(begin_y, 8) + -1.79251e-06 *begin_x*lens_ipow(begin_y, 4)*begin_dx*lens_ipow(begin_dy, 2) + -0.154681 *lens_ipow(begin_dy, 4)*lens_ipow(begin_lambda, 5) + 7.05762e-06 *lens_ipow(begin_y, 4)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 3) + 1.62866e-07 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 3)*begin_dx*begin_dy*begin_lambda + -8.0891e-08 *lens_ipow(begin_x, 4)*lens_ipow(begin_lambda, 5) + -0.0247506 *lens_ipow(begin_lambda, 10) + -4.70624e-07 *lens_ipow(begin_y, 4)*lens_ipow(begin_lambda, 6) + 2.43704e-10 *lens_ipow(begin_x, 7)*begin_y*begin_dx*begin_dy + 4.86351e-11 *lens_ipow(begin_x, 8)*lens_ipow(begin_dx, 2)+0.0f;
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
  out[4] =  + 0.122746  + 0.246012 *begin_lambda + 1.81291e-07 *begin_x + -0.0545617 *lens_ipow(begin_dy, 2) + -0.0577175 *lens_ipow(begin_dx, 2) + -0.00124334 *begin_y*begin_dy + -1.09109e-05 *lens_ipow(begin_y, 2) + -0.00104348 *begin_x*begin_dx + -1.1659e-07 *begin_x*begin_y + -1.1266e-05 *lens_ipow(begin_x, 2) + -0.201334 *lens_ipow(begin_lambda, 3) + -0.0470381 *begin_y*lens_ipow(begin_dx, 2)*begin_dy + -7.83989e-08 *lens_ipow(begin_y, 4) + -0.0469432 *begin_x*begin_dx*lens_ipow(begin_dy, 2) + -0.024014 *begin_x*lens_ipow(begin_dx, 3) + -0.000889007 *begin_x*begin_y*begin_dx*begin_dy + -1.2894e-07 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2) + -6.67779e-08 *lens_ipow(begin_x, 4) + -15.2367 *lens_ipow(begin_dy, 6) + -59.9699 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 4) + -95.6303 *lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 2) + -19.4174 *lens_ipow(begin_dx, 6) + -0.000198097 *lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 3) + -1.23354e-06 *lens_ipow(begin_y, 4)*lens_ipow(begin_dx, 2) + -0.00028129 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 3) + -0.0170956 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 4) + -0.000684665 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dy, 3) + -2.78293e-05 *lens_ipow(begin_x, 3)*lens_ipow(begin_dx, 3) + 7.80966e-07 *lens_ipow(begin_x, 4)*lens_ipow(begin_dy, 2) + -0.16413 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 6) + 1.22364 *begin_x*begin_y*lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 3) + -0.000843017 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + -6.90113e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 4)*lens_ipow(begin_dy, 2) + -9.96442e-09 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2) + -286.151 *begin_x*lens_ipow(begin_dx, 7)*lens_ipow(begin_dy, 2) + -0.0658487 *begin_x*lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 6) + 0.00531327 *begin_x*lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 5)*begin_dy + 0.00489262 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dx, 6)*begin_dy + 1.6838e-09 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 5)*begin_dx*begin_dy + 0.355518 *lens_ipow(begin_lambda, 11);
else
  out[4] = 0.0f;