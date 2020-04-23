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
       + 0.000146105  + 119.392 *begin_dx + -4.05857e-06 *begin_y + 1.247 *begin_x + -0.000155495 *begin_y*begin_dx + -7.6004e-07 *lens_ipow(begin_y, 2) + -0.0470823 *begin_x*begin_lambda + -0.000269889 *begin_x*begin_dy + 0.000205308 *begin_x*begin_dx + -2.41965e-06 *begin_x*begin_y + 315.426 *begin_dx*lens_ipow(begin_dy, 2) + 312.895 *lens_ipow(begin_dx, 3) + 5.80165 *begin_y*begin_dx*begin_dy + 0.025389 *lens_ipow(begin_y, 2)*begin_dx + 2.95132 *begin_x*lens_ipow(begin_dy, 2) + 8.54125 *begin_x*lens_ipow(begin_dx, 2) + 0.0525601 *begin_x*begin_y*begin_dy + 0.00021685 *begin_x*lens_ipow(begin_y, 2) + 0.0772709 *lens_ipow(begin_x, 2)*begin_dx + 0.000215934 *lens_ipow(begin_x, 3) + -0.0139484 *begin_x*begin_dx*lens_ipow(begin_dy, 2) + 0.379694 *begin_x*lens_ipow(begin_dx, 2)*begin_lambda + -0.000975206 *begin_x*begin_y*begin_dy*begin_lambda + 0.00264109 *begin_x*lens_ipow(begin_lambda, 4) + 0.655091 *begin_x*lens_ipow(begin_dy, 4) + -0.0239011 *begin_x*begin_y*lens_ipow(begin_dx, 2)*begin_dy + -10.6975 *begin_dx*lens_ipow(begin_lambda, 5) + 5784.08 *lens_ipow(begin_dx, 5)*lens_ipow(begin_dy, 2) + 0.003015 *begin_x*begin_dy*lens_ipow(begin_lambda, 5) + -5.51715e-05 *lens_ipow(begin_x, 4)*lens_ipow(begin_dx, 3) + 0.0233978 *begin_y*begin_dx*begin_dy*lens_ipow(begin_lambda, 5) + 1115.85 *lens_ipow(begin_dx, 5)*lens_ipow(begin_lambda, 4) + 7923.29 *lens_ipow(begin_dx, 9) + 0.0392799 *begin_x*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 6) + -2.15026e-05 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_lambda, 6) + -6.55007e-14 *lens_ipow(begin_x, 5)*lens_ipow(begin_y, 4) + 1.20246e-08 *lens_ipow(begin_x, 6)*begin_dx*lens_ipow(begin_lambda, 3) + 42.4803 *begin_dx*lens_ipow(begin_lambda, 10) + 403790 *lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 6)*lens_ipow(begin_lambda, 2) + 0.141664 *begin_x*lens_ipow(begin_lambda, 10),
       + 0.000109425  + 119.438 *begin_dy + 0.00155992 *begin_dx + 1.23169 *begin_y + 1.63555e-05 *begin_x + 0.00746521 *lens_ipow(begin_dx, 2) + 0.000168035 *begin_y*begin_dy + -0.000243009 *begin_y*begin_dx + -1.17086e-06 *lens_ipow(begin_x, 2) + 314.444 *lens_ipow(begin_dy, 3) + 312.195 *lens_ipow(begin_dx, 2)*begin_dy + 8.78967 *begin_y*lens_ipow(begin_dy, 2) + 0.00178997 *begin_y*begin_dx*begin_dy + 2.98505 *begin_y*lens_ipow(begin_dx, 2) + 0.0781152 *lens_ipow(begin_y, 2)*begin_dy + 0.00022001 *lens_ipow(begin_y, 3) + 5.76403 *begin_x*begin_dx*begin_dy + 0.0524156 *begin_x*begin_y*begin_dx + -1.65767e-07 *begin_x*lens_ipow(begin_y, 2) + 0.0251438 *lens_ipow(begin_x, 2)*begin_dy + 0.000217566 *lens_ipow(begin_x, 2)*begin_y + -0.0660435 *begin_y*lens_ipow(begin_lambda, 3) + -1.05783e-08 *lens_ipow(begin_x, 3)*begin_y + -0.000238607 *lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 2) + 0.000734854 *begin_x*lens_ipow(begin_y, 2)*begin_dx*begin_dy + 3.20796e-06 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*begin_dy + -12.3294 *begin_dy*lens_ipow(begin_lambda, 5) + 0.0040916 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 4) + 11026.7 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 5) + -4.69828e-05 *lens_ipow(begin_y, 4)*lens_ipow(begin_dy, 3) + -0.00574666 *lens_ipow(begin_y, 2)*begin_dy*lens_ipow(begin_lambda, 5) + -4.14292e-08 *lens_ipow(begin_y, 5)*lens_ipow(begin_lambda, 3) + 42394.4 *lens_ipow(begin_dy, 9) + 1162.52 *lens_ipow(begin_dx, 4)*begin_dy*lens_ipow(begin_lambda, 4) + 9.79145 *begin_y*lens_ipow(begin_dx, 4)*lens_ipow(begin_lambda, 4) + 919.141 *lens_ipow(begin_dy, 5)*lens_ipow(begin_lambda, 5) + 51.7455 *begin_dy*lens_ipow(begin_lambda, 10) + 0.348236 *begin_y*lens_ipow(begin_lambda, 10) + -0.0370741 *lens_ipow(begin_x, 4)*lens_ipow(begin_dy, 7) + -9.08606e-17 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 7)
    };
    const float delta_ap[] = {ap_x - pred_ap[0], ap_y - pred_ap[1]};
    sqr_ap_err = delta_ap[0]*delta_ap[0]+delta_ap[1]*delta_ap[1];
    float dx1_domega0[2][2];
    dx1_domega0[0][0] =  + 119.392  + -0.000155495 *begin_y + 0.000205308 *begin_x + 315.426 *lens_ipow(begin_dy, 2) + 938.685 *lens_ipow(begin_dx, 2) + 5.80165 *begin_y*begin_dy + 0.025389 *lens_ipow(begin_y, 2) + 17.0825 *begin_x*begin_dx + 0.0772709 *lens_ipow(begin_x, 2) + -0.0139484 *begin_x*lens_ipow(begin_dy, 2) + 0.759388 *begin_x*begin_dx*begin_lambda + -0.0478022 *begin_x*begin_y*begin_dx*begin_dy + -10.6975 *lens_ipow(begin_lambda, 5) + 28920.4 *lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 2) + -0.000165514 *lens_ipow(begin_x, 4)*lens_ipow(begin_dx, 2) + 0.0233978 *begin_y*begin_dy*lens_ipow(begin_lambda, 5) + 5579.25 *lens_ipow(begin_dx, 4)*lens_ipow(begin_lambda, 4) + 71309.6 *lens_ipow(begin_dx, 8) + 1.20246e-08 *lens_ipow(begin_x, 6)*lens_ipow(begin_lambda, 3) + 42.4803 *lens_ipow(begin_lambda, 10) + 1.21137e+06 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 6)*lens_ipow(begin_lambda, 2)+0.0f;
    dx1_domega0[0][1] =  + -0.000269889 *begin_x + 630.851 *begin_dx*begin_dy + 5.80165 *begin_y*begin_dx + 5.90264 *begin_x*begin_dy + 0.0525601 *begin_x*begin_y + -0.0278968 *begin_x*begin_dx*begin_dy + -0.000975206 *begin_x*begin_y*begin_lambda + 2.62037 *begin_x*lens_ipow(begin_dy, 3) + -0.0239011 *begin_x*begin_y*lens_ipow(begin_dx, 2) + 11568.2 *lens_ipow(begin_dx, 5)*begin_dy + 0.003015 *begin_x*lens_ipow(begin_lambda, 5) + 0.0233978 *begin_y*begin_dx*lens_ipow(begin_lambda, 5) + 0.0785598 *begin_x*begin_dy*lens_ipow(begin_lambda, 6) + 2.42274e+06 *lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 5)*lens_ipow(begin_lambda, 2)+0.0f;
    dx1_domega0[1][0] =  + 0.00155992  + 0.0149304 *begin_dx + -0.000243009 *begin_y + 624.39 *begin_dx*begin_dy + 0.00178997 *begin_y*begin_dy + 5.97009 *begin_y*begin_dx + 5.76403 *begin_x*begin_dy + 0.0524156 *begin_x*begin_y + -0.000477215 *lens_ipow(begin_y, 3)*begin_dx + 0.000734854 *begin_x*lens_ipow(begin_y, 2)*begin_dy + 22053.4 *begin_dx*lens_ipow(begin_dy, 5) + 4650.06 *lens_ipow(begin_dx, 3)*begin_dy*lens_ipow(begin_lambda, 4) + 39.1658 *begin_y*lens_ipow(begin_dx, 3)*lens_ipow(begin_lambda, 4)+0.0f;
    dx1_domega0[1][1] =  + 119.438  + 0.000168035 *begin_y + 943.331 *lens_ipow(begin_dy, 2) + 312.195 *lens_ipow(begin_dx, 2) + 17.5793 *begin_y*begin_dy + 0.00178997 *begin_y*begin_dx + 0.0781152 *lens_ipow(begin_y, 2) + 5.76403 *begin_x*begin_dx + 0.0251438 *lens_ipow(begin_x, 2) + 0.000734854 *begin_x*lens_ipow(begin_y, 2)*begin_dx + 3.20796e-06 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2) + -12.3294 *lens_ipow(begin_lambda, 5) + 0.0163664 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 3) + 55133.6 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 4) + -0.000140948 *lens_ipow(begin_y, 4)*lens_ipow(begin_dy, 2) + -0.00574666 *lens_ipow(begin_y, 2)*lens_ipow(begin_lambda, 5) + 381550 *lens_ipow(begin_dy, 8) + 1162.52 *lens_ipow(begin_dx, 4)*lens_ipow(begin_lambda, 4) + 4595.71 *lens_ipow(begin_dy, 4)*lens_ipow(begin_lambda, 5) + 51.7455 *lens_ipow(begin_lambda, 10) + -0.259519 *lens_ipow(begin_x, 4)*lens_ipow(begin_dy, 6)+0.0f;
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
    out[0] =  + -7.53719e-05  + -0.0040198 *begin_dy + 399.583 *begin_dx + -2.28631e-05 *begin_y + 2.58883 *begin_x + 0.0282439 *lens_ipow(begin_dx, 2) + 2.41907e-06 *begin_x*begin_y + 5.21373e-06 *lens_ipow(begin_x, 2) + 1.70625 *begin_dx*lens_ipow(begin_lambda, 2) + -241.438 *begin_dx*lens_ipow(begin_dy, 2) + -245.154 *lens_ipow(begin_dx, 3) + -1.18733 *begin_y*begin_dx*begin_dy + 0.0114028 *lens_ipow(begin_y, 2)*begin_dx + 1.26004 *begin_x*lens_ipow(begin_dy, 2) + 0.0414114 *begin_x*begin_y*begin_dy + 0.000350368 *begin_x*lens_ipow(begin_y, 2) + 0.0576284 *lens_ipow(begin_x, 2)*begin_dx + 0.000350798 *lens_ipow(begin_x, 3) + -0.0874267 *begin_x*lens_ipow(begin_lambda, 3) + 0.0227785 *begin_x*lens_ipow(begin_dx, 2)*begin_lambda + 0.0131072 *begin_x*begin_y*begin_dy*begin_lambda + 2.20274e-06 *begin_x*lens_ipow(begin_y, 2)*begin_dy + 343.293 *lens_ipow(begin_dx, 5) + -3.40182e-05 *lens_ipow(begin_x, 2)*lens_ipow(begin_lambda, 3) + 0.00275759 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 3)*begin_lambda + -106.895 *lens_ipow(begin_dx, 3)*lens_ipow(begin_lambda, 4) + 29.8208 *begin_y*lens_ipow(begin_dx, 3)*begin_dy*lens_ipow(begin_lambda, 2) + 1.00066 *begin_x*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 4) + -2.10546e-06 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 3)*begin_dx*begin_dy + 0.0177254 *lens_ipow(begin_x, 3)*lens_ipow(begin_dx, 4) + 0.000103079 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_lambda, 5) + 3.66702e-05 *lens_ipow(begin_x, 4)*begin_dx*lens_ipow(begin_lambda, 3) + -2.0833 *lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 4)*lens_ipow(begin_lambda, 2) + -6.11315 *begin_dx*lens_ipow(begin_lambda, 9) + 0.277279 *begin_x*lens_ipow(begin_lambda, 9) + 2943.8 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 4)*lens_ipow(begin_lambda, 5) + 9.99168e+06 *lens_ipow(begin_dx, 9)*lens_ipow(begin_dy, 2) + 0.719103 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 3)*lens_ipow(begin_lambda, 6) + 4.37915e-07 *lens_ipow(begin_x, 5)*lens_ipow(begin_lambda, 6) + 6.93607e-16 *lens_ipow(begin_x, 7)*lens_ipow(begin_y, 4);
    out[1] =  + -6.14701e-05  + 398.883 *begin_dy + 2.58684 *begin_y + 1.48627e-05 *begin_x + 0.000427641 *begin_y*begin_dy + -2.87632e-06 *begin_x*begin_y + 1.79306e-06 *lens_ipow(begin_x, 2) + 5.16056 *begin_dy*lens_ipow(begin_lambda, 2) + -241.665 *lens_ipow(begin_dy, 3) + -290.011 *lens_ipow(begin_dx, 2)*begin_dy + 0.413533 *begin_y*lens_ipow(begin_dx, 2) + 0.057092 *lens_ipow(begin_y, 2)*begin_dy + 0.000349051 *lens_ipow(begin_y, 3) + -3.05602 *begin_x*begin_dx*begin_dy + 5.3719e-05 *begin_x*begin_y*begin_dy + 0.00877699 *lens_ipow(begin_x, 2)*begin_dy + 0.000301458 *lens_ipow(begin_x, 2)*begin_y + -0.0599654 *begin_y*lens_ipow(begin_lambda, 3) + 0.0330832 *begin_y*lens_ipow(begin_dy, 2)*begin_lambda + -0.019729 *begin_y*lens_ipow(begin_dx, 3) + 0.0821105 *begin_x*begin_y*begin_dx*begin_lambda + 199.497 *lens_ipow(begin_dy, 5) + 5.50017 *begin_x*begin_dx*begin_dy*lens_ipow(begin_lambda, 2) + -0.00051614 *begin_x*lens_ipow(begin_y, 2)*begin_dx*begin_dy + 5.21794e-06 *lens_ipow(begin_x, 3)*begin_y*begin_dx + -128.211 *lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 4) + 7.10664 *begin_y*lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 4) + 160.905 *begin_y*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 4) + 52.3091 *begin_y*lens_ipow(begin_dx, 6) + -0.000133896 *lens_ipow(begin_y, 4)*lens_ipow(begin_dy, 3) + -0.0104561 *begin_x*begin_y*lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 2) + 0.000427871 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_lambda, 4) + 6.65948e-11 *lens_ipow(begin_x, 6)*begin_y + 4.79792e-05 *lens_ipow(begin_y, 4)*begin_dy*lens_ipow(begin_lambda, 3) + 2.88671e-09 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 3)*lens_ipow(begin_lambda, 4) + -40.335 *begin_dy*lens_ipow(begin_lambda, 9) + 1958.2 *lens_ipow(begin_dx, 2)*begin_dy*lens_ipow(begin_lambda, 8) + 5.94329e-07 *lens_ipow(begin_y, 5)*lens_ipow(begin_lambda, 6) + 0.136651 *lens_ipow(begin_x, 2)*begin_dy*lens_ipow(begin_lambda, 8) + 4.95212e-16 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 9);
    out[2] =  + -1.36984e-06  + -3.02386 *begin_dx + -0.0220831 *begin_x + -12.6295 *begin_dx*lens_ipow(begin_dy, 2) + 1.85207 *lens_ipow(begin_dx, 3) + -0.181349 *begin_y*begin_dx*begin_dy + -0.00070382 *lens_ipow(begin_y, 2)*begin_dx + -0.105079 *begin_x*lens_ipow(begin_dy, 2) + 0.00846378 *begin_x*lens_ipow(begin_dx, 2) + -0.0016196 *begin_x*begin_y*begin_dy + -6.67979e-06 *begin_x*lens_ipow(begin_y, 2) + -0.000318603 *lens_ipow(begin_x, 2)*begin_dx + -2.13191e-06 *lens_ipow(begin_x, 3) + 0.000778602 *begin_x*lens_ipow(begin_lambda, 3) + -0.298228 *lens_ipow(begin_dx, 5) + -0.000632288 *begin_y*lens_ipow(begin_dx, 3)*begin_dy + 8.03208e-09 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 4) + -0.00411936 *begin_x*lens_ipow(begin_lambda, 10);
    out[3] =  + -1.39664e-06  + -3.02497 *begin_dy + -5.30828e-06 *begin_dx + -0.0220656 *begin_y + 0.000282626 *lens_ipow(begin_dy, 2) + 2.79523e-06 *begin_y*begin_dx + -2.11973e-06 *begin_x*begin_dy + 1.89904 *lens_ipow(begin_dy, 3) + 16.361 *lens_ipow(begin_dx, 2)*begin_dy + 0.00984536 *begin_y*lens_ipow(begin_dy, 2) + 0.0845495 *begin_y*lens_ipow(begin_dx, 2) + -0.000303261 *lens_ipow(begin_y, 2)*begin_dy + -2.07305e-06 *lens_ipow(begin_y, 3) + 1.20882e-06 *begin_x*lens_ipow(begin_dy, 2) + 0.223632 *begin_x*begin_dx*begin_dy + 0.00106874 *begin_x*begin_y*begin_dx + 0.000690967 *lens_ipow(begin_x, 2)*begin_dy + 2.47702e-06 *lens_ipow(begin_x, 2)*begin_y + 5.98681e-07 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2) + -0.00870921 *begin_x*begin_dx*begin_dy*begin_lambda + -0.000127605 *begin_x*begin_y*begin_dx*begin_lambda + 0.00111463 *begin_y*lens_ipow(begin_lambda, 4) + 0.0171457 *begin_x*begin_dx*lens_ipow(begin_dy, 3) + 2.16915e-09 *begin_x*lens_ipow(begin_y, 3)*begin_dx + -0.840573 *lens_ipow(begin_dx, 2)*begin_dy*lens_ipow(begin_lambda, 3) + 3.8148 *lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 3) + -0.0128798 *begin_y*lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 4) + 0.0016008 *begin_y*lens_ipow(begin_dx, 6) + -6.72972e-05 *lens_ipow(begin_x, 2)*begin_dy*lens_ipow(begin_lambda, 4) + -2.33936 *lens_ipow(begin_dy, 5)*lens_ipow(begin_lambda, 4) + 20.5479 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 4) + -1.4105e-06 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_lambda, 6) + -0.00165138 *begin_dy*lens_ipow(begin_lambda, 9) + 12.7052 *lens_ipow(begin_dx, 4)*begin_dy*lens_ipow(begin_lambda, 5) + -0.00501469 *begin_y*lens_ipow(begin_lambda, 10) + 0.000271308 *lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 4)*lens_ipow(begin_lambda, 4);
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
    domega2_dx0[0][0] =  + -0.0220831  + -0.105079 *lens_ipow(begin_dy, 2) + 0.00846378 *lens_ipow(begin_dx, 2) + -0.0016196 *begin_y*begin_dy + -6.67979e-06 *lens_ipow(begin_y, 2) + -0.000637207 *begin_x*begin_dx + -6.39572e-06 *lens_ipow(begin_x, 2) + 0.000778602 *lens_ipow(begin_lambda, 3) + 2.40962e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 4) + -0.00411936 *lens_ipow(begin_lambda, 10)+0.0f;
    domega2_dx0[0][1] =  + -0.181349 *begin_dx*begin_dy + -0.00140764 *begin_y*begin_dx + -0.0016196 *begin_x*begin_dy + -1.33596e-05 *begin_x*begin_y + -0.000632288 *lens_ipow(begin_dx, 3)*begin_dy + 1.60642e-08 *lens_ipow(begin_x, 3)*begin_y*lens_ipow(begin_dx, 4)+0.0f;
    domega2_dx0[1][0] =  + -2.11973e-06 *begin_dy + 1.20882e-06 *lens_ipow(begin_dy, 2) + 0.223632 *begin_dx*begin_dy + 0.00106874 *begin_y*begin_dx + 0.00138193 *begin_x*begin_dy + 4.95404e-06 *begin_x*begin_y + -0.00870921 *begin_dx*begin_dy*begin_lambda + -0.000127605 *begin_y*begin_dx*begin_lambda + 0.0171457 *begin_dx*lens_ipow(begin_dy, 3) + 2.16915e-09 *lens_ipow(begin_y, 3)*begin_dx + -0.000134594 *begin_x*begin_dy*lens_ipow(begin_lambda, 4) + -2.82099e-06 *begin_x*begin_y*lens_ipow(begin_lambda, 6)+0.0f;
    domega2_dx0[1][1] =  + -0.0220656  + 2.79523e-06 *begin_dx + 0.00984536 *lens_ipow(begin_dy, 2) + 0.0845495 *lens_ipow(begin_dx, 2) + -0.000606523 *begin_y*begin_dy + -6.21915e-06 *lens_ipow(begin_y, 2) + 0.00106874 *begin_x*begin_dx + 2.47702e-06 *lens_ipow(begin_x, 2) + 1.19736e-06 *begin_y*lens_ipow(begin_dx, 2) + -0.000127605 *begin_x*begin_dx*begin_lambda + 0.00111463 *lens_ipow(begin_lambda, 4) + 6.50746e-09 *begin_x*lens_ipow(begin_y, 2)*begin_dx + -0.0128798 *lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 4) + 0.0016008 *lens_ipow(begin_dx, 6) + -1.4105e-06 *lens_ipow(begin_x, 2)*lens_ipow(begin_lambda, 6) + -0.00501469 *lens_ipow(begin_lambda, 10) + 0.000813925 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 4)*lens_ipow(begin_lambda, 4)+0.0f;
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
  out[4] =  + 0.416426  + 0.271139 *begin_lambda + 1.25411e-05 *begin_dy + -9.02847e-05 *begin_dx + 2.92475e-07 *begin_y + -4.35994e-07 *begin_x + -0.143453 *lens_ipow(begin_dy, 2) + -0.134492 *lens_ipow(begin_dx, 2) + -0.00255901 *begin_y*begin_dy + -7.55823e-06 *begin_y*begin_dx + -1.1374e-05 *lens_ipow(begin_y, 2) + -0.00243245 *begin_x*begin_dx + -1.16659e-05 *lens_ipow(begin_x, 2) + -0.225852 *lens_ipow(begin_lambda, 3) + 0.00454149 *lens_ipow(begin_dx, 3) + 6.74301e-05 *begin_y*lens_ipow(begin_dx, 2) + -4.28026e-09 *lens_ipow(begin_x, 2)*begin_y + 0.808027 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + 0.0666357 *begin_y*lens_ipow(begin_dy, 3) + -0.000179446 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2) + -4.81253e-06 *lens_ipow(begin_y, 3)*begin_dy + 0.000739319 *begin_x*lens_ipow(begin_dx, 2)*begin_dy + 0.0388081 *begin_x*lens_ipow(begin_dx, 3) + -5.01914e-06 *begin_x*lens_ipow(begin_y, 2)*begin_dx + 0.000312941 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 2) + -2.94753e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2) + -3.96148e-06 *lens_ipow(begin_x, 3)*begin_dx + 0.105056 *begin_x*begin_dx*lens_ipow(begin_dy, 2)*begin_lambda + 134.378 *lens_ipow(begin_dy, 6) + 38.0779 *lens_ipow(begin_dx, 6) + -6.74165e-11 *lens_ipow(begin_y, 6) + -4.71779e-11 *lens_ipow(begin_x, 6) + -8.04499e-06 *lens_ipow(begin_x, 2)*begin_y*begin_dy*lens_ipow(begin_lambda, 5) + 514.727 *begin_y*lens_ipow(begin_dy, 9) + -3.58428e-14 *lens_ipow(begin_y, 9)*begin_dy + -0.000707711 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 2) + -0.000128013 *lens_ipow(begin_x, 4)*lens_ipow(begin_dy, 6) + -2.22617e-14 *lens_ipow(begin_x, 8)*begin_y*begin_dy + 0.404952 *lens_ipow(begin_lambda, 11) + 15659.6 *lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 4)*lens_ipow(begin_lambda, 3);
else
  out[4] = 0.0f;
