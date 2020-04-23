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
       + 30 *begin_dx + 1 *begin_x,
       + 30 *begin_dy + 1 *begin_y
    };
    const float delta_ap[] = {ap_x - pred_ap[0], ap_y - pred_ap[1]};
    sqr_ap_err = delta_ap[0]*delta_ap[0]+delta_ap[1]*delta_ap[1];
    float dx1_domega0[2][2];
    dx1_domega0[0][0] =  + 30 +0.0f;
    dx1_domega0[0][1] = +0.0f;
    dx1_domega0[1][0] = +0.0f;
    dx1_domega0[1][1] =  + 30 +0.0f;
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
    out[0] =  + 0.000332875  + 38.4939 *begin_dx + 5.52046e-07 *begin_y + 0.808035 *begin_x + -0.00232839 *lens_ipow(begin_dx, 2) + -1.56745e-06 *begin_x*begin_y + 1.13603 *begin_dx*lens_ipow(begin_lambda, 2) + -8.94836 *begin_dx*lens_ipow(begin_dy, 2) + -11.2887 *lens_ipow(begin_dx, 3) + 0.00325261 *lens_ipow(begin_y, 2)*begin_dx + 0.14674 *begin_x*lens_ipow(begin_dy, 2) + 0.0112868 *begin_x*begin_y*begin_dy + 0.000141138 *begin_x*lens_ipow(begin_y, 2) + 0.0121795 *lens_ipow(begin_x, 2)*begin_dx + 0.000129755 *lens_ipow(begin_x, 3) + 0.0221718 *begin_dx*lens_ipow(begin_dy, 3) + -0.0212165 *lens_ipow(begin_dx, 3)*begin_dy + 0.0157898 *begin_x*lens_ipow(begin_lambda, 3) + 0.00665478 *lens_ipow(begin_dx, 4)*begin_dy + 3.20226 *lens_ipow(begin_dx, 5) + -0.0818127 *begin_y*lens_ipow(begin_dx, 3)*begin_dy + -0.00352902 *lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 2) + 0.000239042 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2) + 0.000639126 *lens_ipow(begin_x, 2)*begin_dx*lens_ipow(begin_lambda, 2) + 0.000274387 *lens_ipow(begin_x, 3)*lens_ipow(begin_dx, 2) + 9.3389 *begin_dx*lens_ipow(begin_dy, 6) + 16.4145 *lens_ipow(begin_dx, 5)*lens_ipow(begin_dy, 2) + -1.36058 *begin_dx*lens_ipow(begin_lambda, 7) + -0.70917 *lens_ipow(begin_dx, 3)*lens_ipow(begin_lambda, 5) + 0.0168168 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 7) + 0.00866087 *begin_x*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 6) + 3.03344e-11 *lens_ipow(begin_x, 8)*begin_dx + 3.07216 *lens_ipow(begin_dx, 11) + 0.445589 *begin_y*begin_dx*lens_ipow(begin_dy, 9) + -4.43038e-08 *lens_ipow(begin_y, 6)*lens_ipow(begin_dx, 5) + 3.85742e-14 *lens_ipow(begin_y, 10)*begin_dx + 0.00375095 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dx, 7)*begin_dy + -4.35115e-11 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 5)*lens_ipow(begin_dy, 3) + -6.51474e-12 *lens_ipow(begin_x, 5)*lens_ipow(begin_y, 4)*lens_ipow(begin_dx, 2) + 1.81549e-10 *lens_ipow(begin_x, 7)*lens_ipow(begin_lambda, 4);
    out[1] =  + 0.000800444  + 38.6818 *begin_dy + -0.000396378 *begin_dx + 0.809858 *begin_y + 0.000112878 *begin_y*begin_dy + -4.9063e-05 *begin_x*begin_dy + 0.000166233 *begin_x*begin_dx + -10.6648 *lens_ipow(begin_dy, 3) + -9.49893 *lens_ipow(begin_dx, 2)*begin_dy + 0.141577 *begin_y*lens_ipow(begin_dx, 2) + 0.0109229 *lens_ipow(begin_y, 2)*begin_dy + 0.000115124 *lens_ipow(begin_y, 3) + 0.0113582 *begin_x*begin_y*begin_dx + 0.00352415 *lens_ipow(begin_x, 2)*begin_dy + 0.000148077 *lens_ipow(begin_x, 2)*begin_y + 1.22435 *begin_dy*lens_ipow(begin_lambda, 4) + 4.07587 *lens_ipow(begin_dx, 4)*begin_dy + 0.000111513 *lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 2) + -0.0332369 *begin_x*begin_dx*lens_ipow(begin_dy, 3) + -3.93619e-05 *begin_x*lens_ipow(begin_y, 2)*begin_dx*begin_dy + -0.0037682 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 2)*begin_dy + 0.000233835 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dy, 2) + 4.8376 *lens_ipow(begin_dy, 5)*lens_ipow(begin_lambda, 2) + 17.9149 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 5) + 0.00825841 *begin_y*lens_ipow(begin_dx, 6) + 8.54151e-08 *lens_ipow(begin_y, 5)*lens_ipow(begin_dx, 2) + 0.314717 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 6) + 6.11751 *lens_ipow(begin_dy, 9) + 8.02097 *lens_ipow(begin_dx, 6)*lens_ipow(begin_dy, 3) + 0.0666773 *begin_y*lens_ipow(begin_lambda, 8) + 0.000170091 *lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 4) + 0.0116883 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 7) + 2.43312e-08 *lens_ipow(begin_y, 6)*begin_dy*lens_ipow(begin_lambda, 3) + -10.395 *lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 8) + 0.00108051 *lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 6)*lens_ipow(begin_lambda, 2) + 0.00313065 *begin_x*lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 7) + -7.99358e-05 *lens_ipow(begin_x, 3)*begin_y*lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 4) + -6.593e-12 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 5)*lens_ipow(begin_dy, 2) + -4.07579e-07 *lens_ipow(begin_x, 5)*begin_dx*begin_dy*lens_ipow(begin_lambda, 4) + -1.65797e-11 *lens_ipow(begin_x, 5)*lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 3);
    out[2] =  + 4.01593e-05  + -1.13552 *begin_dx + 1.71298e-09 *begin_y + -0.0492188 *begin_x + -0.000230811 *lens_ipow(begin_dy, 2) + -8.67819e-05 *lens_ipow(begin_dx, 2) + 1.70815e-06 *begin_x*begin_dy + -2.70486e-07 *lens_ipow(begin_x, 2) + -1.98465 *begin_dx*lens_ipow(begin_dy, 2) + -1.25098 *lens_ipow(begin_dx, 3) + -0.0933151 *begin_y*begin_dx*begin_dy + -0.00108289 *lens_ipow(begin_y, 2)*begin_dx + -0.0504185 *begin_x*lens_ipow(begin_dy, 2) + -0.0819684 *begin_x*lens_ipow(begin_dx, 2) + -0.00250089 *begin_x*begin_y*begin_dy + -2.63877e-05 *begin_x*lens_ipow(begin_y, 2) + -0.00192777 *lens_ipow(begin_x, 2)*begin_dx + -1.22139e-05 *lens_ipow(begin_x, 3) + 0.169351 *begin_dx*lens_ipow(begin_lambda, 3) + -0.0027988 *begin_y*lens_ipow(begin_dx, 3)*begin_dy + -2.69104e-05 *begin_y*lens_ipow(begin_dx, 4) + 0.00795007 *begin_x*lens_ipow(begin_lambda, 4) + -0.00104939 *begin_x*lens_ipow(begin_dy, 4) + -1.22159e-05 *lens_ipow(begin_x, 2)*begin_y*begin_dx*begin_dy + 5.53087e-09 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 2) + 6.28782e-06 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 2) + -0.197021 *lens_ipow(begin_dx, 5)*lens_ipow(begin_dy, 2) + 5.19268e-07 *lens_ipow(begin_y, 4)*begin_dx*lens_ipow(begin_dy, 2) + 2.50031e-06 *lens_ipow(begin_x, 3)*lens_ipow(begin_dy, 4) + -1.90139e-05 *lens_ipow(begin_x, 3)*lens_ipow(begin_dx, 4) + -0.080177 *lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 6) + -9.89596e-13 *lens_ipow(begin_y, 8)*begin_dx + -1.09702e-12 *lens_ipow(begin_x, 8)*begin_dx + -0.879262 *begin_dx*lens_ipow(begin_lambda, 10) + -0.137924 *begin_dx*lens_ipow(begin_dy, 10) + -0.224137 *lens_ipow(begin_dx, 11) + -0.000313522 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 9) + -0.0384944 *begin_x*lens_ipow(begin_lambda, 10) + 0.00145553 *begin_x*begin_y*lens_ipow(begin_dx, 6)*begin_dy*lens_ipow(begin_lambda, 2) + 5.33249e-08 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 4);
    out[3] =  + -3.49293e-06  + -1.13523 *begin_dy + -0.0493509 *begin_y + 2.05829e-06 *begin_x + -0.000170224 *lens_ipow(begin_dy, 2) + 1.65345e-07 *lens_ipow(begin_y, 2) + 1.37864e-08 *begin_x*begin_y + -1.19016 *lens_ipow(begin_dy, 3) + -0.444308 *lens_ipow(begin_dx, 2)*begin_dy + 0.000185801 *lens_ipow(begin_dx, 3) + -0.0762816 *begin_y*lens_ipow(begin_dy, 2) + -0.0176375 *begin_y*lens_ipow(begin_dx, 2) + -0.00175014 *lens_ipow(begin_y, 2)*begin_dy + -1.03955e-05 *lens_ipow(begin_y, 3) + 0.00379349 *begin_x*begin_dx*begin_dy + -0.000454984 *begin_x*begin_y*begin_dx + 0.000305372 *lens_ipow(begin_x, 2)*begin_dy + 2.89778e-06 *lens_ipow(begin_x, 2)*begin_y + -0.000653953 *lens_ipow(begin_dx, 3)*begin_dy + 0.265059 *begin_dy*lens_ipow(begin_lambda, 4) + 0.0080274 *begin_y*lens_ipow(begin_lambda, 4) + -0.000148504 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*begin_dy + 0.00103283 *begin_x*begin_y*begin_dx*lens_ipow(begin_dy, 2) + -8.20623e-07 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*begin_dy + -1.00139e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 3) + -0.456593 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 5) + -2.0609e-05 *lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 4) + 0.000710661 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 4)*begin_dy + -5.06788e-07 *lens_ipow(begin_x, 2)*begin_y*begin_dx*lens_ipow(begin_dy, 4) + -0.475785 *lens_ipow(begin_dx, 6)*lens_ipow(begin_dy, 3) + -0.334401 *lens_ipow(begin_dx, 8)*begin_dy + -1.2027e-12 *lens_ipow(begin_y, 8)*begin_dy + 0.00068186 *begin_x*begin_y*lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 4) + -5.29948e-10 *lens_ipow(begin_x, 6)*lens_ipow(begin_dy, 3) + -1.27064 *begin_dy*lens_ipow(begin_lambda, 10) + -0.250263 *lens_ipow(begin_dy, 11) + -0.0386125 *begin_y*lens_ipow(begin_lambda, 10) + -0.00220053 *begin_y*lens_ipow(begin_dx, 10) + 0.000429511 *begin_x*begin_y*begin_dx*lens_ipow(begin_dy, 8) + -1.84131e-15 *lens_ipow(begin_x, 9)*begin_y*begin_dx;
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
    domega2_dx0[0][0] =  + -0.0492188  + 1.70815e-06 *begin_dy + -5.40972e-07 *begin_x + -0.0504185 *lens_ipow(begin_dy, 2) + -0.0819684 *lens_ipow(begin_dx, 2) + -0.00250089 *begin_y*begin_dy + -2.63877e-05 *lens_ipow(begin_y, 2) + -0.00385553 *begin_x*begin_dx + -3.66417e-05 *lens_ipow(begin_x, 2) + 0.00795007 *lens_ipow(begin_lambda, 4) + -0.00104939 *lens_ipow(begin_dy, 4) + -2.44319e-05 *begin_x*begin_y*begin_dx*begin_dy + 1.65926e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2) + 1.25756e-05 *begin_x*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 2) + 7.50092e-06 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 4) + -5.70417e-05 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 4) + -8.77615e-12 *lens_ipow(begin_x, 7)*begin_dx + -0.0384944 *lens_ipow(begin_lambda, 10) + 0.00145553 *begin_y*lens_ipow(begin_dx, 6)*begin_dy*lens_ipow(begin_lambda, 2) + 1.59975e-07 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 4)+0.0f;
    domega2_dx0[0][1] =  + 1.71298e-09  + -0.0933151 *begin_dx*begin_dy + -0.00216578 *begin_y*begin_dx + -0.00250089 *begin_x*begin_dy + -5.27754e-05 *begin_x*begin_y + -0.0027988 *lens_ipow(begin_dx, 3)*begin_dy + -2.69104e-05 *lens_ipow(begin_dx, 4) + -1.22159e-05 *lens_ipow(begin_x, 2)*begin_dx*begin_dy + 1.10617e-08 *lens_ipow(begin_x, 3)*begin_y + 2.07707e-06 *lens_ipow(begin_y, 3)*begin_dx*lens_ipow(begin_dy, 2) + -7.91676e-12 *lens_ipow(begin_y, 7)*begin_dx + -0.000627045 *begin_y*lens_ipow(begin_dx, 9) + 0.00145553 *begin_x*lens_ipow(begin_dx, 6)*begin_dy*lens_ipow(begin_lambda, 2) + 1.0665e-07 *lens_ipow(begin_x, 3)*begin_y*lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 4)+0.0f;
    domega2_dx0[1][0] =  + 2.05829e-06  + 1.37864e-08 *begin_y + 0.00379349 *begin_dx*begin_dy + -0.000454984 *begin_y*begin_dx + 0.000610744 *begin_x*begin_dy + 5.79556e-06 *begin_x*begin_y + 0.00103283 *begin_y*begin_dx*lens_ipow(begin_dy, 2) + -1.64125e-06 *begin_x*lens_ipow(begin_y, 2)*begin_dy + -2.00279e-08 *begin_x*lens_ipow(begin_y, 3) + 0.00142132 *begin_x*lens_ipow(begin_dx, 4)*begin_dy + -1.01358e-06 *begin_x*begin_y*begin_dx*lens_ipow(begin_dy, 4) + 0.00068186 *begin_y*lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 4) + -3.17969e-09 *lens_ipow(begin_x, 5)*lens_ipow(begin_dy, 3) + 0.000429511 *begin_y*begin_dx*lens_ipow(begin_dy, 8) + -1.65718e-14 *lens_ipow(begin_x, 8)*begin_y*begin_dx+0.0f;
    domega2_dx0[1][1] =  + -0.0493509  + 3.3069e-07 *begin_y + 1.37864e-08 *begin_x + -0.0762816 *lens_ipow(begin_dy, 2) + -0.0176375 *lens_ipow(begin_dx, 2) + -0.00350028 *begin_y*begin_dy + -3.11864e-05 *lens_ipow(begin_y, 2) + -0.000454984 *begin_x*begin_dx + 2.89778e-06 *lens_ipow(begin_x, 2) + 0.0080274 *lens_ipow(begin_lambda, 4) + -0.000297008 *begin_y*lens_ipow(begin_dx, 2)*begin_dy + 0.00103283 *begin_x*begin_dx*lens_ipow(begin_dy, 2) + -1.64125e-06 *lens_ipow(begin_x, 2)*begin_y*begin_dy + -3.00418e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2) + -6.18269e-05 *lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 4) + -5.06788e-07 *lens_ipow(begin_x, 2)*begin_dx*lens_ipow(begin_dy, 4) + -9.62162e-12 *lens_ipow(begin_y, 7)*begin_dy + 0.00068186 *begin_x*lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 4) + -0.0386125 *lens_ipow(begin_lambda, 10) + -0.00220053 *lens_ipow(begin_dx, 10) + 0.000429511 *begin_x*begin_dx*lens_ipow(begin_dy, 8) + -1.84131e-15 *lens_ipow(begin_x, 9)*begin_dx+0.0f;
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
  out[4] =  + 0.917368  + -5.66651e-05 *begin_dy + 2.32758e-05 *begin_dx + -0.0368658 *lens_ipow(begin_dy, 2) + -0.0222116 *lens_ipow(begin_dx, 2) + -0.00153216 *begin_y*begin_dy + -1.8982e-05 *lens_ipow(begin_y, 2) + -9.97992e-07 *begin_x*begin_dy + -0.00112255 *begin_x*begin_dx + -8.63335e-06 *lens_ipow(begin_x, 2) + 0.0254024 *lens_ipow(begin_lambda, 3) + -0.000667485 *begin_dx*lens_ipow(begin_dy, 2) + -0.000946979 *begin_dx*lens_ipow(begin_dy, 3) + -0.241391 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + -0.129633 *lens_ipow(begin_dx, 4) + -2.63693e-05 *begin_y*lens_ipow(begin_dx, 3) + 3.60707e-05 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2) + 0.00037023 *begin_x*begin_y*begin_dx*begin_dy + -1.45543e-07 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2) + 0.000126749 *lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 2)*begin_lambda + -0.00679596 *begin_x*lens_ipow(begin_dx, 3)*begin_lambda + -0.212737 *lens_ipow(begin_dy, 6) + -0.084192 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 4) + -0.0896338 *lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 2) + 2.33254e-07 *lens_ipow(begin_y, 4)*lens_ipow(begin_dy, 2) + -2.68231e-05 *lens_ipow(begin_x, 2)*lens_ipow(begin_lambda, 4) + -1.32759e-10 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 4) + 9.16575e-08 *lens_ipow(begin_x, 4)*lens_ipow(begin_dy, 2) + 3.28415e-07 *lens_ipow(begin_x, 4)*lens_ipow(begin_dx, 2) + -2.33314e-10 *lens_ipow(begin_x, 6) + 0.00218107 *begin_y*begin_dy*lens_ipow(begin_lambda, 5) + 2.26227e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 2) + -0.0648339 *lens_ipow(begin_dx, 8) + -6.18196e-13 *lens_ipow(begin_y, 8) + 2.21726e-07 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 4) + -0.243532 *lens_ipow(begin_dx, 4)*lens_ipow(begin_lambda, 5) + -1.68167e-05 *lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 7) + 1.13661e-06 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 2) + -1.12536e-15 *lens_ipow(begin_x, 8)*lens_ipow(begin_y, 2) + -0.121461 *lens_ipow(begin_lambda, 11);
else
  out[4] = 0.0f;
