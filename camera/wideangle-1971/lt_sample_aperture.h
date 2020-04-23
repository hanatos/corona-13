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
       + 4.29068e-05  + 131.941 *begin_dx + 1.47919e-06 *begin_y + 0.760868 *begin_x + -0.00248268 *lens_ipow(begin_dy, 2) + 3.05647e-05 *begin_x*begin_dy + 6.79418e-05 *begin_x*begin_dx + -5.21674e-07 *begin_x*begin_y + 4.33109e-07 *lens_ipow(begin_x, 2) + 31.7351 *begin_dx*lens_ipow(begin_dy, 2) + 32.225 *lens_ipow(begin_dx, 3) + 1.28139 *begin_y*begin_dx*begin_dy + 0.00677532 *lens_ipow(begin_y, 2)*begin_dx + 0.721853 *begin_x*lens_ipow(begin_dy, 2) + 2.03673 *begin_x*lens_ipow(begin_dx, 2) + 0.0146307 *begin_x*begin_y*begin_dy + 6.14833e-05 *begin_x*lens_ipow(begin_y, 2) + 0.0215743 *lens_ipow(begin_x, 2)*begin_dx + 6.19039e-05 *lens_ipow(begin_x, 3) + -0.0115372 *begin_x*lens_ipow(begin_lambda, 3) + -4.70837e-07 *begin_x*lens_ipow(begin_y, 2)*begin_dy + 0.12077 *begin_dx*lens_ipow(begin_lambda, 4) + 4.08016e-05 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 2) + -5.33215e-10 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 2) + -0.359659 *begin_x*lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 4) + 9.11649e-10 *begin_x*lens_ipow(begin_y, 5)*begin_dy + 0.20722 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 2) + 0.10681 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 5) + 0.000134914 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dx, 3)*begin_dy + 6.65254e-08 *lens_ipow(begin_x, 5)*lens_ipow(begin_dy, 2) + 1.41142e-07 *lens_ipow(begin_x, 5)*lens_ipow(begin_dx, 2) + -0.00128574 *begin_x*begin_y*begin_dy*lens_ipow(begin_lambda, 5) + -0.00393656 *lens_ipow(begin_x, 2)*begin_dx*lens_ipow(begin_lambda, 5) + -555.619 *lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 4) + -0.783692 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 7) + -462.046 *lens_ipow(begin_dx, 5)*lens_ipow(begin_lambda, 5) + -2.75031 *begin_x*lens_ipow(begin_dy, 4)*lens_ipow(begin_lambda, 5) + -3.30349e-08 *lens_ipow(begin_x, 5)*lens_ipow(begin_lambda, 5) + -10.5193 *begin_y*lens_ipow(begin_dx, 3)*begin_dy*lens_ipow(begin_lambda, 6) + 0.0659771 *begin_x*lens_ipow(begin_lambda, 10),
       + 4.77201e-05  + 131.943 *begin_dy + 0.760961 *begin_y + 4.46534e-07 *begin_x + 2.72166e-05 *begin_y*begin_dx + -6.26592e-07 *lens_ipow(begin_y, 2) + 3.04161e-05 *begin_x*begin_dx + 4.35454e-07 *begin_x*begin_y + 32.3847 *lens_ipow(begin_dy, 3) + 32.0399 *lens_ipow(begin_dx, 2)*begin_dy + 2.03225 *begin_y*lens_ipow(begin_dy, 2) + 0.725552 *begin_y*lens_ipow(begin_dx, 2) + 0.021516 *lens_ipow(begin_y, 2)*begin_dy + 6.14628e-05 *lens_ipow(begin_y, 3) + 1.29155 *begin_x*begin_dx*begin_dy + 0.0147175 *begin_x*begin_y*begin_dx + 0.00679404 *lens_ipow(begin_x, 2)*begin_dy + 6.12551e-05 *lens_ipow(begin_x, 2)*begin_y + -0.0119604 *begin_y*lens_ipow(begin_lambda, 3) + 3.7335e-05 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2) + 0.0743001 *begin_dy*lens_ipow(begin_lambda, 4) + 7.39082e-05 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dx, 2) + 2.34632e-09 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 3) + -0.30367 *begin_y*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 4) + 0.0510535 *lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 5) + -8.41132 *begin_x*lens_ipow(begin_dx, 5)*begin_dy + -0.00171867 *begin_x*lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 3) + 1.41059e-09 *lens_ipow(begin_x, 5)*begin_y*begin_dx + -0.00326215 *lens_ipow(begin_y, 2)*begin_dy*lens_ipow(begin_lambda, 5) + -0.00247691 *begin_x*begin_y*begin_dx*lens_ipow(begin_lambda, 5) + -7.73829e-05 *lens_ipow(begin_x, 3)*begin_dx*begin_dy*lens_ipow(begin_lambda, 3) + -265.517 *lens_ipow(begin_dy, 5)*lens_ipow(begin_lambda, 4) + -3.7727 *begin_y*lens_ipow(begin_dx, 4)*lens_ipow(begin_lambda, 4) + 4.10599e-06 *lens_ipow(begin_y, 5)*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + -557.907 *lens_ipow(begin_dx, 4)*begin_dy*lens_ipow(begin_lambda, 5) + -4.06703e-05 *lens_ipow(begin_x, 4)*lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 3) + 0.0658052 *begin_y*lens_ipow(begin_lambda, 10) + 1055.34 *begin_y*lens_ipow(begin_dx, 10) + -0.512701 *begin_x*begin_dx*begin_dy*lens_ipow(begin_lambda, 8) + -5.85707e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 3)*lens_ipow(begin_lambda, 6)
    };
    const float delta_ap[] = {ap_x - pred_ap[0], ap_y - pred_ap[1]};
    sqr_ap_err = delta_ap[0]*delta_ap[0]+delta_ap[1]*delta_ap[1];
    float dx1_domega0[2][2];
    dx1_domega0[0][0] =  + 131.941  + 6.79418e-05 *begin_x + 31.7351 *lens_ipow(begin_dy, 2) + 96.675 *lens_ipow(begin_dx, 2) + 1.28139 *begin_y*begin_dy + 0.00677532 *lens_ipow(begin_y, 2) + 4.07347 *begin_x*begin_dx + 0.0215743 *lens_ipow(begin_x, 2) + 0.12077 *lens_ipow(begin_lambda, 4) + -0.719317 *begin_x*begin_dx*lens_ipow(begin_lambda, 4) + 0.62166 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + 0.534052 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 4) + 0.000404743 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dx, 2)*begin_dy + 2.82284e-07 *lens_ipow(begin_x, 5)*begin_dx + -0.00393656 *lens_ipow(begin_x, 2)*lens_ipow(begin_lambda, 5) + -1666.86 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 4) + -5.48584 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 6) + -2310.23 *lens_ipow(begin_dx, 4)*lens_ipow(begin_lambda, 5) + -31.5578 *begin_y*lens_ipow(begin_dx, 2)*begin_dy*lens_ipow(begin_lambda, 6)+0.0f;
    dx1_domega0[0][1] =  + -0.00496536 *begin_dy + 3.05647e-05 *begin_x + 63.4702 *begin_dx*begin_dy + 1.28139 *begin_y*begin_dx + 1.44371 *begin_x*begin_dy + 0.0146307 *begin_x*begin_y + -4.70837e-07 *begin_x*lens_ipow(begin_y, 2) + 8.16031e-05 *begin_x*lens_ipow(begin_y, 2)*begin_dy + 9.11649e-10 *begin_x*lens_ipow(begin_y, 5) + 0.41444 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 3)*begin_dy + 0.000134914 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dx, 3) + 1.33051e-07 *lens_ipow(begin_x, 5)*begin_dy + -0.00128574 *begin_x*begin_y*lens_ipow(begin_lambda, 5) + -1111.24 *lens_ipow(begin_dx, 3)*begin_dy*lens_ipow(begin_lambda, 4) + -11.0012 *begin_x*lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 5) + -10.5193 *begin_y*lens_ipow(begin_dx, 3)*lens_ipow(begin_lambda, 6)+0.0f;
    dx1_domega0[1][0] =  + 2.72166e-05 *begin_y + 3.04161e-05 *begin_x + 64.0799 *begin_dx*begin_dy + 1.4511 *begin_y*begin_dx + 1.29155 *begin_x*begin_dy + 0.0147175 *begin_x*begin_y + 7.46701e-05 *lens_ipow(begin_y, 2)*begin_dx + 0.000147816 *lens_ipow(begin_x, 2)*begin_y*begin_dx + -42.0566 *begin_x*lens_ipow(begin_dx, 4)*begin_dy + -0.00171867 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 3) + 1.41059e-09 *lens_ipow(begin_x, 5)*begin_y + -0.00247691 *begin_x*begin_y*lens_ipow(begin_lambda, 5) + -7.73829e-05 *lens_ipow(begin_x, 3)*begin_dy*lens_ipow(begin_lambda, 3) + -15.0908 *begin_y*lens_ipow(begin_dx, 3)*lens_ipow(begin_lambda, 4) + 8.21199e-06 *lens_ipow(begin_y, 5)*begin_dx*lens_ipow(begin_dy, 2) + -2231.63 *lens_ipow(begin_dx, 3)*begin_dy*lens_ipow(begin_lambda, 5) + 10553.4 *begin_y*lens_ipow(begin_dx, 9) + -0.512701 *begin_x*begin_dy*lens_ipow(begin_lambda, 8)+0.0f;
    dx1_domega0[1][1] =  + 131.943  + 97.1542 *lens_ipow(begin_dy, 2) + 32.0399 *lens_ipow(begin_dx, 2) + 4.06449 *begin_y*begin_dy + 0.021516 *lens_ipow(begin_y, 2) + 1.29155 *begin_x*begin_dx + 0.00679404 *lens_ipow(begin_x, 2) + 0.0743001 *lens_ipow(begin_lambda, 4) + -0.607341 *begin_y*begin_dy*lens_ipow(begin_lambda, 4) + 0.255267 *lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 4) + -8.41132 *begin_x*lens_ipow(begin_dx, 5) + -0.00515601 *begin_x*lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 2) + -0.00326215 *lens_ipow(begin_y, 2)*lens_ipow(begin_lambda, 5) + -7.73829e-05 *lens_ipow(begin_x, 3)*begin_dx*lens_ipow(begin_lambda, 3) + -1327.58 *lens_ipow(begin_dy, 4)*lens_ipow(begin_lambda, 4) + 8.21199e-06 *lens_ipow(begin_y, 5)*lens_ipow(begin_dx, 2)*begin_dy + -557.907 *lens_ipow(begin_dx, 4)*lens_ipow(begin_lambda, 5) + -0.000122011 *lens_ipow(begin_x, 4)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 3) + -0.512701 *begin_x*begin_dx*lens_ipow(begin_lambda, 8)+0.0f;
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
    out[0] =  + 0.00020339  + 97.9694 *begin_dx + -5.56572e-07 *begin_y + -0.15792 *begin_x + 0.000131035 *begin_y*begin_dy + 0.000216653 *begin_x*begin_dx + -2.83842e-06 *begin_x*begin_y + -40.4004 *begin_dx*lens_ipow(begin_dy, 2) + -39.7916 *lens_ipow(begin_dx, 3) + 0.291008 *begin_y*begin_dx*begin_dy + 0.00516604 *lens_ipow(begin_y, 2)*begin_dx + 0.555742 *begin_x*lens_ipow(begin_dy, 2) + 0.825372 *begin_x*lens_ipow(begin_dx, 2) + 0.00822883 *begin_x*begin_y*begin_dy + -2.0305e-05 *begin_x*begin_y*begin_dx + 5.87607e-05 *begin_x*lens_ipow(begin_y, 2) + 0.00976077 *lens_ipow(begin_x, 2)*begin_dx + 5.80578e-05 *lens_ipow(begin_x, 3) + 8.80202 *begin_dx*lens_ipow(begin_lambda, 3) + -0.00638787 *begin_x*lens_ipow(begin_lambda, 3) + -1.72929e-06 *begin_x*lens_ipow(begin_y, 2)*begin_dy + 0.00727222 *lens_ipow(begin_x, 2)*begin_dx*begin_lambda + -0.00173765 *lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_lambda, 2) + -7.36847e-06 *lens_ipow(begin_x, 2)*begin_y*begin_dx*begin_dy + 1.43833e-10 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 2) + -0.898643 *begin_y*begin_dx*begin_dy*lens_ipow(begin_lambda, 3) + 0.105484 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 3)*begin_lambda + 35.2561 *begin_x*lens_ipow(begin_dx, 6) + 0.0537127 *begin_x*begin_y*lens_ipow(begin_dy, 5) + -0.310332 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 2) + 0.000105879 *lens_ipow(begin_x, 4)*begin_dx*lens_ipow(begin_dy, 2) + 8.02862e-07 *lens_ipow(begin_x, 5)*lens_ipow(begin_dy, 2) + 8.088e-07 *lens_ipow(begin_x, 5)*lens_ipow(begin_dx, 2) + -106.691 *begin_dx*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 5) + -195.674 *lens_ipow(begin_dx, 3)*lens_ipow(begin_lambda, 6) + -0.0101361 *begin_x*begin_y*begin_dy*lens_ipow(begin_lambda, 6) + -0.071513 *lens_ipow(begin_x, 2)*begin_dx*lens_ipow(begin_lambda, 7) + -34.997 *begin_dx*lens_ipow(begin_lambda, 10) + -7.11532 *begin_x*lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 8) + -0.0036996 *lens_ipow(begin_x, 3)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 6);
    out[1] =  + 0.000127718  + 98.0203 *begin_dy + -0.00168766 *begin_dx + -0.157835 *begin_y + -1.27501e-05 *begin_x + -0.0221099 *begin_dx*begin_dy + 9.09981e-05 *begin_y*begin_dx + -0.000177701 *begin_x*begin_dy + -42.5469 *lens_ipow(begin_dy, 3) + -45.2442 *lens_ipow(begin_dx, 2)*begin_dy + 0.874717 *begin_y*lens_ipow(begin_dy, 2) + -0.00332038 *begin_y*begin_dx*begin_dy + 0.55939 *begin_y*lens_ipow(begin_dx, 2) + 0.012835 *lens_ipow(begin_y, 2)*begin_dy + 5.95874e-05 *lens_ipow(begin_y, 3) + 0.216174 *begin_x*begin_dx*begin_dy + -2.88804e-05 *begin_x*begin_y*begin_dy + 0.00812763 *begin_x*begin_y*begin_dx + 0.00505774 *lens_ipow(begin_x, 2)*begin_dy + 5.94187e-05 *lens_ipow(begin_x, 2)*begin_y + 8.76415 *begin_dy*lens_ipow(begin_lambda, 3) + -0.00728182 *begin_y*lens_ipow(begin_lambda, 3) + 0.00624187 *begin_y*lens_ipow(begin_dy, 3) + 0.000459158 *lens_ipow(begin_y, 2)*begin_dy*begin_lambda + 0.0116774 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*begin_dy + -0.00101859 *lens_ipow(begin_x, 2)*begin_dy*lens_ipow(begin_lambda, 2) + -1.03325e-06 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*begin_dy + -8.39618e-09 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 3) + -0.943858 *begin_y*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 3) + -0.390514 *begin_x*begin_dx*begin_dy*lens_ipow(begin_lambda, 3) + -5.30612 *begin_y*lens_ipow(begin_dx, 6) + 0.221248 *lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 5) + -2.91859e-09 *begin_x*lens_ipow(begin_y, 5)*begin_dx + -0.00573133 *begin_x*begin_y*begin_dx*lens_ipow(begin_lambda, 6) + 0.155317 *begin_x*begin_y*begin_dx*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 4) + -32.9113 *begin_dy*lens_ipow(begin_lambda, 9) + -1352.65 *lens_ipow(begin_dy, 5)*lens_ipow(begin_lambda, 5) + -1782.48 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 5) + -0.033834 *lens_ipow(begin_y, 2)*begin_dy*lens_ipow(begin_lambda, 8) + -6.96656e-06 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*begin_dy*lens_ipow(begin_lambda, 6);
    out[2] =  + -1.04483 *begin_dx + -0.00852243 *begin_x + -0.103975 *begin_dx*lens_ipow(begin_dy, 2) + 0.46888 *lens_ipow(begin_dx, 3) + 0.00125944 *begin_y*begin_dx*begin_dy + 6.29858e-06 *lens_ipow(begin_y, 2)*begin_dx + -0.00444794 *begin_x*lens_ipow(begin_dy, 2) + -7.06817e-05 *begin_x*begin_y*begin_dy + -2.23936e-07 *begin_x*lens_ipow(begin_y, 2) + -7.42772e-05 *lens_ipow(begin_x, 2)*begin_dx + -2.1259e-07 *lens_ipow(begin_x, 3) + 0.00082206 *begin_x*lens_ipow(begin_lambda, 3) + -0.004355 *begin_x*lens_ipow(begin_lambda, 10);
    out[3] =  + -1.04491 *begin_dy + -0.00852134 *begin_y + 0.466359 *lens_ipow(begin_dy, 3) + 1.01739 *lens_ipow(begin_dx, 2)*begin_dy + -0.0060394 *begin_y*lens_ipow(begin_dx, 2) + -7.35944e-05 *lens_ipow(begin_y, 2)*begin_dy + -2.10827e-07 *lens_ipow(begin_y, 3) + 0.00879264 *begin_x*begin_dx*begin_dy + -8.10567e-05 *begin_x*begin_y*begin_dx + -6.06527e-06 *lens_ipow(begin_x, 2)*begin_dy + -2.08789e-07 *lens_ipow(begin_x, 2)*begin_y + 0.000812826 *begin_y*lens_ipow(begin_lambda, 3) + -0.00428298 *begin_y*lens_ipow(begin_lambda, 10);
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
    domega2_dx0[0][0] =  + -0.00852243  + -0.00444794 *lens_ipow(begin_dy, 2) + -7.06817e-05 *begin_y*begin_dy + -2.23936e-07 *lens_ipow(begin_y, 2) + -0.000148554 *begin_x*begin_dx + -6.37771e-07 *lens_ipow(begin_x, 2) + 0.00082206 *lens_ipow(begin_lambda, 3) + -0.004355 *lens_ipow(begin_lambda, 10)+0.0f;
    domega2_dx0[0][1] =  + 0.00125944 *begin_dx*begin_dy + 1.25972e-05 *begin_y*begin_dx + -7.06817e-05 *begin_x*begin_dy + -4.47871e-07 *begin_x*begin_y+0.0f;
    domega2_dx0[1][0] =  + 0.00879264 *begin_dx*begin_dy + -8.10567e-05 *begin_y*begin_dx + -1.21305e-05 *begin_x*begin_dy + -4.17578e-07 *begin_x*begin_y+0.0f;
    domega2_dx0[1][1] =  + -0.00852134  + -0.0060394 *lens_ipow(begin_dx, 2) + -0.000147189 *begin_y*begin_dy + -6.3248e-07 *lens_ipow(begin_y, 2) + -8.10567e-05 *begin_x*begin_dx + -2.08789e-07 *lens_ipow(begin_x, 2) + 0.000812826 *lens_ipow(begin_lambda, 3) + -0.00428298 *lens_ipow(begin_lambda, 10)+0.0f;
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
  out[4] =  + 0.250918  + 0.300065 *begin_lambda + -3.43167e-05 *begin_dx + 5.94618e-08 *begin_y + -1.76384e-07 *begin_x + -0.0491003 *lens_ipow(begin_dy, 2) + -0.0501791 *lens_ipow(begin_dx, 2) + -0.000258008 *begin_y*begin_dy + -1.41015e-06 *lens_ipow(begin_y, 2) + -0.000202527 *begin_x*begin_dx + -7.67748e-07 *lens_ipow(begin_x, 2) + -0.247432 *lens_ipow(begin_lambda, 3) + -0.0273995 *begin_y*lens_ipow(begin_dy, 3) + -0.0135167 *begin_y*lens_ipow(begin_dx, 2)*begin_dy + -0.000264022 *lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 2) + -9.7489e-07 *lens_ipow(begin_y, 3)*begin_dy + -0.0256067 *begin_x*begin_dx*lens_ipow(begin_dy, 2) + -0.0326559 *begin_x*lens_ipow(begin_dx, 3) + -0.000397084 *begin_x*begin_y*begin_dx*begin_dy + -1.02748e-06 *begin_x*lens_ipow(begin_y, 2)*begin_dx + -7.4779e-05 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 2) + -0.000311522 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 2) + -1.11454e-06 *lens_ipow(begin_x, 2)*begin_y*begin_dy + -9.93961e-07 *lens_ipow(begin_x, 3)*begin_dx + -26.2272 *lens_ipow(begin_dy, 6) + -64.7721 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 4) + -92.9345 *lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 2) + -28.0047 *lens_ipow(begin_dx, 6) + -0.767965 *begin_y*lens_ipow(begin_dx, 4)*begin_dy + -0.0030545 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 4) + -1.58875e-05 *lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 3) + 0.00532933 *begin_x*begin_y*begin_dx*lens_ipow(begin_dy, 3) + -6.12124e-08 *begin_x*lens_ipow(begin_y, 3)*begin_dx*begin_dy + 6.46381e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2) + -1.69391e-11 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 4) + -1.62917e-05 *lens_ipow(begin_x, 3)*begin_dx*lens_ipow(begin_dy, 2) + -2.05514e-05 *lens_ipow(begin_x, 3)*lens_ipow(begin_dx, 3) + -6.76989e-12 *lens_ipow(begin_x, 6) + -2.60409e-06 *lens_ipow(begin_x, 4)*lens_ipow(begin_dy, 4) + 0.438506 *lens_ipow(begin_lambda, 11);
else
  out[4] = 0.0f;
