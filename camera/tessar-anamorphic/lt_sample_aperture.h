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
       + -7.04901e-05  + 88.7549 *begin_dx + -3.59444e-06 *begin_y + 0.746116 *begin_x + 0.00393364 *lens_ipow(begin_dy, 2) + 0.0027943 *lens_ipow(begin_dx, 2) + 5.7837e-05 *begin_x*begin_dy + 61.8435 *begin_dx*lens_ipow(begin_dy, 2) + 62.4798 *lens_ipow(begin_dx, 3) + 1.79321 *begin_y*begin_dx*begin_dy + 0.00947158 *lens_ipow(begin_y, 2)*begin_dx + 1.17879 *begin_x*lens_ipow(begin_dy, 2) + 3.00255 *begin_x*lens_ipow(begin_dx, 2) + 0.025828 *begin_x*begin_y*begin_dy + 0.000118713 *begin_x*lens_ipow(begin_y, 2) + 0.0365142 *lens_ipow(begin_x, 2)*begin_dx + 0.000119619 *lens_ipow(begin_x, 3) + 3.65264 *begin_dx*lens_ipow(begin_lambda, 3) + 0.00209122 *lens_ipow(begin_y, 2)*begin_dx*begin_lambda + 0.0419482 *begin_x*lens_ipow(begin_lambda, 3) + 5.82037e-06 *lens_ipow(begin_x, 3)*lens_ipow(begin_lambda, 2) + 0.000256645 *lens_ipow(begin_x, 3)*lens_ipow(begin_dy, 2)*begin_lambda + 0.0417051 *begin_y*begin_dx*begin_dy*lens_ipow(begin_lambda, 4) + 1.74644e-05 *lens_ipow(begin_y, 4)*begin_dx*lens_ipow(begin_dy, 2) + 1.64312e-05 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_lambda, 4) + 2.23226e-07 *begin_x*lens_ipow(begin_y, 4)*lens_ipow(begin_dy, 2) + 0.147654 *lens_ipow(begin_x, 2)*begin_dx*lens_ipow(begin_dy, 4) + -0.00626952 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dx, 3)*begin_dy + -1.38128e-08 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 2)*begin_dx + -9.19947e-11 *lens_ipow(begin_x, 5)*lens_ipow(begin_y, 2) + -0.00316643 *lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_lambda, 5) + -2.10946 *begin_x*lens_ipow(begin_dy, 6)*begin_lambda + 0.033881 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 3)*lens_ipow(begin_lambda, 3) + -12091.8 *lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 6) + -9173.74 *lens_ipow(begin_dx, 7)*lens_ipow(begin_dy, 2) + -19.1573 *begin_dx*lens_ipow(begin_lambda, 10) + -1567.04 *lens_ipow(begin_dx, 7)*lens_ipow(begin_lambda, 4) + -0.604442 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 7)*lens_ipow(begin_lambda, 2) + -0.226439 *begin_x*lens_ipow(begin_lambda, 10) + 0.00770434 *begin_x*begin_y*begin_dy*lens_ipow(begin_lambda, 8),
       + 0.000111506  + 88.7444 *begin_dy + 0.747015 *begin_y + -9.31942e-06 *begin_x + 0.00496451 *begin_dx*begin_dy + -7.22249e-07 *begin_x*begin_y + 61.5721 *lens_ipow(begin_dy, 3) + 61.457 *lens_ipow(begin_dx, 2)*begin_dy + 3.00018 *begin_y*lens_ipow(begin_dy, 2) + 1.18915 *begin_y*lens_ipow(begin_dx, 2) + 0.0367408 *lens_ipow(begin_y, 2)*begin_dy + 0.000121283 *lens_ipow(begin_y, 3) + 1.79209 *begin_x*begin_dx*begin_dy + 1.52149e-05 *begin_x*begin_y*begin_dy + 0.0259214 *begin_x*begin_y*begin_dx + 2.28374e-07 *begin_x*lens_ipow(begin_y, 2) + 0.0104326 *lens_ipow(begin_x, 2)*begin_dy + 0.000121067 *lens_ipow(begin_x, 2)*begin_y + 3.8145 *begin_dy*lens_ipow(begin_lambda, 3) + 0.000296887 *begin_x*begin_dx*lens_ipow(begin_lambda, 2) + 4.52989e-05 *begin_x*begin_y*begin_dx*begin_dy + 0.068935 *begin_y*lens_ipow(begin_lambda, 4) + 0.0224844 *begin_x*begin_y*begin_dx*lens_ipow(begin_dy, 2)*begin_lambda + 0.184803 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 3) + -0.00264997 *lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 4) + 2.77742e-07 *lens_ipow(begin_y, 5)*lens_ipow(begin_dx, 2) + -9.07379e-09 *lens_ipow(begin_y, 6)*begin_dy + 5.60093e-07 *lens_ipow(begin_x, 4)*begin_y*lens_ipow(begin_dy, 2) + -136.766 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 4) + 0.000789762 *lens_ipow(begin_x, 4)*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 3) + 8.43138e-07 *lens_ipow(begin_x, 4)*begin_y*lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 2) + -1.00385e-13 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 5) + -7.47471e-14 *lens_ipow(begin_x, 8)*begin_y + -15.2218 *begin_y*lens_ipow(begin_dx, 6)*lens_ipow(begin_lambda, 3) + -20.1841 *begin_dy*lens_ipow(begin_lambda, 10) + -28862.2 *lens_ipow(begin_dy, 11) + -89123.6 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 9) + -245931 *lens_ipow(begin_dx, 6)*lens_ipow(begin_dy, 5) + -0.335914 *begin_y*lens_ipow(begin_lambda, 10) + -4.55116e-16 *lens_ipow(begin_y, 11)
    };
    const float delta_ap[] = {ap_x - pred_ap[0], ap_y - pred_ap[1]};
    sqr_ap_err = delta_ap[0]*delta_ap[0]+delta_ap[1]*delta_ap[1];
    float dx1_domega0[2][2];
    dx1_domega0[0][0] =  + 88.7549  + 0.00558861 *begin_dx + 61.8435 *lens_ipow(begin_dy, 2) + 187.439 *lens_ipow(begin_dx, 2) + 1.79321 *begin_y*begin_dy + 0.00947158 *lens_ipow(begin_y, 2) + 6.00511 *begin_x*begin_dx + 0.0365142 *lens_ipow(begin_x, 2) + 3.65264 *lens_ipow(begin_lambda, 3) + 0.00209122 *lens_ipow(begin_y, 2)*begin_lambda + 0.0417051 *begin_y*begin_dy*lens_ipow(begin_lambda, 4) + 1.74644e-05 *lens_ipow(begin_y, 4)*lens_ipow(begin_dy, 2) + 0.147654 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 4) + -0.0188085 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dx, 2)*begin_dy + -1.38128e-08 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 2) + -0.00316643 *lens_ipow(begin_y, 2)*lens_ipow(begin_lambda, 5) + 0.101643 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 3) + -36275.5 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 6) + -64216.2 *lens_ipow(begin_dx, 6)*lens_ipow(begin_dy, 2) + -19.1573 *lens_ipow(begin_lambda, 10) + -10969.3 *lens_ipow(begin_dx, 6)*lens_ipow(begin_lambda, 4) + -4.2311 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 6)*lens_ipow(begin_lambda, 2)+0.0f;
    dx1_domega0[0][1] =  + 0.00786728 *begin_dy + 5.7837e-05 *begin_x + 123.687 *begin_dx*begin_dy + 1.79321 *begin_y*begin_dx + 2.35758 *begin_x*begin_dy + 0.025828 *begin_x*begin_y + 0.000513289 *lens_ipow(begin_x, 3)*begin_dy*begin_lambda + 0.0417051 *begin_y*begin_dx*lens_ipow(begin_lambda, 4) + 3.49289e-05 *lens_ipow(begin_y, 4)*begin_dx*begin_dy + 4.46453e-07 *begin_x*lens_ipow(begin_y, 4)*begin_dy + 0.590615 *lens_ipow(begin_x, 2)*begin_dx*lens_ipow(begin_dy, 3) + -0.00626952 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dx, 3) + -12.6568 *begin_x*lens_ipow(begin_dy, 5)*begin_lambda + -72551 *lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 5) + -18347.5 *lens_ipow(begin_dx, 7)*begin_dy + 0.00770434 *begin_x*begin_y*lens_ipow(begin_lambda, 8)+0.0f;
    dx1_domega0[1][0] =  + 0.00496451 *begin_dy + 122.914 *begin_dx*begin_dy + 2.3783 *begin_y*begin_dx + 1.79209 *begin_x*begin_dy + 0.0259214 *begin_x*begin_y + 0.000296887 *begin_x*lens_ipow(begin_lambda, 2) + 4.52989e-05 *begin_x*begin_y*begin_dy + 0.0224844 *begin_x*begin_y*lens_ipow(begin_dy, 2)*begin_lambda + 0.369606 *lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 3) + 5.55484e-07 *lens_ipow(begin_y, 5)*begin_dx + -273.532 *begin_dx*lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 4) + 0.00157952 *lens_ipow(begin_x, 4)*begin_dx*lens_ipow(begin_dy, 3) + 1.68628e-06 *lens_ipow(begin_x, 4)*begin_y*begin_dx*lens_ipow(begin_lambda, 2) + -91.3306 *begin_y*lens_ipow(begin_dx, 5)*lens_ipow(begin_lambda, 3) + -178247 *begin_dx*lens_ipow(begin_dy, 9) + -1.47558e+06 *lens_ipow(begin_dx, 5)*lens_ipow(begin_dy, 5)+0.0f;
    dx1_domega0[1][1] =  + 88.7444  + 0.00496451 *begin_dx + 184.716 *lens_ipow(begin_dy, 2) + 61.457 *lens_ipow(begin_dx, 2) + 6.00036 *begin_y*begin_dy + 0.0367408 *lens_ipow(begin_y, 2) + 1.79209 *begin_x*begin_dx + 1.52149e-05 *begin_x*begin_y + 0.0104326 *lens_ipow(begin_x, 2) + 3.8145 *lens_ipow(begin_lambda, 3) + 4.52989e-05 *begin_x*begin_y*begin_dx + 0.0449689 *begin_x*begin_y*begin_dx*begin_dy*begin_lambda + 0.554409 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + -0.0105999 *lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 3) + -9.07379e-09 *lens_ipow(begin_y, 6) + 1.12019e-06 *lens_ipow(begin_x, 4)*begin_y*begin_dy + -410.297 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 4) + 0.00236929 *lens_ipow(begin_x, 4)*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + -20.1841 *lens_ipow(begin_lambda, 10) + -317485 *lens_ipow(begin_dy, 10) + -802113 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 8) + -1.22965e+06 *lens_ipow(begin_dx, 6)*lens_ipow(begin_dy, 4)+0.0f;
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
    out[0] =  + -0.000251644  + -0.00169379 *begin_dy + 138.55 *begin_dx + -3.17676e-05 *begin_y + 0.520192 *begin_x + 0.000180213 *begin_y*begin_dx + 3.07875e-06 *lens_ipow(begin_y, 2) + -2.79049e-06 *begin_x*begin_y + -15.2536 *begin_dx*lens_ipow(begin_dy, 2) + -23.5386 *lens_ipow(begin_dx, 3) + 1.42507 *begin_y*begin_dx*begin_dy + 0.0154736 *lens_ipow(begin_y, 2)*begin_dx + 0.547231 *begin_x*lens_ipow(begin_dy, 2) + 2.5213 *begin_x*lens_ipow(begin_dx, 2) + 0.0262714 *begin_x*begin_y*begin_dy + -8.58505e-06 *begin_x*begin_y*begin_dx + 0.000188386 *begin_x*lens_ipow(begin_y, 2) + 0.050176 *lens_ipow(begin_x, 2)*begin_dx + 0.000218817 *lens_ipow(begin_x, 3) + -8.46306 *begin_dx*lens_ipow(begin_lambda, 3) + 0.000190627 *lens_ipow(begin_y, 2)*begin_dx*begin_dy + 3.72883e-06 *begin_x*lens_ipow(begin_y, 2)*begin_dx + 0.000154273 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 2) + -0.0922462 *begin_x*lens_ipow(begin_lambda, 4) + 0.742721 *lens_ipow(begin_dx, 4)*lens_ipow(begin_lambda, 2) + 0.0602294 *begin_y*lens_ipow(begin_dy, 5) + 1881.05 *lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 4) + 0.0034023 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + 0.0021023 *lens_ipow(begin_x, 3)*lens_ipow(begin_dy, 4) + 3.16763e-07 *lens_ipow(begin_x, 5)*lens_ipow(begin_dx, 2) + -9.42875e-11 *lens_ipow(begin_x, 5)*lens_ipow(begin_y, 2) + -0.0807065 *begin_x*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 5) + 2050.32 *begin_dx*lens_ipow(begin_dy, 8) + 11495.3 *lens_ipow(begin_dx, 7)*lens_ipow(begin_dy, 2) + 5062.85 *lens_ipow(begin_dx, 9) + 3.83813e-12 *lens_ipow(begin_y, 8)*begin_dx + 8.44829e-08 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 2)*lens_ipow(begin_lambda, 4) + 5.27165e-06 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 4) + 44.979 *begin_dx*lens_ipow(begin_lambda, 10) + 0.442883 *begin_x*lens_ipow(begin_lambda, 10);
    out[1] =  + 0.000727971  + 102.059 *begin_dy + -0.0022679 *begin_dx + 0.283204 *begin_y + -2.81477e-05 *begin_x + 0.0348028 *begin_dx*begin_dy + 0.000360639 *begin_y*begin_dx + -2.24295e-06 *lens_ipow(begin_y, 2) + 0.000610253 *begin_x*begin_dy + -1.40208e-06 *lens_ipow(begin_x, 2) + -0.00187788 *lens_ipow(begin_lambda, 3) + -4.61012 *lens_ipow(begin_dy, 3) + 0.0677415 *begin_dx*lens_ipow(begin_dy, 2) + 0.486623 *lens_ipow(begin_dx, 2)*begin_dy + 1.64922 *begin_y*lens_ipow(begin_dy, 2) + 1.06526 *begin_y*lens_ipow(begin_dx, 2) + 0.0316329 *lens_ipow(begin_y, 2)*begin_dy + 0.000139082 *lens_ipow(begin_y, 3) + 1.24182 *begin_x*begin_dx*begin_dy + 0.0286956 *begin_x*begin_y*begin_dx + 1.51945e-07 *begin_x*lens_ipow(begin_y, 2) + 0.0110697 *lens_ipow(begin_x, 2)*begin_dy + 0.000168025 *lens_ipow(begin_x, 2)*begin_y + 1.02608 *begin_dy*lens_ipow(begin_lambda, 3) + 0.000336615 *begin_x*begin_y*begin_dx*begin_dy + 3.12224e-08 *begin_x*lens_ipow(begin_y, 3) + -8.89702e-07 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*begin_dy + 9.27523e-05 *lens_ipow(begin_x, 3)*begin_dx*lens_ipow(begin_dy, 2) + 787.263 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 5) + 639.45 *lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 3) + 347.714 *lens_ipow(begin_dx, 6)*begin_dy + -0.162013 *begin_y*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 4) + 5.03266e-07 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 2) + 5.38559e-10 *lens_ipow(begin_y, 7)*lens_ipow(begin_dy, 2) + -2.69159e-13 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 5) + -4.88984e-08 *lens_ipow(begin_x, 5)*begin_y*begin_dx*lens_ipow(begin_dy, 2) + -6.52236 *begin_dy*lens_ipow(begin_lambda, 10) + -1499.79 *lens_ipow(begin_dy, 7)*lens_ipow(begin_lambda, 4) + 42482.7 *lens_ipow(begin_dy, 11) + -0.000318588 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 5)*lens_ipow(begin_lambda, 2);
    out[2] =  + -1.14332 *begin_dx + -0.0115217 *begin_x + 0.623924 *begin_dx*lens_ipow(begin_dy, 2) + 0.675525 *lens_ipow(begin_dx, 3) + -2.85112e-05 *lens_ipow(begin_y, 2)*begin_dx + 0.000652111 *begin_x*lens_ipow(begin_dy, 2) + -0.00350055 *begin_x*lens_ipow(begin_dx, 2) + -8.36284e-05 *begin_x*begin_y*begin_dy + -2.83412e-07 *begin_x*lens_ipow(begin_y, 2) + -0.000189131 *lens_ipow(begin_x, 2)*begin_dx + -6.03404e-07 *lens_ipow(begin_x, 3) + 0.0680919 *begin_dx*lens_ipow(begin_lambda, 3) + -1.72969 *lens_ipow(begin_dx, 7) + -0.369181 *begin_dx*lens_ipow(begin_lambda, 10);
    out[3] =  + -0.908274 *begin_dy + -0.0123092 *begin_y + 0.614547 *lens_ipow(begin_dy, 3) + 1.05366 *lens_ipow(begin_dx, 2)*begin_dy + 0.00443453 *begin_y*lens_ipow(begin_dy, 2) + -0.00136161 *begin_y*lens_ipow(begin_dx, 2) + -2.68644e-05 *lens_ipow(begin_y, 2)*begin_dy + 8.06763e-08 *lens_ipow(begin_y, 3) + 0.0114794 *begin_x*begin_dx*begin_dy + -6.73759e-05 *begin_x*begin_y*begin_dx + 2.09426e-05 *lens_ipow(begin_x, 2)*begin_dy + -8.54798e-08 *lens_ipow(begin_x, 2)*begin_y + 9.34548e-05 *begin_y*lens_ipow(begin_lambda, 4);
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
    domega2_dx0[0][0] =  + -0.0115217  + 0.000652111 *lens_ipow(begin_dy, 2) + -0.00350055 *lens_ipow(begin_dx, 2) + -8.36284e-05 *begin_y*begin_dy + -2.83412e-07 *lens_ipow(begin_y, 2) + -0.000378262 *begin_x*begin_dx + -1.81021e-06 *lens_ipow(begin_x, 2)+0.0f;
    domega2_dx0[0][1] =  + -5.70223e-05 *begin_y*begin_dx + -8.36284e-05 *begin_x*begin_dy + -5.66825e-07 *begin_x*begin_y+0.0f;
    domega2_dx0[1][0] =  + 0.0114794 *begin_dx*begin_dy + -6.73759e-05 *begin_y*begin_dx + 4.18852e-05 *begin_x*begin_dy + -1.7096e-07 *begin_x*begin_y+0.0f;
    domega2_dx0[1][1] =  + -0.0123092  + 0.00443453 *lens_ipow(begin_dy, 2) + -0.00136161 *lens_ipow(begin_dx, 2) + -5.37288e-05 *begin_y*begin_dy + 2.42029e-07 *lens_ipow(begin_y, 2) + -6.73759e-05 *begin_x*begin_dx + -8.54798e-08 *lens_ipow(begin_x, 2) + 9.34548e-05 *lens_ipow(begin_lambda, 4)+0.0f;
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
  out[4] =  + 0.374643  + 0.275601 *begin_lambda + -4.11061e-06 *begin_dx + 3.29286e-07 *begin_y + -8.57117e-07 *begin_x + -0.0996656 *lens_ipow(begin_dy, 2) + 0.00109695 *begin_dx*begin_dy + -0.101182 *lens_ipow(begin_dx, 2) + -0.00183375 *begin_y*begin_dy + 1.73247e-05 *begin_y*begin_dx + -1.06392e-05 *lens_ipow(begin_y, 2) + 7.5642e-06 *begin_x*begin_dy + -0.00180111 *begin_x*begin_dx + 1.96049e-07 *begin_x*begin_y + -1.04885e-05 *lens_ipow(begin_x, 2) + -0.228709 *lens_ipow(begin_lambda, 3) + 0.000563113 *lens_ipow(begin_dy, 3) + -3.3104e-07 *begin_x*begin_y*begin_dx + -8.57044e-07 *lens_ipow(begin_y, 3)*begin_dy + -1.36808e-06 *begin_x*lens_ipow(begin_y, 2)*begin_dx + -1.52614e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2) + -1.27153e-06 *lens_ipow(begin_x, 3)*begin_dx + -1.25362e-06 *lens_ipow(begin_x, 2)*begin_dx*lens_ipow(begin_lambda, 2) + -2.03711e-06 *lens_ipow(begin_x, 2)*begin_y*begin_dy*begin_lambda + -4.04302 *lens_ipow(begin_dy, 6) + -13.0805 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 4) + -13.8032 *lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 2) + -4.28798 *lens_ipow(begin_dx, 6) + 4.87319e-08 *lens_ipow(begin_y, 4)*lens_ipow(begin_dy, 2) + -1.55187e-11 *lens_ipow(begin_y, 6) + -1.66807e-07 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2) + 5.68674e-07 *lens_ipow(begin_x, 3)*lens_ipow(begin_dy, 3) + -2.08832e-11 *lens_ipow(begin_x, 6) + -5.19397e-10 *begin_x*lens_ipow(begin_y, 4)*lens_ipow(begin_dx, 2) + 2.19224e-06 *begin_x*lens_ipow(begin_y, 3)*begin_dx*lens_ipow(begin_dy, 3) + -0.000158287 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 5) + -24.3205 *begin_x*lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 6) + -0.00470834 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dx, 6)*begin_dy + 0.409594 *lens_ipow(begin_lambda, 11) + -1.23177e-13 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 6)*lens_ipow(begin_lambda, 3);
else
  out[4] = 0.0f;
