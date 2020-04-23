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
       + 0.000135438  + 102.533 *begin_dx + -6.98227e-06 *begin_y + 0.852921 *begin_x + -0.000167245 *begin_y*begin_dx + 0.000124562 *begin_x*begin_dx + -2.71144e-06 *begin_x*begin_y + 93.3833 *begin_dx*lens_ipow(begin_dy, 2) + 93.9526 *lens_ipow(begin_dx, 3) + 2.77587 *begin_y*begin_dx*begin_dy + 0.0140524 *lens_ipow(begin_y, 2)*begin_dx + 1.30551 *begin_x*lens_ipow(begin_dy, 2) + 4.19041 *begin_x*lens_ipow(begin_dx, 2) + 0.0274398 *begin_x*begin_y*begin_dy + 3.07261e-05 *begin_x*lens_ipow(begin_y, 2) + 0.0423831 *lens_ipow(begin_x, 2)*begin_dx + 3.40673e-05 *lens_ipow(begin_x, 3) + -0.0148532 *begin_x*lens_ipow(begin_lambda, 3) + -0.637815 *begin_dx*lens_ipow(begin_lambda, 4) + 0.000421555 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 2) + 0.000503239 *lens_ipow(begin_x, 2)*begin_y*begin_dx*begin_dy + 0.000148145 *lens_ipow(begin_x, 3)*lens_ipow(begin_dy, 2) + 2.4049e-06 *lens_ipow(begin_x, 3)*begin_y*begin_dy + -5.91966 *begin_x*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2)*begin_lambda + -0.00269212 *lens_ipow(begin_x, 2)*begin_dx*lens_ipow(begin_lambda, 3) + -3385 *lens_ipow(begin_dx, 5)*lens_ipow(begin_dy, 2) + 0.634199 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 2) + -14.2498 *begin_x*lens_ipow(begin_dy, 6) + -9.02729 *begin_x*lens_ipow(begin_dx, 4)*lens_ipow(begin_lambda, 2) + 5.06598e-09 *begin_x*lens_ipow(begin_y, 5)*begin_dy + 7.75447e-07 *lens_ipow(begin_x, 5)*lens_ipow(begin_dx, 2) + -40421.2 *lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 6) + 0.760604 *begin_x*begin_y*lens_ipow(begin_dx, 4)*begin_dy*lens_ipow(begin_lambda, 2) + 1.50012e-09 *lens_ipow(begin_x, 5)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2) + -6.29746e-13 *lens_ipow(begin_x, 5)*lens_ipow(begin_y, 4) + -1.51618e-13 *lens_ipow(begin_x, 9) + -9863.75 *lens_ipow(begin_dx, 7)*lens_ipow(begin_lambda, 4) + -7.32514e-14 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 8)*begin_dx + 7.85729e-05 *lens_ipow(begin_x, 3)*lens_ipow(begin_lambda, 8) + 1.99145e-07 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 2)*lens_ipow(begin_lambda, 6),
       + 5.51359e-06  + 102.517 *begin_dy + 0.853298 *begin_y + -0.000190306 *begin_y*begin_dy + -2.89701e-06 *lens_ipow(begin_y, 2) + 94.8048 *lens_ipow(begin_dy, 3) + 95.0504 *lens_ipow(begin_dx, 2)*begin_dy + 4.14295 *begin_y*lens_ipow(begin_dy, 2) + 1.32585 *begin_y*lens_ipow(begin_dx, 2) + 0.0420256 *lens_ipow(begin_y, 2)*begin_dy + 3.41119e-05 *lens_ipow(begin_y, 3) + 2.799 *begin_x*begin_dx*begin_dy + 0.0277656 *begin_x*begin_y*begin_dx + 0.014278 *lens_ipow(begin_x, 2)*begin_dy + 3.06988e-05 *lens_ipow(begin_x, 2)*begin_y + -0.0174171 *begin_y*lens_ipow(begin_lambda, 3) + -0.786584 *begin_dy*lens_ipow(begin_lambda, 4) + 0.000499816 *begin_x*lens_ipow(begin_y, 2)*begin_dx*begin_dy + 0.000103369 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dy, 2) + -1290.57 *lens_ipow(begin_dx, 6)*begin_dy + 0.386262 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 3) + 1.11921e-06 *lens_ipow(begin_y, 5)*lens_ipow(begin_dy, 2) + -0.00518935 *lens_ipow(begin_x, 3)*lens_ipow(begin_dx, 3)*begin_dy + 1.88569e-06 *lens_ipow(begin_x, 4)*begin_y*lens_ipow(begin_dx, 2) + -6.00433 *begin_y*lens_ipow(begin_dy, 4)*lens_ipow(begin_lambda, 3) + -138.657 *begin_y*lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 2) + -178.463 *begin_y*lens_ipow(begin_dx, 8) + -0.0702456 *lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 6) + -1.49787e-13 *lens_ipow(begin_y, 9) + 2.60456 *begin_x*begin_y*begin_dx*lens_ipow(begin_dy, 4)*lens_ipow(begin_lambda, 2) + -1.38786e-07 *begin_x*lens_ipow(begin_y, 5)*begin_dx*lens_ipow(begin_dy, 2) + 5.40927e-05 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_lambda, 6) + -3.44904e-09 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 4)*lens_ipow(begin_dx, 2)*begin_dy + -6.85293e-13 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 5) + 4.39955e-11 *lens_ipow(begin_x, 7)*begin_y*begin_dx + -9.55545e-12 *lens_ipow(begin_x, 8)*begin_dy + -16582.9 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 5)*lens_ipow(begin_lambda, 3) + -48491.8 *lens_ipow(begin_dy, 9)*lens_ipow(begin_lambda, 2) + 0.000152731 *lens_ipow(begin_y, 3)*lens_ipow(begin_lambda, 8) + -7.96795e-11 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 4)*begin_dy*lens_ipow(begin_lambda, 2)
    };
    const float delta_ap[] = {ap_x - pred_ap[0], ap_y - pred_ap[1]};
    sqr_ap_err = delta_ap[0]*delta_ap[0]+delta_ap[1]*delta_ap[1];
    float dx1_domega0[2][2];
    dx1_domega0[0][0] =  + 102.533  + -0.000167245 *begin_y + 0.000124562 *begin_x + 93.3833 *lens_ipow(begin_dy, 2) + 281.858 *lens_ipow(begin_dx, 2) + 2.77587 *begin_y*begin_dy + 0.0140524 *lens_ipow(begin_y, 2) + 8.38082 *begin_x*begin_dx + 0.0423831 *lens_ipow(begin_x, 2) + -0.637815 *lens_ipow(begin_lambda, 4) + 0.000503239 *lens_ipow(begin_x, 2)*begin_y*begin_dy + -11.8393 *begin_x*begin_dx*lens_ipow(begin_dy, 2)*begin_lambda + -0.00269212 *lens_ipow(begin_x, 2)*lens_ipow(begin_lambda, 3) + -16925 *lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 2) + 1.9026 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + -36.1092 *begin_x*lens_ipow(begin_dx, 3)*lens_ipow(begin_lambda, 2) + 1.55089e-06 *lens_ipow(begin_x, 5)*begin_dx + -121264 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 6) + 3.04241 *begin_x*begin_y*lens_ipow(begin_dx, 3)*begin_dy*lens_ipow(begin_lambda, 2) + 3.00023e-09 *lens_ipow(begin_x, 5)*lens_ipow(begin_y, 2)*begin_dx + -69046.2 *lens_ipow(begin_dx, 6)*lens_ipow(begin_lambda, 4) + -7.32514e-14 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 8)+0.0f;
    dx1_domega0[0][1] =  + 186.767 *begin_dx*begin_dy + 2.77587 *begin_y*begin_dx + 2.61103 *begin_x*begin_dy + 0.0274398 *begin_x*begin_y + 0.000843109 *begin_x*lens_ipow(begin_y, 2)*begin_dy + 0.000503239 *lens_ipow(begin_x, 2)*begin_y*begin_dx + 0.00029629 *lens_ipow(begin_x, 3)*begin_dy + 2.4049e-06 *lens_ipow(begin_x, 3)*begin_y + -11.8393 *begin_x*lens_ipow(begin_dx, 2)*begin_dy*begin_lambda + -6770.01 *lens_ipow(begin_dx, 5)*begin_dy + 1.2684 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 3)*begin_dy + -85.4989 *begin_x*lens_ipow(begin_dy, 5) + 5.06598e-09 *begin_x*lens_ipow(begin_y, 5) + -242527 *lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 5) + 0.760604 *begin_x*begin_y*lens_ipow(begin_dx, 4)*lens_ipow(begin_lambda, 2)+0.0f;
    dx1_domega0[1][0] =  + 190.101 *begin_dx*begin_dy + 2.65171 *begin_y*begin_dx + 2.799 *begin_x*begin_dy + 0.0277656 *begin_x*begin_y + 0.000499816 *begin_x*lens_ipow(begin_y, 2)*begin_dy + -7743.45 *lens_ipow(begin_dx, 5)*begin_dy + 0.772525 *lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 3) + -0.015568 *lens_ipow(begin_x, 3)*lens_ipow(begin_dx, 2)*begin_dy + 3.77137e-06 *lens_ipow(begin_x, 4)*begin_y*begin_dx + -554.63 *begin_y*lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 2) + -1427.71 *begin_y*lens_ipow(begin_dx, 7) + 2.60456 *begin_x*begin_y*lens_ipow(begin_dy, 4)*lens_ipow(begin_lambda, 2) + -1.38786e-07 *begin_x*lens_ipow(begin_y, 5)*lens_ipow(begin_dy, 2) + -6.89809e-09 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 4)*begin_dx*begin_dy + 4.39955e-11 *lens_ipow(begin_x, 7)*begin_y + -33165.7 *begin_dx*lens_ipow(begin_dy, 5)*lens_ipow(begin_lambda, 3)+0.0f;
    dx1_domega0[1][1] =  + 102.517  + -0.000190306 *begin_y + 284.414 *lens_ipow(begin_dy, 2) + 95.0504 *lens_ipow(begin_dx, 2) + 8.28589 *begin_y*begin_dy + 0.0420256 *lens_ipow(begin_y, 2) + 2.799 *begin_x*begin_dx + 0.014278 *lens_ipow(begin_x, 2) + -0.786584 *lens_ipow(begin_lambda, 4) + 0.000499816 *begin_x*lens_ipow(begin_y, 2)*begin_dx + 0.000206738 *lens_ipow(begin_x, 2)*begin_y*begin_dy + -1290.57 *lens_ipow(begin_dx, 6) + 1.15879 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + 2.23843e-06 *lens_ipow(begin_y, 5)*begin_dy + -0.00518935 *lens_ipow(begin_x, 3)*lens_ipow(begin_dx, 3) + -24.0173 *begin_y*lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 3) + -277.315 *begin_y*lens_ipow(begin_dx, 4)*begin_dy*lens_ipow(begin_lambda, 2) + -0.421474 *lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 5) + 10.4182 *begin_x*begin_y*begin_dx*lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 2) + -2.77572e-07 *begin_x*lens_ipow(begin_y, 5)*begin_dx*begin_dy + -3.44904e-09 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 4)*lens_ipow(begin_dx, 2) + -9.55545e-12 *lens_ipow(begin_x, 8) + -82914.3 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 4)*lens_ipow(begin_lambda, 3) + -436427 *lens_ipow(begin_dy, 8)*lens_ipow(begin_lambda, 2) + -7.96795e-11 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 4)*lens_ipow(begin_lambda, 2)+0.0f;
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
    out[0] =  + 0.000254069  + 70.2621 *begin_dx + -2.46007e-05 *begin_y + -0.653785 *begin_x + -0.0235688 *begin_dx*begin_dy + -0.000406482 *begin_y*begin_dx + -4.64527e-06 *begin_x*begin_y + -19.2785 *begin_dx*lens_ipow(begin_dy, 2) + -19.8202 *lens_ipow(begin_dx, 3) + 0.271372 *begin_y*begin_dx*begin_dy + 0.00259451 *lens_ipow(begin_y, 2)*begin_dx + 0.241786 *begin_x*lens_ipow(begin_dy, 2) + 0.483651 *begin_x*lens_ipow(begin_dx, 2) + 0.00548664 *begin_x*begin_y*begin_dy + -0.000142526 *begin_x*lens_ipow(begin_y, 2) + 0.00788821 *lens_ipow(begin_x, 2)*begin_dx + -0.000143924 *lens_ipow(begin_x, 3) + 0.000122817 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2) + -0.0935941 *begin_x*lens_ipow(begin_lambda, 3) + -3.0171e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2) + 4.14781e-06 *lens_ipow(begin_y, 4)*begin_dx + -0.0798335 *begin_x*begin_y*lens_ipow(begin_dx, 2)*begin_dy + 0.000656544 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2) + 1.7139e-05 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*begin_dx + 0.000298764 *lens_ipow(begin_x, 3)*lens_ipow(begin_dy, 2) + 0.000338906 *lens_ipow(begin_x, 3)*lens_ipow(begin_dx, 2) + -5.89805e-05 *lens_ipow(begin_y, 4)*begin_dx*lens_ipow(begin_dy, 2) + -0.0934717 *begin_x*begin_y*lens_ipow(begin_dy, 5) + 0.773903 *lens_ipow(begin_x, 2)*begin_dx*lens_ipow(begin_dy, 4) + -0.00625468 *begin_x*begin_y*begin_dy*lens_ipow(begin_lambda, 5) + -2.47424 *begin_dx*lens_ipow(begin_lambda, 8) + 6776.95 *lens_ipow(begin_dx, 5)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 2) + 3.73333 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 7) + 7.79127e-11 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 5)*begin_dy + -6.40074e-09 *lens_ipow(begin_x, 6)*begin_y*begin_dx*begin_dy + 4.56499e-11 *lens_ipow(begin_x, 8)*begin_dx + -0.27136 *lens_ipow(begin_x, 2)*begin_dx*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 5) + -5.13597e-05 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_lambda, 5) + 0.440288 *begin_x*lens_ipow(begin_lambda, 10) + 1.02652 *begin_x*begin_y*lens_ipow(begin_dx, 2)*begin_dy*lens_ipow(begin_lambda, 6);
    out[1] =  + 0.000211885  + 70.217 *begin_dy + -0.00130993 *begin_dx + -0.654816 *begin_y + 0.0249546 *begin_dx*begin_dy + -4.52927e-06 *lens_ipow(begin_y, 2) + -2.06797e-06 *begin_x*begin_y + -15.9776 *lens_ipow(begin_dy, 3) + -18.8596 *lens_ipow(begin_dx, 2)*begin_dy + 0.661458 *begin_y*lens_ipow(begin_dy, 2) + 0.288332 *begin_y*lens_ipow(begin_dx, 2) + 0.0101164 *lens_ipow(begin_y, 2)*begin_dy + -0.000131313 *lens_ipow(begin_y, 3) + 0.3147 *begin_x*begin_dx*begin_dy + 0.00536003 *begin_x*begin_y*begin_dx + -2.17119e-07 *begin_x*lens_ipow(begin_y, 2) + 0.00272707 *lens_ipow(begin_x, 2)*begin_dy + -0.0001419 *lens_ipow(begin_x, 2)*begin_y + -0.0993731 *begin_y*lens_ipow(begin_lambda, 3) + 0.000963957 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dy, 2) + 1.87393e-05 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*begin_dy + 3.89315e-06 *lens_ipow(begin_x, 4)*begin_dy + -0.0193599 *lens_ipow(begin_y, 2)*begin_dy*lens_ipow(begin_lambda, 4) + -5.31401e-05 *lens_ipow(begin_x, 4)*lens_ipow(begin_dx, 2)*begin_dy + -1.35538e-07 *lens_ipow(begin_y, 5)*lens_ipow(begin_lambda, 3) + -0.0055037 *begin_x*begin_y*begin_dx*lens_ipow(begin_lambda, 5) + -2.27588 *begin_dy*lens_ipow(begin_lambda, 8) + 407.071 *lens_ipow(begin_dx, 4)*begin_dy*lens_ipow(begin_lambda, 4) + 2.25484e-11 *lens_ipow(begin_y, 8)*begin_dy + -1.96566e-05 *begin_x*lens_ipow(begin_y, 4)*lens_ipow(begin_dx, 3)*begin_dy + 6.16017 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 7) + 1.45081e-05 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 4) + -0.00206149 *lens_ipow(begin_x, 3)*begin_y*begin_dx*lens_ipow(begin_dy, 4) + 5.79454e-11 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 5)*begin_dx + 3.49016e-11 *lens_ipow(begin_x, 7)*begin_y*begin_dx + -13.0414 *begin_x*begin_dx*lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 5) + -2.49575e-05 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*begin_dy*lens_ipow(begin_lambda, 5) + -1071.41 *lens_ipow(begin_dy, 5)*lens_ipow(begin_lambda, 6) + 0.561412 *begin_y*lens_ipow(begin_lambda, 10) + -6.76347 *begin_y*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 8);
    out[2] =  + -0.257776 *begin_dx + -0.0118238 *begin_x + 0.21722 *begin_dx*lens_ipow(begin_dy, 2) + 0.222862 *lens_ipow(begin_dx, 3) + 0.00220759 *begin_y*begin_dx*begin_dy + 2.75933e-05 *lens_ipow(begin_y, 2)*begin_dx + 0.00054091 *begin_x*lens_ipow(begin_dy, 2) + 0.00289383 *begin_x*lens_ipow(begin_dx, 2) + 8.75148e-06 *begin_x*begin_y*begin_dy + 8.65239e-07 *begin_x*lens_ipow(begin_y, 2) + 3.32157e-05 *lens_ipow(begin_x, 2)*begin_dx + 9.092e-07 *lens_ipow(begin_x, 3) + -0.00679466 *begin_dx*lens_ipow(begin_lambda, 3) + 0.00030849 *begin_x*lens_ipow(begin_lambda, 3) + -1.68694e-06 *lens_ipow(begin_x, 2)*begin_y*begin_dx*begin_dy + 5.74153e-15 *lens_ipow(begin_x, 5)*lens_ipow(begin_y, 4) + -0.00187362 *begin_x*lens_ipow(begin_lambda, 10);
    out[3] =  + -0.257627 *begin_dy + -0.0118217 *begin_y + 0.220846 *lens_ipow(begin_dy, 3) + 0.215936 *lens_ipow(begin_dx, 2)*begin_dy + 0.00282102 *begin_y*lens_ipow(begin_dy, 2) + 0.00020588 *begin_y*lens_ipow(begin_dx, 2) + 3.24937e-05 *lens_ipow(begin_y, 2)*begin_dy + 9.08127e-07 *lens_ipow(begin_y, 3) + 0.0023944 *begin_x*begin_dx*begin_dy + 2.02132e-05 *lens_ipow(begin_x, 2)*begin_dy + 9.04483e-07 *lens_ipow(begin_x, 2)*begin_y + -0.00690557 *begin_dy*lens_ipow(begin_lambda, 3) + 0.000306695 *begin_y*lens_ipow(begin_lambda, 3) + -1.85119e-06 *begin_x*lens_ipow(begin_y, 2)*begin_dx*begin_dy + 6.21095e-15 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 5) + -0.00185606 *begin_y*lens_ipow(begin_lambda, 10);
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
    domega2_dx0[0][0] =  + -0.0118238  + 0.00054091 *lens_ipow(begin_dy, 2) + 0.00289383 *lens_ipow(begin_dx, 2) + 8.75148e-06 *begin_y*begin_dy + 8.65239e-07 *lens_ipow(begin_y, 2) + 6.64314e-05 *begin_x*begin_dx + 2.7276e-06 *lens_ipow(begin_x, 2) + 0.00030849 *lens_ipow(begin_lambda, 3) + -3.37388e-06 *begin_x*begin_y*begin_dx*begin_dy + 2.87077e-14 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 4) + -0.00187362 *lens_ipow(begin_lambda, 10)+0.0f;
    domega2_dx0[0][1] =  + 0.00220759 *begin_dx*begin_dy + 5.51866e-05 *begin_y*begin_dx + 8.75148e-06 *begin_x*begin_dy + 1.73048e-06 *begin_x*begin_y + -1.68694e-06 *lens_ipow(begin_x, 2)*begin_dx*begin_dy + 2.29661e-14 *lens_ipow(begin_x, 5)*lens_ipow(begin_y, 3)+0.0f;
    domega2_dx0[1][0] =  + 0.0023944 *begin_dx*begin_dy + 4.04264e-05 *begin_x*begin_dy + 1.80897e-06 *begin_x*begin_y + -1.85119e-06 *lens_ipow(begin_y, 2)*begin_dx*begin_dy + 2.48438e-14 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 5)+0.0f;
    domega2_dx0[1][1] =  + -0.0118217  + 0.00282102 *lens_ipow(begin_dy, 2) + 0.00020588 *lens_ipow(begin_dx, 2) + 6.49875e-05 *begin_y*begin_dy + 2.72438e-06 *lens_ipow(begin_y, 2) + 9.04483e-07 *lens_ipow(begin_x, 2) + 0.000306695 *lens_ipow(begin_lambda, 3) + -3.70238e-06 *begin_x*begin_y*begin_dx*begin_dy + 3.10547e-14 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 4) + -0.00185606 *lens_ipow(begin_lambda, 10)+0.0f;
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
  out[4] =  + 0.0478789  + 0.200766 *begin_lambda + -3.1003e-06 *begin_dx + 1.67699e-07 *begin_y + -1.46095e-07 *begin_x + -0.0293785 *lens_ipow(begin_dy, 2) + -0.0292605 *lens_ipow(begin_dx, 2) + -0.000296007 *begin_y*begin_dy + -2.94687e-06 *lens_ipow(begin_y, 2) + -0.000278624 *begin_x*begin_dx + -3.10664e-06 *lens_ipow(begin_x, 2) + -0.161495 *lens_ipow(begin_lambda, 3) + -1.43251e-05 *begin_x*begin_dx*begin_dy + 0.349992 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + -0.0187968 *begin_y*lens_ipow(begin_dy, 3) + -0.000397756 *lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 2) + -3.08e-06 *lens_ipow(begin_y, 3)*begin_dy + -0.00989689 *begin_x*begin_dx*lens_ipow(begin_dy, 2) + -0.0179809 *begin_x*lens_ipow(begin_dx, 3) + -0.000393196 *begin_x*begin_y*begin_dx*begin_dy + -2.14122e-06 *begin_x*lens_ipow(begin_y, 2)*begin_dx + -0.000390689 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 2) + -2.26345e-06 *lens_ipow(begin_x, 2)*begin_y*begin_dy + -2.02248e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2) + -2.98377e-06 *lens_ipow(begin_x, 3)*begin_dx + -3.59085 *lens_ipow(begin_dy, 6) + -2.37733 *lens_ipow(begin_dx, 6) + -0.216175 *begin_y*lens_ipow(begin_dx, 4)*begin_dy + -2.18891e-07 *lens_ipow(begin_y, 4)*lens_ipow(begin_dx, 2) + -2.09708e-11 *lens_ipow(begin_y, 6) + 3.83117e-12 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 4) + -3.09072e-07 *lens_ipow(begin_x, 4)*lens_ipow(begin_dy, 2) + -3.60435e-11 *lens_ipow(begin_x, 6)*begin_lambda + -226.248 *lens_ipow(begin_dx, 6)*lens_ipow(begin_dy, 2) + -7.43462e-10 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 4)*lens_ipow(begin_dy, 2) + -5.39322e-10 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2) + -2.2056e-11 *lens_ipow(begin_y, 6)*lens_ipow(begin_lambda, 3) + -4722.67 *lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 6) + -0.16562 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 8) + 0.278312 *lens_ipow(begin_lambda, 11);
else
  out[4] = 0.0f;
