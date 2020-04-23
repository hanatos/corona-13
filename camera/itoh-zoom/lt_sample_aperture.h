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
       + 0.00012656  + 0.00134211 *begin_dy + 188.258 *begin_dx + 0.763924 *begin_x + -2.6138e-06 *begin_x*begin_y + 32.2433 *begin_dx*lens_ipow(begin_dy, 2) + 36.1036 *lens_ipow(begin_dx, 3) + 2.01647 *begin_y*begin_dx*begin_dy + -0.00176571 *begin_y*lens_ipow(begin_dx, 2) + -2.19382e-05 *lens_ipow(begin_y, 2)*begin_dy + 0.0101 *lens_ipow(begin_y, 2)*begin_dx + 0.517272 *begin_x*lens_ipow(begin_dy, 2) + 2.67182 *begin_x*lens_ipow(begin_dx, 2) + 0.0162914 *begin_x*begin_y*begin_dy + 7.2713e-05 *begin_x*lens_ipow(begin_y, 2) + 0.0267636 *lens_ipow(begin_x, 2)*begin_dx + 8.0486e-05 *lens_ipow(begin_x, 3) + -15.854 *begin_dx*lens_ipow(begin_lambda, 3) + -0.000189382 *lens_ipow(begin_y, 2)*begin_dx*begin_dy + 0.000224024 *begin_x*begin_y*lens_ipow(begin_dy, 2) + -0.179877 *begin_x*lens_ipow(begin_lambda, 4) + 0.00219533 *begin_x*begin_y*begin_dx*lens_ipow(begin_dy, 2) + 0.000367953 *lens_ipow(begin_x, 2)*begin_y*begin_dx*begin_dy + 1.08213e-09 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 3) + -2.44349e-05 *lens_ipow(begin_x, 3)*lens_ipow(begin_lambda, 2) + 5.41245e-09 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 2) + -1598.94 *lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 2) + -1.15876 *begin_x*lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 4) + 4.26096e-05 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 3) + 5.25322e-05 *lens_ipow(begin_x, 4)*lens_ipow(begin_dx, 3) + -1004.87 *lens_ipow(begin_dx, 5)*lens_ipow(begin_lambda, 3) + -0.416428 *begin_y*begin_dx*begin_dy*lens_ipow(begin_lambda, 5) + -0.00390952 *begin_x*begin_y*begin_dy*lens_ipow(begin_lambda, 5) + 0.258253 *begin_x*begin_y*lens_ipow(begin_dx, 2)*begin_dy*lens_ipow(begin_lambda, 3) + 245261 *lens_ipow(begin_dx, 7)*lens_ipow(begin_dy, 2) + -1.57854 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 5)*lens_ipow(begin_lambda, 2) + 1.43705e-11 *lens_ipow(begin_x, 5)*lens_ipow(begin_y, 3)*begin_dy + 83.2242 *begin_dx*lens_ipow(begin_lambda, 10) + 0.885528 *begin_x*lens_ipow(begin_lambda, 10) + -5.6281e-05 *lens_ipow(begin_x, 4)*begin_dx*lens_ipow(begin_lambda, 6),
       + 0.000391942  + 188.207 *begin_dy + -0.00104105 *begin_dx + 0.764316 *begin_y + 5.52489e-05 *begin_x + -2.88233e-06 *lens_ipow(begin_y, 2) + 0.000120519 *begin_x*begin_dx + -1.32712e-06 *lens_ipow(begin_x, 2) + 33.607 *lens_ipow(begin_dy, 3) + 39.2076 *lens_ipow(begin_dx, 2)*begin_dy + 2.69578 *begin_y*lens_ipow(begin_dy, 2) + 0.564695 *begin_y*lens_ipow(begin_dx, 2) + 0.027093 *lens_ipow(begin_y, 2)*begin_dy + 7.42741e-05 *lens_ipow(begin_y, 3) + 2.05443 *begin_x*begin_dx*begin_dy + 0.0164373 *begin_x*begin_y*begin_dx + 0.0102186 *lens_ipow(begin_x, 2)*begin_dy + 7.43968e-05 *lens_ipow(begin_x, 2)*begin_y + -2.79695e-07 *lens_ipow(begin_x, 3) + -15.5452 *begin_dy*lens_ipow(begin_lambda, 3) + 0.000256818 *begin_x*begin_y*begin_dx*begin_dy + -0.188057 *begin_y*lens_ipow(begin_lambda, 4) + 0.0429834 *lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 3) + 0.000577833 *begin_x*lens_ipow(begin_y, 2)*begin_dx*begin_dy + -3.88665e-09 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 3) + 0.151664 *begin_x*begin_y*begin_dx*lens_ipow(begin_dy, 2)*begin_lambda + -1603.66 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 2) + -754.094 *lens_ipow(begin_dx, 4)*begin_dy*lens_ipow(begin_lambda, 2) + -1.35088 *begin_y*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 4) + -0.00800449 *lens_ipow(begin_y, 2)*begin_dy*lens_ipow(begin_lambda, 4) + 7.55687e-07 *lens_ipow(begin_y, 5)*lens_ipow(begin_dy, 2) + -44.4892 *begin_y*lens_ipow(begin_dx, 6)*begin_lambda + -0.000915322 *lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 3) + -0.816715 *begin_x*begin_dx*begin_dy*lens_ipow(begin_lambda, 5) + 1.2861e-09 *begin_x*lens_ipow(begin_y, 5)*begin_dx*begin_lambda + -0.716688 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 5)*begin_lambda + -1767.98 *lens_ipow(begin_dy, 5)*lens_ipow(begin_lambda, 4) + -0.00905914 *begin_x*begin_y*begin_dx*lens_ipow(begin_lambda, 6) + 82.1062 *begin_dy*lens_ipow(begin_lambda, 10) + 0.908817 *begin_y*lens_ipow(begin_lambda, 10)
    };
    const float delta_ap[] = {ap_x - pred_ap[0], ap_y - pred_ap[1]};
    sqr_ap_err = delta_ap[0]*delta_ap[0]+delta_ap[1]*delta_ap[1];
    float dx1_domega0[2][2];
    dx1_domega0[0][0] =  + 188.258  + 32.2433 *lens_ipow(begin_dy, 2) + 108.311 *lens_ipow(begin_dx, 2) + 2.01647 *begin_y*begin_dy + -0.00353142 *begin_y*begin_dx + 0.0101 *lens_ipow(begin_y, 2) + 5.34365 *begin_x*begin_dx + 0.0267636 *lens_ipow(begin_x, 2) + -15.854 *lens_ipow(begin_lambda, 3) + -0.000189382 *lens_ipow(begin_y, 2)*begin_dy + 0.00219533 *begin_x*begin_y*lens_ipow(begin_dy, 2) + 0.000367953 *lens_ipow(begin_x, 2)*begin_y*begin_dy + -4796.81 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 2) + -2.31752 *begin_x*begin_dx*lens_ipow(begin_lambda, 4) + 0.000127829 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2) + 0.000157597 *lens_ipow(begin_x, 4)*lens_ipow(begin_dx, 2) + -5024.37 *lens_ipow(begin_dx, 4)*lens_ipow(begin_lambda, 3) + -0.416428 *begin_y*begin_dy*lens_ipow(begin_lambda, 5) + 0.516505 *begin_x*begin_y*begin_dx*begin_dy*lens_ipow(begin_lambda, 3) + 1.71683e+06 *lens_ipow(begin_dx, 6)*lens_ipow(begin_dy, 2) + -7.89271 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 4)*lens_ipow(begin_lambda, 2) + 83.2242 *lens_ipow(begin_lambda, 10) + -5.6281e-05 *lens_ipow(begin_x, 4)*lens_ipow(begin_lambda, 6)+0.0f;
    dx1_domega0[0][1] =  + 0.00134211  + 64.4866 *begin_dx*begin_dy + 2.01647 *begin_y*begin_dx + -2.19382e-05 *lens_ipow(begin_y, 2) + 1.03454 *begin_x*begin_dy + 0.0162914 *begin_x*begin_y + -0.000189382 *lens_ipow(begin_y, 2)*begin_dx + 0.000448047 *begin_x*begin_y*begin_dy + 0.00439065 *begin_x*begin_y*begin_dx*begin_dy + 0.000367953 *lens_ipow(begin_x, 2)*begin_y*begin_dx + -3197.87 *lens_ipow(begin_dx, 3)*begin_dy*lens_ipow(begin_lambda, 2) + -0.416428 *begin_y*begin_dx*lens_ipow(begin_lambda, 5) + -0.00390952 *begin_x*begin_y*lens_ipow(begin_lambda, 5) + 0.258253 *begin_x*begin_y*lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 3) + 490523 *lens_ipow(begin_dx, 7)*begin_dy + 1.43705e-11 *lens_ipow(begin_x, 5)*lens_ipow(begin_y, 3)+0.0f;
    dx1_domega0[1][0] =  + -0.00104105  + 0.000120519 *begin_x + 78.4152 *begin_dx*begin_dy + 1.12939 *begin_y*begin_dx + 2.05443 *begin_x*begin_dy + 0.0164373 *begin_x*begin_y + 0.000256818 *begin_x*begin_y*begin_dy + 0.000577833 *begin_x*lens_ipow(begin_y, 2)*begin_dy + 0.151664 *begin_x*begin_y*lens_ipow(begin_dy, 2)*begin_lambda + -3207.33 *begin_dx*lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 2) + -3016.38 *lens_ipow(begin_dx, 3)*begin_dy*lens_ipow(begin_lambda, 2) + -266.935 *begin_y*lens_ipow(begin_dx, 5)*begin_lambda + -0.00183064 *lens_ipow(begin_y, 3)*begin_dx*lens_ipow(begin_lambda, 3) + -0.816715 *begin_x*begin_dy*lens_ipow(begin_lambda, 5) + 1.2861e-09 *begin_x*lens_ipow(begin_y, 5)*begin_lambda + -0.00905914 *begin_x*begin_y*lens_ipow(begin_lambda, 6)+0.0f;
    dx1_domega0[1][1] =  + 188.207  + 100.821 *lens_ipow(begin_dy, 2) + 39.2076 *lens_ipow(begin_dx, 2) + 5.39155 *begin_y*begin_dy + 0.027093 *lens_ipow(begin_y, 2) + 2.05443 *begin_x*begin_dx + 0.0102186 *lens_ipow(begin_x, 2) + -15.5452 *lens_ipow(begin_lambda, 3) + 0.000256818 *begin_x*begin_y*begin_dx + 0.12895 *lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 2) + 0.000577833 *begin_x*lens_ipow(begin_y, 2)*begin_dx + 0.303327 *begin_x*begin_y*begin_dx*begin_dy*begin_lambda + -4810.99 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 2) + -754.094 *lens_ipow(begin_dx, 4)*lens_ipow(begin_lambda, 2) + -2.70175 *begin_y*begin_dy*lens_ipow(begin_lambda, 4) + -0.00800449 *lens_ipow(begin_y, 2)*lens_ipow(begin_lambda, 4) + 1.51137e-06 *lens_ipow(begin_y, 5)*begin_dy + -0.816715 *begin_x*begin_dx*lens_ipow(begin_lambda, 5) + -3.58344 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 4)*begin_lambda + -8839.92 *lens_ipow(begin_dy, 4)*lens_ipow(begin_lambda, 4) + 82.1062 *lens_ipow(begin_lambda, 10)+0.0f;
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
    out[0] =  + 0.000853931  + 97.9961 *begin_dx + 1.97042e-05 *begin_y + -0.37657 *begin_x + -0.0470299 *lens_ipow(begin_dy, 2) + -0.0600743 *lens_ipow(begin_dx, 2) + -2.53798e-06 *begin_x*begin_y + 1.87661 *begin_dx*lens_ipow(begin_lambda, 2) + 27.8227 *begin_dx*lens_ipow(begin_dy, 2) + 31.7923 *lens_ipow(begin_dx, 3) + 0.931397 *begin_y*begin_dx*begin_dy + 0.00703035 *lens_ipow(begin_y, 2)*begin_dx + 0.643104 *begin_x*lens_ipow(begin_dy, 2) + 2.12778 *begin_x*lens_ipow(begin_dx, 2) + 0.0132105 *begin_x*begin_y*begin_dy + 6.56453e-05 *begin_x*lens_ipow(begin_y, 2) + 0.0201332 *lens_ipow(begin_x, 2)*begin_dx + 6.78108e-05 *lens_ipow(begin_x, 3) + -0.100499 *begin_x*lens_ipow(begin_lambda, 3) + 0.0364028 *lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 2) + 0.232962 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 3) + 0.00175695 *lens_ipow(begin_x, 2)*begin_y*begin_dx*begin_dy + -9.79743e-09 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 2) + -3.20268 *begin_x*lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 3) + -7936.74 *lens_ipow(begin_dx, 7) + 118.81 *begin_y*lens_ipow(begin_dx, 5)*begin_dy + 0.83254 *begin_x*begin_y*lens_ipow(begin_dx, 2)*begin_dy*lens_ipow(begin_lambda, 2) + 0.000425213 *lens_ipow(begin_x, 4)*begin_dx*lens_ipow(begin_dy, 2) + 4.47562e-06 *lens_ipow(begin_x, 5)*lens_ipow(begin_dx, 2) + -1.51328 *begin_y*begin_dx*begin_dy*lens_ipow(begin_lambda, 5) + -0.0187963 *begin_x*begin_y*begin_dy*lens_ipow(begin_lambda, 5) + -183344 *lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 4)*lens_ipow(begin_lambda, 2) + -0.0460233 *lens_ipow(begin_x, 2)*begin_dx*lens_ipow(begin_lambda, 6) + -0.0176533 *lens_ipow(begin_x, 3)*lens_ipow(begin_dy, 6) + -7416.37 *lens_ipow(begin_dx, 5)*lens_ipow(begin_lambda, 5) + -3.80988e-10 *lens_ipow(begin_x, 7)*lens_ipow(begin_lambda, 3) + -0.59839 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 3)*lens_ipow(begin_lambda, 6) + 0.534569 *begin_x*lens_ipow(begin_lambda, 10) + -69.3283 *begin_x*lens_ipow(begin_dy, 4)*lens_ipow(begin_lambda, 6) + -1.3621 *lens_ipow(begin_x, 2)*begin_dx*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 6);
    out[1] =  + 0.000137611  + 97.9425 *begin_dy + -0.000754665 *begin_dx + -0.374459 *begin_y + 0.000321517 *begin_y*begin_dy + 0.000273961 *begin_y*begin_dx + 4.24285e-06 *begin_x*begin_y + 2.18591 *begin_dy*lens_ipow(begin_lambda, 2) + 27.8545 *lens_ipow(begin_dy, 3) + 33.2953 *lens_ipow(begin_dx, 2)*begin_dy + 1.67327 *begin_y*lens_ipow(begin_dy, 2) + 0.65598 *begin_y*lens_ipow(begin_dx, 2) + 0.019324 *lens_ipow(begin_y, 2)*begin_dy + 6.06604e-05 *lens_ipow(begin_y, 3) + 0.958745 *begin_x*begin_dx*begin_dy + 0.0132771 *begin_x*begin_y*begin_dx + 2.38106e-07 *begin_x*lens_ipow(begin_y, 2) + 0.00690698 *lens_ipow(begin_x, 2)*begin_dy + 6.48198e-05 *lens_ipow(begin_x, 2)*begin_y + -0.111376 *begin_y*lens_ipow(begin_lambda, 3) + 0.00170548 *begin_x*lens_ipow(begin_y, 2)*begin_dx*begin_dy + -6.10793e-09 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 3) + 0.538828 *lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 3)*begin_lambda + -26711.1 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 5) + 7.99012e-06 *lens_ipow(begin_y, 5)*lens_ipow(begin_dy, 2) + -0.996917 *begin_x*begin_dx*begin_dy*lens_ipow(begin_lambda, 4) + 0.71917 *begin_x*begin_y*begin_dx*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 2) + 3.09036 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 3) + -0.017372 *begin_x*begin_y*begin_dx*lens_ipow(begin_lambda, 5) + -106612 *lens_ipow(begin_dy, 9) + -180.329 *begin_y*lens_ipow(begin_dx, 6)*lens_ipow(begin_lambda, 2) + 7.24363e-07 *lens_ipow(begin_y, 6)*lens_ipow(begin_dx, 2)*begin_dy + -10205.4 *lens_ipow(begin_dy, 5)*lens_ipow(begin_lambda, 5) + -0.0877426 *lens_ipow(begin_y, 2)*begin_dy*lens_ipow(begin_lambda, 7) + -499466 *lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 4) + 0.617364 *begin_y*lens_ipow(begin_lambda, 10) + -13.6221 *begin_y*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 8) + -0.00546888 *lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 6) + 2.96072e-13 *lens_ipow(begin_y, 10)*begin_dy + -7.57184 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 5)*lens_ipow(begin_lambda, 4);
    out[2] =  + -5.49498e-07  + -0.349126 *begin_dx + -0.00884588 *begin_x + 0.846317 *begin_dx*lens_ipow(begin_dy, 2) + 0.82841 *lens_ipow(begin_dx, 3) + 0.00840566 *begin_y*begin_dx*begin_dy + 2.65844e-05 *lens_ipow(begin_y, 2)*begin_dx + 0.00297861 *begin_x*lens_ipow(begin_dy, 2) + 0.0124968 *begin_x*lens_ipow(begin_dx, 2) + 3.17635e-05 *begin_x*begin_y*begin_dy + 3.14338e-07 *begin_x*lens_ipow(begin_y, 2) + 3.02677e-05 *lens_ipow(begin_x, 2)*begin_dx + 0.00039821 *begin_x*lens_ipow(begin_lambda, 3) + 7.27783e-07 *lens_ipow(begin_x, 3)*begin_lambda + 0.0195371 *begin_dx*lens_ipow(begin_lambda, 4) + 0.000921139 *lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 2) + 2.6712e-06 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 2) + 0.00211755 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 3) + 8.03988e-06 *lens_ipow(begin_x, 2)*begin_y*begin_dx*begin_dy + -0.0123753 *begin_y*begin_dx*begin_dy*lens_ipow(begin_lambda, 3) + -56.6621 *lens_ipow(begin_dx, 7) + -4.7201e-05 *lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_lambda, 4) + 1.12933e-08 *lens_ipow(begin_y, 5)*begin_dx*begin_dy + -0.0317466 *begin_x*lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 4) + -0.00014451 *begin_x*begin_y*begin_dy*lens_ipow(begin_lambda, 4) + 0.0293829 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 2) + -8.97827e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 3) + 6.2369e-08 *lens_ipow(begin_x, 5)*lens_ipow(begin_dx, 2) + -5064.61 *lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 6) + 0.0113469 *begin_x*begin_y*lens_ipow(begin_dx, 2)*begin_dy*lens_ipow(begin_lambda, 4) + -9.15396e-12 *lens_ipow(begin_x, 7)*lens_ipow(begin_dy, 2) + 8.64732e-13 *lens_ipow(begin_x, 8)*begin_dx + -71.0647 *begin_dx*lens_ipow(begin_dy, 4)*lens_ipow(begin_lambda, 5) + -50.6656 *lens_ipow(begin_dx, 5)*lens_ipow(begin_lambda, 5) + -2976.11 *lens_ipow(begin_dx, 5)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 4) + -0.023315 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 5)*lens_ipow(begin_lambda, 4) + -0.804972 *begin_x*lens_ipow(begin_dy, 4)*lens_ipow(begin_lambda, 6) + -0.000776015 *lens_ipow(begin_x, 2)*begin_dx*lens_ipow(begin_lambda, 8) + -0.0041045 *lens_ipow(begin_x, 2)*begin_dx*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 6) + -6.4823e-06 *lens_ipow(begin_x, 3)*lens_ipow(begin_lambda, 8);
    out[3] =  + 2.40458e-06  + -0.349228 *begin_dy + -0.00885275 *begin_y + -2.73429e-07 *begin_x + 2.44939e-06 *begin_y*begin_dy + 1.80382e-08 *begin_x*begin_y + 0.840591 *lens_ipow(begin_dy, 3) + 0.905061 *lens_ipow(begin_dx, 2)*begin_dy + 0.013117 *begin_y*lens_ipow(begin_dy, 2) + 0.00302914 *begin_y*lens_ipow(begin_dx, 2) + 3.4731e-05 *lens_ipow(begin_y, 2)*begin_dy + 0.00785333 *begin_x*begin_dx*begin_dy + 2.48971e-05 *begin_x*begin_y*begin_dx + 2.13674e-05 *lens_ipow(begin_x, 2)*begin_dy + 3.44812e-07 *lens_ipow(begin_x, 2)*begin_y + 0.000419495 *begin_y*lens_ipow(begin_lambda, 3) + 7.40496e-07 *lens_ipow(begin_y, 3)*begin_lambda + 0.0163011 *begin_dy*lens_ipow(begin_lambda, 4) + 0.00210992 *lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 3) + 1.3325e-05 *begin_x*lens_ipow(begin_y, 2)*begin_dx*begin_dy + 1.85897e-11 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 2) + -57.1146 *lens_ipow(begin_dx, 6)*begin_dy + -0.0326856 *begin_y*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 4) + 5.34934e-08 *lens_ipow(begin_y, 5)*lens_ipow(begin_dy, 2) + -0.00825709 *begin_x*begin_dx*begin_dy*lens_ipow(begin_lambda, 4) + -0.00011687 *begin_x*begin_y*begin_dx*lens_ipow(begin_lambda, 4) + -0.0025329 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 5) + 4.07096e-06 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*begin_dy + 0.0114687 *begin_x*begin_y*begin_dx*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 3) + -1323.42 *lens_ipow(begin_dy, 9) + 6.28052e-09 *lens_ipow(begin_y, 6)*lens_ipow(begin_dx, 2)*begin_dy + -50.7618 *lens_ipow(begin_dy, 5)*lens_ipow(begin_lambda, 5) + 1.42775e-07 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 2)*begin_dx*begin_dy*lens_ipow(begin_lambda, 3) + -2611.24 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 5)*lens_ipow(begin_lambda, 4) + -0.616453 *begin_y*lens_ipow(begin_dx, 4)*lens_ipow(begin_lambda, 6) + 28.2314 *begin_y*lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 4) + -0.000774292 *lens_ipow(begin_y, 2)*begin_dy*lens_ipow(begin_lambda, 8) + -6.7578e-06 *lens_ipow(begin_y, 3)*lens_ipow(begin_lambda, 8) + 1.87035e-15 *lens_ipow(begin_y, 10)*begin_dy + -0.00363419 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 6);
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
    domega2_dx0[0][0] =  + -0.00884588  + 0.00297861 *lens_ipow(begin_dy, 2) + 0.0124968 *lens_ipow(begin_dx, 2) + 3.17635e-05 *begin_y*begin_dy + 3.14338e-07 *lens_ipow(begin_y, 2) + 6.05354e-05 *begin_x*begin_dx + 0.00039821 *lens_ipow(begin_lambda, 3) + 2.18335e-06 *lens_ipow(begin_x, 2)*begin_lambda + 2.6712e-06 *lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 2) + 0.0042351 *begin_x*lens_ipow(begin_dx, 3) + 1.60798e-05 *begin_x*begin_y*begin_dx*begin_dy + -0.0317466 *lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 4) + -0.00014451 *begin_y*begin_dy*lens_ipow(begin_lambda, 4) + 0.0587658 *begin_x*lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 2) + -1.79565e-07 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 3) + 3.11845e-07 *lens_ipow(begin_x, 4)*lens_ipow(begin_dx, 2) + 0.0113469 *begin_y*lens_ipow(begin_dx, 2)*begin_dy*lens_ipow(begin_lambda, 4) + -6.40777e-11 *lens_ipow(begin_x, 6)*lens_ipow(begin_dy, 2) + 6.91786e-12 *lens_ipow(begin_x, 7)*begin_dx + -0.804972 *lens_ipow(begin_dy, 4)*lens_ipow(begin_lambda, 6) + -0.00155203 *begin_x*begin_dx*lens_ipow(begin_lambda, 8) + -0.00820899 *begin_x*begin_dx*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 6) + -1.94469e-05 *lens_ipow(begin_x, 2)*lens_ipow(begin_lambda, 8)+0.0f;
    domega2_dx0[0][1] =  + 0.00840566 *begin_dx*begin_dy + 5.31689e-05 *begin_y*begin_dx + 3.17635e-05 *begin_x*begin_dy + 6.28676e-07 *begin_x*begin_y + 0.00184228 *begin_y*begin_dx*lens_ipow(begin_dy, 2) + 5.3424e-06 *begin_x*begin_y*lens_ipow(begin_dy, 2) + 8.03988e-06 *lens_ipow(begin_x, 2)*begin_dx*begin_dy + -0.0123753 *begin_dx*begin_dy*lens_ipow(begin_lambda, 3) + -9.4402e-05 *begin_y*begin_dx*lens_ipow(begin_lambda, 4) + 5.64665e-08 *lens_ipow(begin_y, 4)*begin_dx*begin_dy + -0.00014451 *begin_x*begin_dy*lens_ipow(begin_lambda, 4) + -1.79565e-07 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dx, 3) + 0.0113469 *begin_x*lens_ipow(begin_dx, 2)*begin_dy*lens_ipow(begin_lambda, 4) + -0.04663 *begin_y*lens_ipow(begin_dx, 5)*lens_ipow(begin_lambda, 4)+0.0f;
    domega2_dx0[1][0] =  + -2.73429e-07  + 1.80382e-08 *begin_y + 0.00785333 *begin_dx*begin_dy + 2.48971e-05 *begin_y*begin_dx + 4.27347e-05 *begin_x*begin_dy + 6.89624e-07 *begin_x*begin_y + 1.3325e-05 *lens_ipow(begin_y, 2)*begin_dx*begin_dy + 5.57691e-11 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2) + -0.00825709 *begin_dx*begin_dy*lens_ipow(begin_lambda, 4) + -0.00011687 *begin_y*begin_dx*lens_ipow(begin_lambda, 4) + -0.0050658 *begin_x*lens_ipow(begin_dy, 5) + 8.14192e-06 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*begin_dy + 0.0114687 *begin_y*begin_dx*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 3) + 4.28324e-07 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*begin_dx*begin_dy*lens_ipow(begin_lambda, 3) + -0.00726839 *begin_x*lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 6)+0.0f;
    domega2_dx0[1][1] =  + -0.00885275  + 2.44939e-06 *begin_dy + 1.80382e-08 *begin_x + 0.013117 *lens_ipow(begin_dy, 2) + 0.00302914 *lens_ipow(begin_dx, 2) + 6.9462e-05 *begin_y*begin_dy + 2.48971e-05 *begin_x*begin_dx + 3.44812e-07 *lens_ipow(begin_x, 2) + 0.000419495 *lens_ipow(begin_lambda, 3) + 2.22149e-06 *lens_ipow(begin_y, 2)*begin_lambda + 0.00421983 *begin_y*lens_ipow(begin_dy, 3) + 2.66499e-05 *begin_x*begin_y*begin_dx*begin_dy + 3.71794e-11 *lens_ipow(begin_x, 3)*begin_y + -0.0326856 *lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 4) + 2.67467e-07 *lens_ipow(begin_y, 4)*lens_ipow(begin_dy, 2) + -0.00011687 *begin_x*begin_dx*lens_ipow(begin_lambda, 4) + 8.14192e-06 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dx, 2)*begin_dy + 0.0114687 *begin_x*begin_dx*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 3) + 3.76831e-08 *lens_ipow(begin_y, 5)*lens_ipow(begin_dx, 2)*begin_dy + 2.85549e-07 *lens_ipow(begin_x, 3)*begin_y*begin_dx*begin_dy*lens_ipow(begin_lambda, 3) + -0.616453 *lens_ipow(begin_dx, 4)*lens_ipow(begin_lambda, 6) + 28.2314 *lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 4) + -0.00154858 *begin_y*begin_dy*lens_ipow(begin_lambda, 8) + -2.02734e-05 *lens_ipow(begin_y, 2)*lens_ipow(begin_lambda, 8) + 1.87035e-14 *lens_ipow(begin_y, 9)*begin_dy+0.0f;
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
  out[4] =  + 0.148064  + 0.301995 *begin_lambda + -2.0734e-05 *begin_dx + 1.70581e-07 *begin_y + 3.97549e-08 *begin_x + -0.0112984 *lens_ipow(begin_dy, 2) + -0.013445 *lens_ipow(begin_dx, 2) + -1.30829e-05 *begin_y*begin_dy + -2.1516e-07 *lens_ipow(begin_y, 2) + -1.71305e-05 *begin_x*begin_dx + 9.71968e-09 *begin_x*begin_y + -1.0729e-07 *lens_ipow(begin_x, 2) + -0.246483 *lens_ipow(begin_lambda, 3) + 0.00133384 *lens_ipow(begin_dx, 2)*begin_dy + 6.88755e-05 *begin_y*begin_dx*begin_dy + -2.77381 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + -0.0238315 *begin_y*lens_ipow(begin_dy, 3) + -0.0272569 *begin_y*lens_ipow(begin_dx, 2)*begin_dy + -0.000224861 *lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 2) + -7.86473e-05 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2) + -1.00716e-06 *lens_ipow(begin_y, 3)*begin_dy + -0.0260068 *begin_x*begin_dx*lens_ipow(begin_dy, 2) + -0.0232514 *begin_x*lens_ipow(begin_dx, 3) + -0.000324813 *begin_x*begin_y*begin_dx*begin_dy + -1.11621e-06 *begin_x*lens_ipow(begin_y, 2)*begin_dx + -8.03731e-05 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 2) + -0.000228197 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 2) + -1.09579e-06 *lens_ipow(begin_x, 2)*begin_y*begin_dy + -4.38463e-09 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2) + -1.06634e-06 *lens_ipow(begin_x, 3)*begin_dx + -37.5481 *lens_ipow(begin_dy, 6) + -31.8755 *lens_ipow(begin_dx, 6) + -1.16464e-08 *lens_ipow(begin_y, 4)*lens_ipow(begin_dx, 2) + -4.61178e-12 *lens_ipow(begin_y, 6) + -0.0365375 *begin_x*begin_dx*lens_ipow(begin_dy, 4) + -6.25054e-12 *lens_ipow(begin_x, 6) + 0.0742931 *begin_dx*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 5) + -104.61 *begin_y*lens_ipow(begin_dy, 9) + 0.940835 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 8) + 0.430603 *lens_ipow(begin_lambda, 11);
else
  out[4] = 0.0f;
