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
       + 1.47104e-05  + 25.2552 *begin_dx + 1.64108e-06 *begin_y + 0.88529 *begin_x + 5.31843e-07 *begin_x*begin_y + -1.83617 *begin_dx*lens_ipow(begin_dy, 2) + -1.52153 *lens_ipow(begin_dx, 3) + 0.147112 *begin_y*begin_dx*begin_dy + 0.00609101 *lens_ipow(begin_y, 2)*begin_dx + -0.0554282 *begin_x*lens_ipow(begin_dy, 2) + 0.137339 *begin_x*lens_ipow(begin_dx, 2) + 7.14484e-05 *begin_x*lens_ipow(begin_y, 2) + 0.0100807 *lens_ipow(begin_x, 2)*begin_dx + 0.000107588 *lens_ipow(begin_x, 3) + 0.874804 *begin_dx*lens_ipow(begin_lambda, 3) + 0.00494302 *begin_x*begin_y*begin_dy*begin_lambda + 0.29696 *lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 2) + 0.0355258 *begin_x*lens_ipow(begin_lambda, 4) + -3.69121e-05 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 2) + 1.24453e-06 *lens_ipow(begin_x, 3)*begin_y*begin_dy + 0.0354535 *begin_y*begin_dx*begin_dy*lens_ipow(begin_lambda, 3) + 0.000170763 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_lambda, 3) + 0.000683737 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 5) + -2.09067e-05 *lens_ipow(begin_y, 3)*begin_dx*lens_ipow(begin_dy, 3) + 0.0926377 *begin_x*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 4) + 0.0256521 *begin_x*lens_ipow(begin_dy, 6) + 1.45029e-06 *lens_ipow(begin_x, 4)*lens_ipow(begin_dx, 3) + 0.44816 *begin_dx*lens_ipow(begin_dy, 8) + 0.0372671 *lens_ipow(begin_dx, 9) + 0.00160921 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 3)*lens_ipow(begin_lambda, 4) + 6.7244e-09 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 2)*begin_dy + 3.01533e-08 *lens_ipow(begin_x, 5)*lens_ipow(begin_dy, 4) + -1.33148e-12 *lens_ipow(begin_x, 8)*begin_dx + -4.65791 *begin_dx*lens_ipow(begin_lambda, 10) + 3.50038 *begin_dx*lens_ipow(begin_dy, 4)*lens_ipow(begin_lambda, 6) + -0.189056 *begin_x*lens_ipow(begin_lambda, 10) + -0.00423704 *begin_x*begin_y*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 5)*lens_ipow(begin_lambda, 2) + -0.000459711 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_lambda, 8) + 0.000101431 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dx, 7)*begin_dy + 6.8061e-13 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 6)*lens_ipow(begin_lambda, 2),
       + -1.66982e-05  + 25.2549 *begin_dy + 0.884941 *begin_y + 2.46235e-06 *begin_x + -0.000487611 *begin_dx*begin_dy + -1.46287 *lens_ipow(begin_dy, 3) + -1.45263 *lens_ipow(begin_dx, 2)*begin_dy + 0.144255 *begin_y*lens_ipow(begin_dy, 2) + -0.0348381 *begin_y*lens_ipow(begin_dx, 2) + 0.010328 *lens_ipow(begin_y, 2)*begin_dy + 0.000109554 *lens_ipow(begin_y, 3) + 0.180011 *begin_x*begin_dx*begin_dy + 0.00367345 *begin_x*begin_y*begin_dx + 0.00654957 *lens_ipow(begin_x, 2)*begin_dy + 0.000108976 *lens_ipow(begin_x, 2)*begin_y + 0.797561 *begin_dy*lens_ipow(begin_lambda, 3) + -6.53929e-05 *begin_x*lens_ipow(begin_dy, 3) + 0.0362703 *begin_y*lens_ipow(begin_lambda, 4) + -0.0176623 *begin_x*begin_dx*lens_ipow(begin_dy, 3) + 9.38447e-07 *begin_x*lens_ipow(begin_y, 3)*begin_dx + -1.53514e-05 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dx, 2) + 0.0188592 *begin_y*lens_ipow(begin_dx, 6) + 1.55345e-06 *lens_ipow(begin_y, 4)*lens_ipow(begin_dy, 3) + 0.000537731 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 5) + -5.72727e-05 *lens_ipow(begin_x, 3)*lens_ipow(begin_dx, 3)*begin_dy + 0.0886898 *lens_ipow(begin_dy, 5)*lens_ipow(begin_lambda, 3) + 0.0879175 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 5)*begin_lambda + 0.016735 *lens_ipow(begin_dy, 9) + 0.183099 *lens_ipow(begin_dx, 8)*begin_dy + 2.19235e-08 *lens_ipow(begin_y, 5)*lens_ipow(begin_dx, 4) + -1.50573e-12 *lens_ipow(begin_y, 8)*begin_dy + -0.0117741 *begin_x*lens_ipow(begin_dx, 3)*begin_dy*lens_ipow(begin_lambda, 4) + 8.16325e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 3)*lens_ipow(begin_lambda, 4) + 4.36515e-09 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 3)*begin_dx*lens_ipow(begin_dy, 2) + -3.22089e-12 *lens_ipow(begin_x, 8)*begin_dy + -4.37325 *begin_dy*lens_ipow(begin_lambda, 10) + -0.662976 *lens_ipow(begin_dx, 6)*lens_ipow(begin_dy, 5) + -0.183804 *begin_y*lens_ipow(begin_lambda, 10) + 9.05996e-05 *begin_x*lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 7) + -1.27567e-14 *lens_ipow(begin_x, 6)*lens_ipow(begin_y, 3)*lens_ipow(begin_lambda, 2)
    };
    const float delta_ap[] = {ap_x - pred_ap[0], ap_y - pred_ap[1]};
    sqr_ap_err = delta_ap[0]*delta_ap[0]+delta_ap[1]*delta_ap[1];
    float dx1_domega0[2][2];
    dx1_domega0[0][0] =  + 25.2552  + -1.83617 *lens_ipow(begin_dy, 2) + -4.56459 *lens_ipow(begin_dx, 2) + 0.147112 *begin_y*begin_dy + 0.00609101 *lens_ipow(begin_y, 2) + 0.274678 *begin_x*begin_dx + 0.0100807 *lens_ipow(begin_x, 2) + 0.874804 *lens_ipow(begin_lambda, 3) + 0.89088 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + 0.0354535 *begin_y*begin_dy*lens_ipow(begin_lambda, 3) + 0.00341868 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 4) + -2.09067e-05 *lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 3) + 4.35087e-06 *lens_ipow(begin_x, 4)*lens_ipow(begin_dx, 2) + 0.44816 *lens_ipow(begin_dy, 8) + 0.335404 *lens_ipow(begin_dx, 8) + 0.00482762 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 4) + 1.34488e-08 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 3)*begin_dx*begin_dy + -1.33148e-12 *lens_ipow(begin_x, 8) + -4.65791 *lens_ipow(begin_lambda, 10) + 3.50038 *lens_ipow(begin_dy, 4)*lens_ipow(begin_lambda, 6) + -0.00847409 *begin_x*begin_y*begin_dx*lens_ipow(begin_dy, 5)*lens_ipow(begin_lambda, 2) + 0.000710015 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dx, 6)*begin_dy+0.0f;
    dx1_domega0[0][1] =  + -3.67234 *begin_dx*begin_dy + 0.147112 *begin_y*begin_dx + -0.110856 *begin_x*begin_dy + 0.00494302 *begin_x*begin_y*begin_lambda + 0.59392 *lens_ipow(begin_dx, 3)*begin_dy + -7.38242e-05 *begin_x*lens_ipow(begin_y, 2)*begin_dy + 1.24453e-06 *lens_ipow(begin_x, 3)*begin_y + 0.0354535 *begin_y*begin_dx*lens_ipow(begin_lambda, 3) + -6.272e-05 *lens_ipow(begin_y, 3)*begin_dx*lens_ipow(begin_dy, 2) + 0.185275 *begin_x*begin_dy*lens_ipow(begin_lambda, 4) + 0.153912 *begin_x*lens_ipow(begin_dy, 5) + 3.58528 *begin_dx*lens_ipow(begin_dy, 7) + 6.7244e-09 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 2) + 1.20613e-07 *lens_ipow(begin_x, 5)*lens_ipow(begin_dy, 3) + 14.0015 *begin_dx*lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 6) + -0.0211852 *begin_x*begin_y*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 4)*lens_ipow(begin_lambda, 2) + 0.000101431 *lens_ipow(begin_x, 2)*begin_y*lens_ipow(begin_dx, 7)+0.0f;
    dx1_domega0[1][0] =  + -0.000487611 *begin_dy + -2.90525 *begin_dx*begin_dy + -0.0696762 *begin_y*begin_dx + 0.180011 *begin_x*begin_dy + 0.00367345 *begin_x*begin_y + -0.0176623 *begin_x*lens_ipow(begin_dy, 3) + 9.38447e-07 *begin_x*lens_ipow(begin_y, 3) + -3.07028e-05 *lens_ipow(begin_x, 2)*begin_y*begin_dx + 0.113155 *begin_y*lens_ipow(begin_dx, 5) + -0.000171818 *lens_ipow(begin_x, 3)*lens_ipow(begin_dx, 2)*begin_dy + 0.175835 *begin_dx*lens_ipow(begin_dy, 5)*begin_lambda + 1.4648 *lens_ipow(begin_dx, 7)*begin_dy + 8.76939e-08 *lens_ipow(begin_y, 5)*lens_ipow(begin_dx, 3) + -0.0353224 *begin_x*lens_ipow(begin_dx, 2)*begin_dy*lens_ipow(begin_lambda, 4) + 4.36515e-09 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 2) + -3.97786 *lens_ipow(begin_dx, 5)*lens_ipow(begin_dy, 5) + 9.05996e-05 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 7)+0.0f;
    dx1_domega0[1][1] =  + 25.2549  + -0.000487611 *begin_dx + -4.3886 *lens_ipow(begin_dy, 2) + -1.45263 *lens_ipow(begin_dx, 2) + 0.288511 *begin_y*begin_dy + 0.010328 *lens_ipow(begin_y, 2) + 0.180011 *begin_x*begin_dx + 0.00654957 *lens_ipow(begin_x, 2) + 0.797561 *lens_ipow(begin_lambda, 3) + -0.000196179 *begin_x*lens_ipow(begin_dy, 2) + -0.0529869 *begin_x*begin_dx*lens_ipow(begin_dy, 2) + 4.66035e-06 *lens_ipow(begin_y, 4)*lens_ipow(begin_dy, 2) + 0.00268865 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 4) + -5.72727e-05 *lens_ipow(begin_x, 3)*lens_ipow(begin_dx, 3) + 0.443449 *lens_ipow(begin_dy, 4)*lens_ipow(begin_lambda, 3) + 0.439588 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 4)*begin_lambda + 0.150615 *lens_ipow(begin_dy, 8) + 0.183099 *lens_ipow(begin_dx, 8) + -1.50573e-12 *lens_ipow(begin_y, 8) + -0.0117741 *begin_x*lens_ipow(begin_dx, 3)*lens_ipow(begin_lambda, 4) + 8.7303e-09 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 3)*begin_dx*begin_dy + -3.22089e-12 *lens_ipow(begin_x, 8) + -4.37325 *lens_ipow(begin_lambda, 10) + -3.31488 *lens_ipow(begin_dx, 6)*lens_ipow(begin_dy, 4) + 0.000634197 *begin_x*lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 6)+0.0f;
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
    out[0] =  + 2.01925e-05  + 30.421 *begin_dx + 4.18341e-05 *begin_y + 0.494773 *begin_x + 2.1368e-06 *begin_x*begin_y + 0.0035086 *lens_ipow(begin_dy, 3) + -12.3511 *begin_dx*lens_ipow(begin_dy, 2) + -13.7462 *lens_ipow(begin_dx, 3) + 0.00431615 *lens_ipow(begin_y, 2)*begin_dx + 0.00761577 *begin_x*begin_y*begin_dy + 5.13835e-05 *begin_x*lens_ipow(begin_y, 2) + 0.0144417 *lens_ipow(begin_x, 2)*begin_dx + 6.9513e-05 *lens_ipow(begin_x, 3) + 3.27288 *begin_dx*lens_ipow(begin_lambda, 4) + 17.3946 *lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 2) + 6.06409 *lens_ipow(begin_dx, 5) + 0.204007 *begin_y*lens_ipow(begin_dx, 3)*begin_dy + -0.00331509 *lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 2) + -0.0338336 *begin_x*lens_ipow(begin_dy, 4) + 0.268621 *begin_x*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + 0.000394869 *lens_ipow(begin_x, 3)*lens_ipow(begin_dx, 2) + 0.0964549 *begin_x*lens_ipow(begin_lambda, 5) + 5.74207 *begin_dx*lens_ipow(begin_dy, 6) + 7.95751e-05 *lens_ipow(begin_x, 3)*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 2) + -5.18691 *lens_ipow(begin_dx, 5)*lens_ipow(begin_dy, 2)*begin_lambda + -1.58279 *lens_ipow(begin_dx, 3)*lens_ipow(begin_dy, 6) + -0.749862 *lens_ipow(begin_dx, 9) + 2.75395e-07 *begin_x*lens_ipow(begin_y, 4)*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + 1.31971e-11 *lens_ipow(begin_x, 6)*lens_ipow(begin_y, 2)*begin_dx + 5.3637e-11 *lens_ipow(begin_x, 8)*begin_dx + -10.377 *begin_dx*lens_ipow(begin_lambda, 10) + 1.88946e-08 *lens_ipow(begin_y, 6)*lens_ipow(begin_dx, 5) + -0.0927918 *begin_x*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 8) + -0.00679733 *begin_x*begin_y*lens_ipow(begin_dy, 9) + 0.0102953 *begin_x*begin_y*lens_ipow(begin_dx, 4)*begin_dy*lens_ipow(begin_lambda, 4) + 3.69317e-08 *begin_x*lens_ipow(begin_y, 4)*lens_ipow(begin_lambda, 6) + -0.000436406 *lens_ipow(begin_x, 3)*lens_ipow(begin_lambda, 8) + -1.2961e-05 *lens_ipow(begin_x, 3)*begin_y*lens_ipow(begin_dx, 6)*begin_dy + 3.44419e-11 *lens_ipow(begin_x, 8)*begin_dx*lens_ipow(begin_dy, 2) + 1.7835e-12 *lens_ipow(begin_x, 9)*lens_ipow(begin_lambda, 2);
    out[1] =  + 0.000134552  + 30.4512 *begin_dy + 0.495423 *begin_y + -13.4389 *lens_ipow(begin_dy, 3) + -11.7804 *lens_ipow(begin_dx, 2)*begin_dy + 0.0139055 *lens_ipow(begin_y, 2)*begin_dy + 7.20459e-05 *lens_ipow(begin_y, 3) + 0.00775022 *begin_x*begin_y*begin_dx + 0.00390391 *lens_ipow(begin_x, 2)*begin_dy + 3.99703e-05 *lens_ipow(begin_x, 2)*begin_y + 2.65911 *begin_dy*lens_ipow(begin_lambda, 4) + 5.25799 *lens_ipow(begin_dy, 5) + -0.0243531 *begin_y*lens_ipow(begin_dx, 4) + 0.000282917 *lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 2) + -0.00974809 *begin_x*begin_y*begin_dx*lens_ipow(begin_dy, 2) + 16.5918 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 3)*begin_lambda + 0.0793546 *begin_y*lens_ipow(begin_lambda, 5) + -2.89862e-07 *lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 3) + -3.76053 *lens_ipow(begin_dx, 2)*begin_dy*lens_ipow(begin_lambda, 4) + 5.27996 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 5) + 4.60202 *lens_ipow(begin_dx, 6)*begin_dy + 0.000135874 *lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 2) + 0.000683088 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 5) + -1.38802e-05 *lens_ipow(begin_x, 3)*begin_y*begin_dx*lens_ipow(begin_lambda, 2) + -9.81338e-06 *lens_ipow(begin_x, 4)*lens_ipow(begin_dx, 2)*begin_dy + -0.280351 *lens_ipow(begin_dy, 9) + 10.9156 *lens_ipow(begin_dx, 6)*lens_ipow(begin_dy, 3) + -7.2789e-13 *lens_ipow(begin_y, 9) + -3.99257e-13 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 7) + -5.60589e-10 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 2) + -6.16018e-13 *lens_ipow(begin_x, 8)*begin_y + -5.10244 *begin_dy*lens_ipow(begin_lambda, 10) + -34.5504 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 5)*lens_ipow(begin_lambda, 4) + 0.367653 *begin_y*lens_ipow(begin_dx, 8)*lens_ipow(begin_dy, 2) + 0.0245688 *begin_x*begin_y*begin_dx*lens_ipow(begin_lambda, 8) + -0.035008 *begin_x*begin_y*lens_ipow(begin_dx, 5)*lens_ipow(begin_dy, 4) + -0.00375033 *begin_x*begin_y*lens_ipow(begin_dx, 9) + 0.000523625 *begin_x*lens_ipow(begin_y, 2)*begin_dx*lens_ipow(begin_dy, 7) + -2.70559e-09 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 4)*begin_dx*lens_ipow(begin_dy, 3) + -8.67083e-07 *lens_ipow(begin_x, 5)*begin_dx*begin_dy*lens_ipow(begin_lambda, 4);
    out[2] =  + -8.5188e-07  + 2.87302e-05 *begin_dy + -1.02651 *begin_dx + 1.19235e-06 *begin_y + -0.0492502 *begin_x + -0.902276 *begin_dx*lens_ipow(begin_dy, 2) + -0.553257 *lens_ipow(begin_dx, 3) + -0.0468609 *begin_y*begin_dx*begin_dy + -0.000567729 *lens_ipow(begin_y, 2)*begin_dx + -0.0215981 *begin_x*lens_ipow(begin_dy, 2) + -0.0411757 *begin_x*lens_ipow(begin_dx, 2) + -0.00103873 *begin_x*begin_y*begin_dy + -4.47595e-06 *begin_x*lens_ipow(begin_y, 2) + -0.00101613 *lens_ipow(begin_x, 2)*begin_dx + -5.62063e-07 *lens_ipow(begin_x, 3) + 0.113201 *begin_dx*lens_ipow(begin_lambda, 3) + 3.9992e-07 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2) + -0.000290235 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 3) + 0.0061882 *begin_x*lens_ipow(begin_lambda, 4) + 0.000480974 *begin_x*begin_y*lens_ipow(begin_dx, 2)*begin_dy + -0.000268528 *lens_ipow(begin_x, 2)*begin_dx*lens_ipow(begin_dy, 2) + 0.0630923 *lens_ipow(begin_dx, 5)*lens_ipow(begin_dy, 2) + 2.66538e-07 *begin_x*lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 3) + 1.23711e-05 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 5) + 9.80647e-09 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 3)*begin_dx*begin_dy + -9.6384e-09 *lens_ipow(begin_x, 5)*lens_ipow(begin_dy, 2) + 5.66917e-11 *lens_ipow(begin_y, 7)*begin_dx*begin_dy + -2.99423e-14 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 6) + -3.5794e-07 *lens_ipow(begin_x, 4)*begin_dx*lens_ipow(begin_dy, 4) + -2.94361e-11 *lens_ipow(begin_x, 5)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2) + 1.7055e-05 *lens_ipow(begin_y, 3)*begin_dx*begin_dy*lens_ipow(begin_lambda, 5) + -0.573908 *begin_dx*lens_ipow(begin_lambda, 10) + 0.0828084 *begin_dx*lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 8) + -0.0209417 *lens_ipow(begin_dx, 3)*lens_ipow(begin_lambda, 8) + -0.004286 *begin_y*begin_dx*lens_ipow(begin_dy, 9) + 3.14721e-12 *lens_ipow(begin_y, 8)*begin_dx*lens_ipow(begin_lambda, 2) + -0.0286649 *begin_x*lens_ipow(begin_lambda, 10) + 0.00390668 *begin_x*lens_ipow(begin_dy, 10) + -1.36586e-08 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 2)*lens_ipow(begin_lambda, 6) + -2.56995e-09 *lens_ipow(begin_x, 5)*begin_y*lens_ipow(begin_dx, 4)*begin_dy;
    out[3] =  + -9.46651e-06  + -1.02313 *begin_dy + -0.0492014 *begin_y + -0.570195 *lens_ipow(begin_dy, 3) + -0.203656 *lens_ipow(begin_dx, 2)*begin_dy + -0.0421999 *begin_y*lens_ipow(begin_dy, 2) + -3.10808e-06 *begin_y*begin_dx*begin_dy + -0.00967712 *begin_y*lens_ipow(begin_dx, 2) + -0.0010564 *lens_ipow(begin_y, 2)*begin_dy + -1.05146e-06 *lens_ipow(begin_y, 3) + -0.000276425 *begin_x*begin_y*begin_dx + -2.29079e-05 *lens_ipow(begin_x, 2)*begin_dy + 4.5662e-06 *lens_ipow(begin_x, 2)*begin_y + 0.115364 *begin_dy*lens_ipow(begin_lambda, 3) + -1.00689e-08 *lens_ipow(begin_y, 3)*begin_dy + 0.00626454 *begin_y*lens_ipow(begin_lambda, 4) + -0.0150238 *begin_x*begin_dx*lens_ipow(begin_dy, 3) + -0.000318517 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 3) + 0.000162642 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 2)*begin_dy + 1.35945e-07 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*begin_dy + -0.155475 *lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 3) + 0.00695589 *begin_y*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 4) + -6.36973e-07 *lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 4) + -0.00327368 *begin_x*lens_ipow(begin_dx, 5)*begin_dy + -2.20063e-07 *lens_ipow(begin_x, 4)*lens_ipow(begin_dy, 3) + -5.75323e-11 *lens_ipow(begin_x, 6)*begin_dy + 0.0138669 *lens_ipow(begin_dy, 9) + 0.00323231 *begin_y*lens_ipow(begin_dx, 8) + -1.4082e-10 *lens_ipow(begin_y, 6)*lens_ipow(begin_dx, 2)*begin_dy + -9.43962e-05 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 7) + 2.46735e-12 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 5)*lens_ipow(begin_dx, 2) + -0.594001 *begin_dy*lens_ipow(begin_lambda, 10) + -0.0881703 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 3)*lens_ipow(begin_lambda, 6) + -0.0293551 *begin_y*lens_ipow(begin_lambda, 10) + -0.00010206 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 8)*begin_dy + -1.36781e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 3)*lens_ipow(begin_lambda, 6) + -8.83314e-06 *lens_ipow(begin_x, 3)*begin_dx*lens_ipow(begin_dy, 7) + -2.99129e-14 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 5)*lens_ipow(begin_dy, 2) + -4.44069e-11 *lens_ipow(begin_x, 6)*begin_y*lens_ipow(begin_dx, 4) + -2.94901e-12 *lens_ipow(begin_x, 8)*lens_ipow(begin_dx, 2)*begin_dy;
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
    domega2_dx0[0][0] =  + -0.0492502  + -0.0215981 *lens_ipow(begin_dy, 2) + -0.0411757 *lens_ipow(begin_dx, 2) + -0.00103873 *begin_y*begin_dy + -4.47595e-06 *lens_ipow(begin_y, 2) + -0.00203225 *begin_x*begin_dx + -1.68619e-06 *lens_ipow(begin_x, 2) + 0.0061882 *lens_ipow(begin_lambda, 4) + 0.000480974 *begin_y*lens_ipow(begin_dx, 2)*begin_dy + -0.000537056 *begin_x*begin_dx*lens_ipow(begin_dy, 2) + 2.66538e-07 *lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 3) + 2.47423e-05 *begin_x*lens_ipow(begin_dx, 5) + 1.96129e-08 *begin_x*lens_ipow(begin_y, 3)*begin_dx*begin_dy + -4.8192e-08 *lens_ipow(begin_x, 4)*lens_ipow(begin_dy, 2) + -8.9827e-14 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 6) + -1.43176e-06 *lens_ipow(begin_x, 3)*begin_dx*lens_ipow(begin_dy, 4) + -1.47181e-10 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2) + -0.0286649 *lens_ipow(begin_lambda, 10) + 0.00390668 *lens_ipow(begin_dy, 10) + -4.09758e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*lens_ipow(begin_lambda, 6) + -1.28498e-08 *lens_ipow(begin_x, 4)*begin_y*lens_ipow(begin_dx, 4)*begin_dy+0.0f;
    domega2_dx0[0][1] =  + 1.19235e-06  + -0.0468609 *begin_dx*begin_dy + -0.00113546 *begin_y*begin_dx + -0.00103873 *begin_x*begin_dy + -8.95191e-06 *begin_x*begin_y + 7.99839e-07 *begin_y*lens_ipow(begin_dx, 2) + -0.000580469 *begin_y*lens_ipow(begin_dx, 3) + 0.000480974 *begin_x*lens_ipow(begin_dx, 2)*begin_dy + 7.99613e-07 *begin_x*lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 3) + 2.94194e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*begin_dx*begin_dy + 3.96842e-10 *lens_ipow(begin_y, 6)*begin_dx*begin_dy + -1.79654e-13 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 5) + -5.88723e-11 *lens_ipow(begin_x, 5)*begin_y*lens_ipow(begin_dx, 2) + 5.11651e-05 *lens_ipow(begin_y, 2)*begin_dx*begin_dy*lens_ipow(begin_lambda, 5) + -0.004286 *begin_dx*lens_ipow(begin_dy, 9) + 2.51777e-11 *lens_ipow(begin_y, 7)*begin_dx*lens_ipow(begin_lambda, 2) + -2.73172e-08 *lens_ipow(begin_x, 3)*begin_y*lens_ipow(begin_lambda, 6) + -2.56995e-09 *lens_ipow(begin_x, 5)*lens_ipow(begin_dx, 4)*begin_dy+0.0f;
    domega2_dx0[1][0] =  + -0.000276425 *begin_y*begin_dx + -4.58158e-05 *begin_x*begin_dy + 9.13239e-06 *begin_x*begin_y + -0.0150238 *begin_dx*lens_ipow(begin_dy, 3) + -0.000637033 *begin_x*lens_ipow(begin_dy, 3) + 0.000325284 *begin_x*lens_ipow(begin_dx, 2)*begin_dy + 2.71891e-07 *begin_x*lens_ipow(begin_y, 2)*begin_dy + -0.00327368 *lens_ipow(begin_dx, 5)*begin_dy + -8.80251e-07 *lens_ipow(begin_x, 3)*lens_ipow(begin_dy, 3) + -3.45194e-10 *lens_ipow(begin_x, 5)*begin_dy + -0.000188792 *begin_x*lens_ipow(begin_dy, 7) + 4.9347e-12 *begin_x*lens_ipow(begin_y, 5)*lens_ipow(begin_dx, 2) + -2.73561e-08 *begin_x*lens_ipow(begin_y, 3)*lens_ipow(begin_lambda, 6) + -2.64994e-05 *lens_ipow(begin_x, 2)*begin_dx*lens_ipow(begin_dy, 7) + -1.19652e-13 *lens_ipow(begin_x, 3)*lens_ipow(begin_y, 5)*lens_ipow(begin_dy, 2) + -2.66441e-10 *lens_ipow(begin_x, 5)*begin_y*lens_ipow(begin_dx, 4) + -2.35921e-11 *lens_ipow(begin_x, 7)*lens_ipow(begin_dx, 2)*begin_dy+0.0f;
    domega2_dx0[1][1] =  + -0.0492014  + -0.0421999 *lens_ipow(begin_dy, 2) + -3.10808e-06 *begin_dx*begin_dy + -0.00967712 *lens_ipow(begin_dx, 2) + -0.00211281 *begin_y*begin_dy + -3.15439e-06 *lens_ipow(begin_y, 2) + -0.000276425 *begin_x*begin_dx + 4.5662e-06 *lens_ipow(begin_x, 2) + -3.02066e-08 *lens_ipow(begin_y, 2)*begin_dy + 0.00626454 *lens_ipow(begin_lambda, 4) + 2.71891e-07 *lens_ipow(begin_x, 2)*begin_y*begin_dy + 0.00695589 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 4) + -1.91092e-06 *lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 4) + 0.00323231 *lens_ipow(begin_dx, 8) + -8.44922e-10 *lens_ipow(begin_y, 5)*lens_ipow(begin_dx, 2)*begin_dy + 1.23368e-11 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 4)*lens_ipow(begin_dx, 2) + -0.0293551 *lens_ipow(begin_lambda, 10) + -0.00020412 *begin_y*lens_ipow(begin_dx, 8)*begin_dy + -4.10342e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*lens_ipow(begin_lambda, 6) + -1.49564e-13 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 4)*lens_ipow(begin_dy, 2) + -4.44069e-11 *lens_ipow(begin_x, 6)*lens_ipow(begin_dx, 4)+0.0f;
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
  out[4] =  + 0.810678  + 0.0952603 *begin_lambda + -9.98441e-06 *begin_dy + -4.54686e-07 *begin_x + -0.0368051 *lens_ipow(begin_dy, 2) + -0.0368304 *lens_ipow(begin_dx, 2) + -0.00152821 *begin_y*begin_dy + 1.15717e-06 *begin_y*begin_dx + -2.61916e-05 *lens_ipow(begin_y, 2) + -9.82938e-07 *begin_x*begin_dy + -0.00153956 *begin_x*begin_dx + -2.85621e-05 *lens_ipow(begin_x, 2) + -0.0798538 *lens_ipow(begin_lambda, 3) + -0.101475 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 2) + 3.87742e-05 *begin_x*begin_y*begin_dx*begin_dy + 1.32615e-05 *lens_ipow(begin_x, 2)*lens_ipow(begin_dy, 2) + 2.43746e-05 *lens_ipow(begin_y, 2)*lens_ipow(begin_dy, 2)*begin_lambda + 2.63961e-05 *lens_ipow(begin_y, 2)*lens_ipow(begin_dx, 2)*begin_lambda + -1.5318e-06 *lens_ipow(begin_x, 3)*begin_dx*begin_lambda + -0.0967881 *lens_ipow(begin_dy, 6) + -0.0761375 *lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 4) + -0.0444828 *lens_ipow(begin_dx, 4)*lens_ipow(begin_dy, 2) + -0.0924571 *lens_ipow(begin_dx, 6) + -1.82861e-10 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 2) + -0.00957324 *lens_ipow(begin_dy, 2)*lens_ipow(begin_lambda, 5) + -0.0133643 *lens_ipow(begin_dx, 2)*lens_ipow(begin_lambda, 5) + -1.31181e-07 *lens_ipow(begin_y, 4)*lens_ipow(begin_dx, 4) + -5.52301e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 2)*lens_ipow(begin_lambda, 5) + -0.00319302 *begin_y*lens_ipow(begin_dy, 9) + -0.00631391 *begin_y*lens_ipow(begin_dx, 2)*lens_ipow(begin_dy, 7) + -0.00179212 *begin_y*lens_ipow(begin_dx, 8)*begin_dy + -5.04147e-16 *lens_ipow(begin_y, 10) + -0.00306392 *begin_x*lens_ipow(begin_dx, 5)*lens_ipow(begin_dy, 4) + 0.000165949 *lens_ipow(begin_x, 2)*lens_ipow(begin_dx, 8) + 1.36554e-08 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 3)*lens_ipow(begin_dy, 5) + 1.63207e-16 *lens_ipow(begin_x, 2)*lens_ipow(begin_y, 8) + 3.41374e-11 *lens_ipow(begin_x, 4)*lens_ipow(begin_y, 3)*lens_ipow(begin_dx, 2)*begin_dy + 2.90361e-12 *lens_ipow(begin_x, 6)*lens_ipow(begin_dy, 4) + -7.11866e-13 *lens_ipow(begin_x, 8)*lens_ipow(begin_dx, 2) + 0.159353 *lens_ipow(begin_lambda, 11);
else
  out[4] = 0.0f;