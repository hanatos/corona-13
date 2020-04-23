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
  pred_x =  + 30 *begin_dx + 1 *begin_x;
  pred_y =  + 30 *begin_dy + 1 *begin_y;
  pred_dx =  + 1 *begin_dx;
  pred_dy =  + 1 *begin_dy;
  float dx1_domega0[2][2];
  dx1_domega0[0][0] =  + 30 +0.0f;
  dx1_domega0[0][1] = +0.0f;
  dx1_domega0[1][0] = +0.0f;
  dx1_domega0[1][1] =  + 30 +0.0f;
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
