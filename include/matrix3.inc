
static inline void funcname(set_zero)(mreal *restrict m)
{
  memset(m, 0, sizeof(mreal)*9);
}

static inline void funcname(set_identity)(mreal *restrict m)
{
  funcname(set_zero)(m);
  m[0] = m[4] = m[8] = 1.0;
}

static inline void funcname(add)(
    const mreal *const restrict a,
    const mreal *const restrict b,
    mreal *const restrict res)
{
  for(int k=0;k<9;k++) res[k] = a[k] + b[k];
}

static inline void funcname(sub)(
    const mreal *const restrict a,
    const mreal *const restrict b,
    mreal *const restrict res)
{
  for(int k=0;k<9;k++) res[k] = a[k] - b[k];
}

static inline void funcname(mul)(
    const mreal *const restrict a,
    const mreal *const restrict b,
    mreal *const restrict res)
{
  funcname(set_zero)(res);
  for(int j=0;j<3;j++)
    for(int i=0;i<3;i++)
      for(int k=0;k<3;k++)
        res[i+3*j] += a[3*j+k] * b[3*k+i];
}

static inline void funcname(mulv)(
    const mreal *const restrict a,
    const mreal *const restrict v,
    mreal *const restrict res)
{
  res[0] = res[1] = res[2] = 0.0f;
  for(int j=0;j<3;j++)
    for(int i=0;i<3;i++)
        res[j] += a[3*j+i] * v[i];
}

static inline void funcname(transpose)(const mreal *const restrict m, mreal *const restrict res)
{
  for(int i=0;i<3;i++) for(int k=0;k<3;k++) res[3*k+i] = m[3*i+k];
}

static inline void funcname(transpose_sub2)(const mreal *const restrict m, mreal *const restrict res)
{
  for(int i=0;i<3;i++) for(int k=0;k<3;k++){
    if (i < 2 || k < 2) 
      res[3*k+i] = m[3*i+k];
    else 
      res[3*k+i] = m[3*k+i];
  }
}

static inline mreal funcname(det)(const mreal *const restrict a)
{
#define A(y, x) a[(y - 1) * 3 + (x - 1)]
  return
    A(1, 1) * (A(3, 3) * A(2, 2) - A(3, 2) * A(2, 3)) -
    A(2, 1) * (A(3, 3) * A(1, 2) - A(3, 2) * A(1, 3)) +
    A(3, 1) * (A(2, 3) * A(1, 2) - A(2, 2) * A(1, 3));
}

static inline mreal funcname(invert)(const mreal *const restrict a, mreal *const restrict inv)
{
  const mreal det = funcname(det)(a);
  if(!(det != (mreal)0.0)) return (mreal)0.0;
  const mreal invdet = (mreal)1.0 / det;

  inv[3*0+0] =  invdet * (A(3, 3) * A(2, 2) - A(3, 2) * A(2, 3));
  inv[3*0+1] = -invdet * (A(3, 3) * A(1, 2) - A(3, 2) * A(1, 3));
  inv[3*0+2] =  invdet * (A(2, 3) * A(1, 2) - A(2, 2) * A(1, 3));

  inv[3*1+0] = -invdet * (A(3, 3) * A(2, 1) - A(3, 1) * A(2, 3));
  inv[3*1+1] =  invdet * (A(3, 3) * A(1, 1) - A(3, 1) * A(1, 3));
  inv[3*1+2] = -invdet * (A(2, 3) * A(1, 1) - A(2, 1) * A(1, 3));

  inv[3*2+0] =  invdet * (A(3, 2) * A(2, 1) - A(3, 1) * A(2, 2));
  inv[3*2+1] = -invdet * (A(3, 2) * A(1, 1) - A(3, 1) * A(1, 2));
  inv[3*2+2] =  invdet * (A(2, 2) * A(1, 1) - A(2, 1) * A(1, 2));
  return det;
#undef A
}

static inline mreal funcname(invert_sub2)(const mreal *const restrict a, mreal *const restrict inv)
{
  const mreal det = a[0]*a[4] - a[1]*a[3];
  if(!(det != (mreal)0.0)) return (mreal)0.0;
  inv[0] = a[4]/det;
  inv[1] = -a[1]/det;
  inv[3] = -a[3]/det;
  inv[4] = a[0]/det;

  inv[2] = inv[5] = inv[6] = inv[7] = (mreal)0.0;
  inv[8] = (mreal)1.0;
  return det;
}

static inline void funcname(rotate)(
	const mreal *const axis, 
	mreal angle, // in degree
	mreal *const restrict res)
{
#define A(y, x) res[y * 3 + x]
	angle = angle/((mreal)180)*M_PI;
	mreal sin_theta = (mreal)sinf(angle);
	mreal cos_theta = (mreal)cosf(angle);

	A(0, 0) = axis[0]*axis[0] + (1 - axis[0]*axis[0]) * cos_theta;
	A(0, 1) = axis[0]*axis[1] * (1 - cos_theta) - axis[2] * sin_theta;
	A(0, 2) = axis[0]*axis[2] * (1 - cos_theta) + axis[1] * sin_theta;

	A(1, 0) = axis[1]*axis[0] * (1 - cos_theta) + axis[2] * sin_theta;
	A(1, 1) = axis[1]*axis[1] + (1 - axis[1]*axis[1]) * cos_theta;
	A(1, 2) = axis[1]*axis[2] * (1 - cos_theta) - axis[0] * sin_theta;
	
	A(2, 0) = axis[2]*axis[0] * (1 - cos_theta) - axis[1] * sin_theta;
	A(2, 1) = axis[2]*axis[1] * (1 - cos_theta) + axis[0] * sin_theta;
	A(2, 2) = axis[2]*axis[2] + (1 - axis[2]*axis[2]) * cos_theta;
#undef A
}
