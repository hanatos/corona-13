
static inline void funcname(set_zero)(mreal *m)
{
  m[0] = m[1] = m[2] = m[3] = 0.0f;
}

static inline void funcname(set_identity)(mreal *m)
{
  m[0] = 1.0f; m[1] = 0.0f; m[2] = 0.0f; m[3] = 1.0f;
}

static inline void funcname(transpose)(const mreal *a, mreal *res)
{
  res[0] = a[0];
  res[1] = a[2];
  res[2] = a[1];
  res[3] = a[3];
}

static inline void funcname(add)(const mreal *a, const mreal *b, mreal *res)
{
  for(int k=0;k<4;k++) res[k] = a[k] + b[k];
}

static inline void funcname(sub)(const mreal *a, const mreal *b, mreal *res)
{
  for(int k=0;k<4;k++) res[k] = a[k] - b[k];
}

static inline void funcname(mul)(const mreal *a, const mreal *b, mreal *res)
{
  funcname(set_zero)(res);
  for(int j=0;j<2;j++)
    for(int i=0;i<2;i++)
      for(int k=0;k<2;k++)
        res[i+2*j] += a[2*j+k] * b[2*k+i];
}

static inline void funcname(mulv)(const mreal *a, const mreal *v, mreal *res)
{
  res[0] = res[1] = 0.0f;
  for(int j=0;j<2;j++)
    for(int i=0;i<2;i++)
        res[j] += a[2*j+i] * v[i];
}

static inline mreal funcname(det)(const mreal *a)
{
  return a[0]*a[3] - a[1]*a[2];
}

static inline mreal funcname(invert)(const mreal *a, mreal *inv)
{
  const mreal det = funcname(det)(a);
  if(!(det != 0.0f)) return 0.0f;
  inv[0] = a[3]/det;
  inv[1] = -a[1]/det;
  inv[2] = -a[2]/det;
  inv[3] = a[0]/det;
  return det;
}

