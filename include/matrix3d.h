#ifndef _CORONA_MAT_3D_H
#define _CORONA_MAT_3D_H

#ifdef __cplusplus
#define restrict
#endif

#include <string.h>

#define mreal double
#define funcname(name) mat3d_ ## name
#include "matrix3.inc"
#undef funcname
#undef mreal

#ifdef __cplusplus
#undef restrict
#endif
#endif
