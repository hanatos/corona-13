#ifndef _CORONA_MAT_3F_H
#define _CORONA_MAT_3F_H

#ifdef __cplusplus
#define restrict
#endif

#include <string.h>

#define mreal float
#define funcname(name) mat3_ ## name
#include "matrix3.inc"
#undef funcname
#undef mreal

#ifdef __cplusplus
#undef restrict
#endif
#endif
