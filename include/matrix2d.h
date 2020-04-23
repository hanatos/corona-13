#ifndef _CORONA_MAT_2D_H
#define _CORONA_MAT_2D_H

#define mreal double
#define funcname(name) mat2d_ ## name
#include "matrix2.inc"
#undef funcname
#undef mreal

#endif
