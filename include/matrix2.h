#ifndef _CORONA_MAT_2F_H
#define _CORONA_MAT_2F_H

#define mreal float
#define funcname(name) mat2_ ## name
#include "matrix2.inc"
#undef funcname
#undef mreal

#endif

