const float out_x =  + 0.000146105  + 119.392 *dx + -4.05857e-06 *y + 1.247 *x + -0.000155495 *y*dx + -7.6004e-07 *lens_ipow(y, 2) + -0.0470823 *x*lambda + -0.000269889 *x*dy + 0.000205308 *x*dx + -2.41965e-06 *x*y + 315.426 *dx*lens_ipow(dy, 2) + 312.895 *lens_ipow(dx, 3) + 5.80165 *y*dx*dy + 0.025389 *lens_ipow(y, 2)*dx + 2.95132 *x*lens_ipow(dy, 2) + 8.54125 *x*lens_ipow(dx, 2) + 0.0525601 *x*y*dy + 0.00021685 *x*lens_ipow(y, 2) + 0.0772709 *lens_ipow(x, 2)*dx + 0.000215934 *lens_ipow(x, 3) + -0.0139484 *x*dx*lens_ipow(dy, 2) + 0.379694 *x*lens_ipow(dx, 2)*lambda + -0.000975206 *x*y*dy*lambda + 0.00264109 *x*lens_ipow(lambda, 4) + 0.655091 *x*lens_ipow(dy, 4) + -0.0239011 *x*y*lens_ipow(dx, 2)*dy + -10.6975 *dx*lens_ipow(lambda, 5) + 5784.08 *lens_ipow(dx, 5)*lens_ipow(dy, 2) + 0.003015 *x*dy*lens_ipow(lambda, 5) + -5.51715e-05 *lens_ipow(x, 4)*lens_ipow(dx, 3) + 0.0233978 *y*dx*dy*lens_ipow(lambda, 5) + 1115.85 *lens_ipow(dx, 5)*lens_ipow(lambda, 4) + 7923.29 *lens_ipow(dx, 9) + 0.0392799 *x*lens_ipow(dy, 2)*lens_ipow(lambda, 6) + -2.15026e-05 *x*lens_ipow(y, 2)*lens_ipow(lambda, 6) + -6.55007e-14 *lens_ipow(x, 5)*lens_ipow(y, 4) + 1.20246e-08 *lens_ipow(x, 6)*dx*lens_ipow(lambda, 3) + 42.4803 *dx*lens_ipow(lambda, 10) + 403790 *lens_ipow(dx, 3)*lens_ipow(dy, 6)*lens_ipow(lambda, 2) + 0.141664 *x*lens_ipow(lambda, 10);
const float out_y =  + 0.000109425  + 119.438 *dy + 0.00155992 *dx + 1.23169 *y + 1.63555e-05 *x + 0.00746521 *lens_ipow(dx, 2) + 0.000168035 *y*dy + -0.000243009 *y*dx + -1.17086e-06 *lens_ipow(x, 2) + 314.444 *lens_ipow(dy, 3) + 312.195 *lens_ipow(dx, 2)*dy + 8.78967 *y*lens_ipow(dy, 2) + 0.00178997 *y*dx*dy + 2.98505 *y*lens_ipow(dx, 2) + 0.0781152 *lens_ipow(y, 2)*dy + 0.00022001 *lens_ipow(y, 3) + 5.76403 *x*dx*dy + 0.0524156 *x*y*dx + -1.65767e-07 *x*lens_ipow(y, 2) + 0.0251438 *lens_ipow(x, 2)*dy + 0.000217566 *lens_ipow(x, 2)*y + -0.0660435 *y*lens_ipow(lambda, 3) + -1.05783e-08 *lens_ipow(x, 3)*y + -0.000238607 *lens_ipow(y, 3)*lens_ipow(dx, 2) + 0.000734854 *x*lens_ipow(y, 2)*dx*dy + 3.20796e-06 *lens_ipow(x, 2)*lens_ipow(y, 2)*dy + -12.3294 *dy*lens_ipow(lambda, 5) + 0.0040916 *lens_ipow(x, 2)*lens_ipow(dy, 4) + 11026.7 *lens_ipow(dx, 2)*lens_ipow(dy, 5) + -4.69828e-05 *lens_ipow(y, 4)*lens_ipow(dy, 3) + -0.00574666 *lens_ipow(y, 2)*dy*lens_ipow(lambda, 5) + -4.14292e-08 *lens_ipow(y, 5)*lens_ipow(lambda, 3) + 42394.4 *lens_ipow(dy, 9) + 1162.52 *lens_ipow(dx, 4)*dy*lens_ipow(lambda, 4) + 9.79145 *y*lens_ipow(dx, 4)*lens_ipow(lambda, 4) + 919.141 *lens_ipow(dy, 5)*lens_ipow(lambda, 5) + 51.7455 *dy*lens_ipow(lambda, 10) + 0.348236 *y*lens_ipow(lambda, 10) + -0.0370741 *lens_ipow(x, 4)*lens_ipow(dy, 7) + -9.08606e-17 *lens_ipow(x, 4)*lens_ipow(y, 7);
const float out_dx =  + 1.87894e-07  + 0.872308 *dx + 0.000640927 *x + -5.45618 *dx*lens_ipow(dy, 2) + -5.40499 *lens_ipow(dx, 3) + -0.073883 *y*dx*dy + -0.000220612 *lens_ipow(y, 2)*dx + -0.0359731 *x*lens_ipow(dy, 2) + -0.109628 *x*lens_ipow(dx, 2) + -0.000432771 *x*y*dy + -9.23302e-07 *x*lens_ipow(y, 2) + -0.000653513 *lens_ipow(x, 2)*dx + -9.22903e-07 *lens_ipow(x, 3) + 0.0892395 *dx*lens_ipow(lambda, 3) + -1.55338 *lens_ipow(dx, 5) + 0.00077113 *x*lens_ipow(lambda, 4) + 3.68011e-11 *lens_ipow(x, 3)*lens_ipow(y, 2) + 7.69945e-06 *lens_ipow(x, 4)*lens_ipow(dx, 5) + -0.485149 *dx*lens_ipow(lambda, 10) + -0.00382843 *x*lens_ipow(lambda, 10);
const float out_dy =  + 7.38203e-08  + 0.872227 *dy + 0.000639659 *y + -5.40825 *lens_ipow(dy, 3) + -5.46406 *lens_ipow(dx, 2)*dy + -0.109758 *y*lens_ipow(dy, 2) + -0.0359932 *y*lens_ipow(dx, 2) + -0.000655186 *lens_ipow(y, 2)*dy + -9.28129e-07 *lens_ipow(y, 3) + -0.0739962 *x*dx*dy + -0.000432971 *x*y*dx + -0.000220468 *lens_ipow(x, 2)*dy + -9.12588e-07 *lens_ipow(x, 2)*y + 0.0901749 *dy*lens_ipow(lambda, 3) + -1.64243 *lens_ipow(dy, 5) + 0.000790711 *y*lens_ipow(lambda, 4) + 8.42114e-06 *lens_ipow(y, 4)*lens_ipow(dy, 5) + -0.496129 *dy*lens_ipow(lambda, 10) + -0.00400237 *y*lens_ipow(lambda, 10);
const float out_transmittance =  + 0.723565  + 0.135528 *lambda + 9.6303e-06 *dy + -4.86495e-06 *dx + 2.99687e-07 *y + -8.19086e-08 *x + -0.105107 *lens_ipow(dy, 2) + -0.161621 *lens_ipow(dx, 2) + -0.00214702 *y*dy + -1.13636e-05 *lens_ipow(y, 2) + -0.00224038 *x*dx + 3.2561e-09 *x*y + -1.1997e-05 *lens_ipow(x, 2) + -0.11374 *lens_ipow(lambda, 3) + 5.21372e-05 *y*lens_ipow(dx, 2) + -3.62746e-09 *lens_ipow(x, 2)*y + 23.6446 *lens_ipow(dx, 2)*lens_ipow(dy, 2) + 5.73641 *lens_ipow(dx, 4) + 0.192988 *y*lens_ipow(dx, 2)*dy + -0.00196333 *lens_ipow(y, 2)*lens_ipow(dy, 2) + 3.9132e-05 *lens_ipow(y, 2)*lens_ipow(dx, 2) + -1.80697e-05 *lens_ipow(y, 3)*dy + 0.186359 *x*dx*lens_ipow(dy, 2) + -1.20639e-05 *x*lens_ipow(y, 2)*dx + -8.51877e-07 *lens_ipow(x, 2)*lens_ipow(dy, 2) + -0.00199798 *lens_ipow(x, 2)*lens_ipow(dx, 2) + -1.22748e-05 *lens_ipow(x, 2)*y*dy + -9.97164e-08 *lens_ipow(x, 2)*lens_ipow(y, 2) + -1.92459e-05 *lens_ipow(x, 3)*dx + -0.0247724 *lens_ipow(dy, 2)*lens_ipow(lambda, 4) + -2.55308 *y*lens_ipow(dy, 5) + 1.74833e-06 *lens_ipow(y, 2)*lens_ipow(lambda, 4) + -1.91172e-10 *lens_ipow(y, 6) + -1.78277e-10 *lens_ipow(x, 6) + -27242 *lens_ipow(dy, 10) + -1.51599e-07 *lens_ipow(y, 6)*lens_ipow(dy, 4) + -1.59765e-13 *lens_ipow(y, 9)*dy + 2.27518e-06 *lens_ipow(x, 5)*lens_ipow(dx, 5) + -9.49538e-14 *lens_ipow(x, 9)*dx + 0.206428 *lens_ipow(lambda, 11);