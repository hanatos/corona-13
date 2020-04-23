const float dx00 =  + 0.88529  + 5.31843e-07 *y + -0.0554282 *lens_ipow(dy, 2) + 0.137339 *lens_ipow(dx, 2) + 7.14484e-05 *lens_ipow(y, 2) + 0.0201614 *x*dx + 0.000322765 *lens_ipow(x, 2) + 0.00494302 *y*dy*lambda + 0.0355258 *lens_ipow(lambda, 4) + -3.69121e-05 *lens_ipow(y, 2)*lens_ipow(dy, 2) + 3.7336e-06 *lens_ipow(x, 2)*y*dy + 0.000170763 *lens_ipow(y, 2)*lens_ipow(lambda, 3) + 0.0926377 *lens_ipow(dy, 2)*lens_ipow(lambda, 4) + 0.0256521 *lens_ipow(dy, 6) + 5.80116e-06 *lens_ipow(x, 3)*lens_ipow(dx, 3) + 2.01732e-08 *lens_ipow(x, 2)*lens_ipow(y, 3)*lens_ipow(dx, 2)*dy + 1.50767e-07 *lens_ipow(x, 4)*lens_ipow(dy, 4) + -1.06518e-11 *lens_ipow(x, 7)*dx + -0.189056 *lens_ipow(lambda, 10) + -0.00423704 *y*lens_ipow(dx, 2)*lens_ipow(dy, 5)*lens_ipow(lambda, 2) + -0.000459711 *lens_ipow(y, 2)*lens_ipow(lambda, 8) + 0.000202861 *x*y*lens_ipow(dx, 7)*dy + 2.04183e-12 *lens_ipow(x, 2)*lens_ipow(y, 6)*lens_ipow(lambda, 2)+0.0f;
const float dx01 =  + 1.64108e-06  + 5.31843e-07 *x + 0.147112 *dx*dy + 0.012182 *y*dx + 0.000142897 *x*y + 0.00494302 *x*dy*lambda + -7.38242e-05 *x*y*lens_ipow(dy, 2) + 1.24453e-06 *lens_ipow(x, 3)*dy + 0.0354535 *dx*dy*lens_ipow(lambda, 3) + 0.000341527 *x*y*lens_ipow(lambda, 3) + 0.00136747 *y*lens_ipow(dx, 5) + -6.272e-05 *lens_ipow(y, 2)*dx*lens_ipow(dy, 3) + 0.00321842 *y*lens_ipow(dx, 3)*lens_ipow(lambda, 4) + 2.01732e-08 *lens_ipow(x, 3)*lens_ipow(y, 2)*lens_ipow(dx, 2)*dy + -0.00423704 *x*lens_ipow(dx, 2)*lens_ipow(dy, 5)*lens_ipow(lambda, 2) + -0.000919421 *x*y*lens_ipow(lambda, 8) + 0.000101431 *lens_ipow(x, 2)*lens_ipow(dx, 7)*dy + 4.08366e-12 *lens_ipow(x, 3)*lens_ipow(y, 5)*lens_ipow(lambda, 2)+0.0f;
const float dx02 =  + 25.2552  + -1.83617 *lens_ipow(dy, 2) + -4.56459 *lens_ipow(dx, 2) + 0.147112 *y*dy + 0.00609101 *lens_ipow(y, 2) + 0.274678 *x*dx + 0.0100807 *lens_ipow(x, 2) + 0.874804 *lens_ipow(lambda, 3) + 0.89088 *lens_ipow(dx, 2)*lens_ipow(dy, 2) + 0.0354535 *y*dy*lens_ipow(lambda, 3) + 0.00341868 *lens_ipow(y, 2)*lens_ipow(dx, 4) + -2.09067e-05 *lens_ipow(y, 3)*lens_ipow(dy, 3) + 4.35087e-06 *lens_ipow(x, 4)*lens_ipow(dx, 2) + 0.44816 *lens_ipow(dy, 8) + 0.335404 *lens_ipow(dx, 8) + 0.00482762 *lens_ipow(y, 2)*lens_ipow(dx, 2)*lens_ipow(lambda, 4) + 1.34488e-08 *lens_ipow(x, 3)*lens_ipow(y, 3)*dx*dy + -1.33148e-12 *lens_ipow(x, 8) + -4.65791 *lens_ipow(lambda, 10) + 3.50038 *lens_ipow(dy, 4)*lens_ipow(lambda, 6) + -0.00847409 *x*y*dx*lens_ipow(dy, 5)*lens_ipow(lambda, 2) + 0.000710015 *lens_ipow(x, 2)*y*lens_ipow(dx, 6)*dy+0.0f;
const float dx03 =  + -3.67234 *dx*dy + 0.147112 *y*dx + -0.110856 *x*dy + 0.00494302 *x*y*lambda + 0.59392 *lens_ipow(dx, 3)*dy + -7.38242e-05 *x*lens_ipow(y, 2)*dy + 1.24453e-06 *lens_ipow(x, 3)*y + 0.0354535 *y*dx*lens_ipow(lambda, 3) + -6.272e-05 *lens_ipow(y, 3)*dx*lens_ipow(dy, 2) + 0.185275 *x*dy*lens_ipow(lambda, 4) + 0.153912 *x*lens_ipow(dy, 5) + 3.58528 *dx*lens_ipow(dy, 7) + 6.7244e-09 *lens_ipow(x, 3)*lens_ipow(y, 3)*lens_ipow(dx, 2) + 1.20613e-07 *lens_ipow(x, 5)*lens_ipow(dy, 3) + 14.0015 *dx*lens_ipow(dy, 3)*lens_ipow(lambda, 6) + -0.0211852 *x*y*lens_ipow(dx, 2)*lens_ipow(dy, 4)*lens_ipow(lambda, 2) + 0.000101431 *lens_ipow(x, 2)*y*lens_ipow(dx, 7)+0.0f;
const float dx04 =  + 2.62441 *dx*lens_ipow(lambda, 2) + 0.00494302 *x*y*dy + 0.142103 *x*lens_ipow(lambda, 3) + 0.10636 *y*dx*dy*lens_ipow(lambda, 2) + 0.00051229 *x*lens_ipow(y, 2)*lens_ipow(lambda, 2) + 0.370551 *x*lens_ipow(dy, 2)*lens_ipow(lambda, 3) + 0.00643683 *lens_ipow(y, 2)*lens_ipow(dx, 3)*lens_ipow(lambda, 3) + -46.5791 *dx*lens_ipow(lambda, 9) + 21.0023 *dx*lens_ipow(dy, 4)*lens_ipow(lambda, 5) + -1.89056 *x*lens_ipow(lambda, 9) + -0.00847409 *x*y*lens_ipow(dx, 2)*lens_ipow(dy, 5)*lambda + -0.00367769 *x*lens_ipow(y, 2)*lens_ipow(lambda, 7) + 1.36122e-12 *lens_ipow(x, 3)*lens_ipow(y, 6)*lambda+0.0f;
const float dx10 =  + 2.46235e-06  + 0.180011 *dx*dy + 0.00367345 *y*dx + 0.0130991 *x*dy + 0.000217952 *x*y + -6.53929e-05 *lens_ipow(dy, 3) + -0.0176623 *dx*lens_ipow(dy, 3) + 9.38447e-07 *lens_ipow(y, 3)*dx + -3.07028e-05 *x*y*lens_ipow(dx, 2) + 0.00107546 *x*lens_ipow(dy, 5) + -0.000171818 *lens_ipow(x, 2)*lens_ipow(dx, 3)*dy + -0.0117741 *lens_ipow(dx, 3)*dy*lens_ipow(lambda, 4) + 1.63265e-07 *x*lens_ipow(y, 3)*lens_ipow(lambda, 4) + 1.30954e-08 *lens_ipow(x, 2)*lens_ipow(y, 3)*dx*lens_ipow(dy, 2) + -2.57672e-11 *lens_ipow(x, 7)*dy + 9.05996e-05 *lens_ipow(y, 2)*dx*lens_ipow(dy, 7) + -7.65402e-14 *lens_ipow(x, 5)*lens_ipow(y, 3)*lens_ipow(lambda, 2)+0.0f;
const float dx11 =  + 0.884941  + 0.144255 *lens_ipow(dy, 2) + -0.0348381 *lens_ipow(dx, 2) + 0.0206559 *y*dy + 0.000328663 *lens_ipow(y, 2) + 0.00367345 *x*dx + 0.000108976 *lens_ipow(x, 2) + 0.0362703 *lens_ipow(lambda, 4) + 2.81534e-06 *x*lens_ipow(y, 2)*dx + -1.53514e-05 *lens_ipow(x, 2)*lens_ipow(dx, 2) + 0.0188592 *lens_ipow(dx, 6) + 6.2138e-06 *lens_ipow(y, 3)*lens_ipow(dy, 3) + 1.09617e-07 *lens_ipow(y, 4)*lens_ipow(dx, 4) + -1.20459e-11 *lens_ipow(y, 7)*dy + 2.44898e-07 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(lambda, 4) + 1.30954e-08 *lens_ipow(x, 3)*lens_ipow(y, 2)*dx*lens_ipow(dy, 2) + -0.183804 *lens_ipow(lambda, 10) + 0.000181199 *x*y*dx*lens_ipow(dy, 7) + -3.82701e-14 *lens_ipow(x, 6)*lens_ipow(y, 2)*lens_ipow(lambda, 2)+0.0f;
const float dx12 =  + -0.000487611 *dy + -2.90525 *dx*dy + -0.0696762 *y*dx + 0.180011 *x*dy + 0.00367345 *x*y + -0.0176623 *x*lens_ipow(dy, 3) + 9.38447e-07 *x*lens_ipow(y, 3) + -3.07028e-05 *lens_ipow(x, 2)*y*dx + 0.113155 *y*lens_ipow(dx, 5) + -0.000171818 *lens_ipow(x, 3)*lens_ipow(dx, 2)*dy + 0.175835 *dx*lens_ipow(dy, 5)*lambda + 1.4648 *lens_ipow(dx, 7)*dy + 8.76939e-08 *lens_ipow(y, 5)*lens_ipow(dx, 3) + -0.0353224 *x*lens_ipow(dx, 2)*dy*lens_ipow(lambda, 4) + 4.36515e-09 *lens_ipow(x, 3)*lens_ipow(y, 3)*lens_ipow(dy, 2) + -3.97786 *lens_ipow(dx, 5)*lens_ipow(dy, 5) + 9.05996e-05 *x*lens_ipow(y, 2)*lens_ipow(dy, 7)+0.0f;
const float dx13 =  + 25.2549  + -0.000487611 *dx + -4.3886 *lens_ipow(dy, 2) + -1.45263 *lens_ipow(dx, 2) + 0.288511 *y*dy + 0.010328 *lens_ipow(y, 2) + 0.180011 *x*dx + 0.00654957 *lens_ipow(x, 2) + 0.797561 *lens_ipow(lambda, 3) + -0.000196179 *x*lens_ipow(dy, 2) + -0.0529869 *x*dx*lens_ipow(dy, 2) + 4.66035e-06 *lens_ipow(y, 4)*lens_ipow(dy, 2) + 0.00268865 *lens_ipow(x, 2)*lens_ipow(dy, 4) + -5.72727e-05 *lens_ipow(x, 3)*lens_ipow(dx, 3) + 0.443449 *lens_ipow(dy, 4)*lens_ipow(lambda, 3) + 0.439588 *lens_ipow(dx, 2)*lens_ipow(dy, 4)*lambda + 0.150615 *lens_ipow(dy, 8) + 0.183099 *lens_ipow(dx, 8) + -1.50573e-12 *lens_ipow(y, 8) + -0.0117741 *x*lens_ipow(dx, 3)*lens_ipow(lambda, 4) + 8.7303e-09 *lens_ipow(x, 3)*lens_ipow(y, 3)*dx*dy + -3.22089e-12 *lens_ipow(x, 8) + -4.37325 *lens_ipow(lambda, 10) + -3.31488 *lens_ipow(dx, 6)*lens_ipow(dy, 4) + 0.000634197 *x*lens_ipow(y, 2)*dx*lens_ipow(dy, 6)+0.0f;
const float dx14 =  + 2.39268 *dy*lens_ipow(lambda, 2) + 0.145081 *y*lens_ipow(lambda, 3) + 0.26607 *lens_ipow(dy, 5)*lens_ipow(lambda, 2) + 0.0879175 *lens_ipow(dx, 2)*lens_ipow(dy, 5) + -0.0470966 *x*lens_ipow(dx, 3)*dy*lens_ipow(lambda, 3) + 3.2653e-07 *lens_ipow(x, 2)*lens_ipow(y, 3)*lens_ipow(lambda, 3) + -43.7325 *dy*lens_ipow(lambda, 9) + -1.83804 *y*lens_ipow(lambda, 9) + -2.55134e-14 *lens_ipow(x, 6)*lens_ipow(y, 3)*lambda+0.0f;
const float dx20 =  + -0.0151942  + 1.29215e-06 *dy + -0.0161896 *lens_ipow(dy, 2) + -0.0435965 *lens_ipow(dx, 2) + -0.000377757 *y*dy + -8.28865e-06 *lens_ipow(y, 2) + -0.00122188 *x*dx + -2.13871e-05 *lens_ipow(x, 2) + 0.00244608 *lens_ipow(lambda, 3) + 0.000598336 *y*lens_ipow(dx, 2)*dy + -2.42909e-05 *lens_ipow(x, 2)*lens_ipow(dx, 2) + 7.37244e-07 *lens_ipow(x, 2)*y*dy + -0.00058324 *lens_ipow(dx, 2)*lens_ipow(lambda, 4) + 3.12934e-08 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(lambda, 2) + -3.21535e-06 *lens_ipow(x, 3)*dx*lens_ipow(dy, 2) + 0.00765051 *lens_ipow(dy, 6)*lambda + 7.13623e-09 *lens_ipow(y, 4)*lens_ipow(lambda, 4) + 1.07101e-09 *lens_ipow(y, 5)*lens_ipow(dy, 3) + -0.000270254 *x*dx*lens_ipow(dy, 2)*lens_ipow(lambda, 4) + -0.000793589 *x*dx*lens_ipow(dy, 6) + -1.99501e-10 *lens_ipow(x, 2)*lens_ipow(y, 4)*lens_ipow(dx, 2) + 0.000866781 *y*lens_ipow(dy, 3)*lens_ipow(lambda, 6) + 3.12448e-11 *lens_ipow(x, 6)*y*lens_ipow(dy, 3) + 3.10405e-14 *lens_ipow(x, 8)*lens_ipow(lambda, 2)+0.0f;
const float dx21 =  + 6.74336e-07  + -0.0298412 *dx*dy + -0.000626643 *y*dx + -0.000377757 *x*dy + -1.65773e-05 *x*y + -0.00024782 *y*lens_ipow(dx, 3) + -1.77481e-05 *lens_ipow(y, 2)*dx*dy + 0.000598336 *x*lens_ipow(dx, 2)*dy + 2.45748e-07 *lens_ipow(x, 3)*dy + 2.08622e-08 *lens_ipow(x, 3)*y*lens_ipow(lambda, 2) + 2.85449e-08 *x*lens_ipow(y, 3)*lens_ipow(lambda, 4) + 5.35503e-09 *x*lens_ipow(y, 4)*lens_ipow(dy, 3) + -2.66002e-10 *lens_ipow(x, 3)*lens_ipow(y, 3)*lens_ipow(dx, 2) + -3.14482e-06 *lens_ipow(y, 3)*lens_ipow(dx, 7) + 7.34132e-13 *lens_ipow(y, 8)*dx*dy + 0.000866781 *x*lens_ipow(dy, 3)*lens_ipow(lambda, 6) + 4.46354e-12 *lens_ipow(x, 7)*lens_ipow(dy, 3)+0.0f;
const float dx22 =  + 0.697653  + -0.739846 *lens_ipow(dy, 2) + -2.08595 *lens_ipow(dx, 2) + -0.0298412 *y*dy + -0.000313322 *lens_ipow(y, 2) + -0.0871929 *x*dx + -0.000610939 *lens_ipow(x, 2) + 0.0636618 *lens_ipow(lambda, 4) + -0.412956 *lens_ipow(dx, 2)*lens_ipow(dy, 2) + -0.304067 *lens_ipow(dx, 4) + -0.000371729 *lens_ipow(y, 2)*lens_ipow(dx, 2) + -5.91603e-06 *lens_ipow(y, 3)*dy + 0.00119667 *x*y*dx*dy + -1.6194e-05 *lens_ipow(x, 3)*dx + 0.0293939 *lens_ipow(dy, 6) + -0.00116648 *x*dx*lens_ipow(lambda, 4) + -8.03838e-07 *lens_ipow(x, 4)*lens_ipow(dy, 2) + -0.000135127 *lens_ipow(x, 2)*lens_ipow(dy, 2)*lens_ipow(lambda, 4) + -0.000396794 *lens_ipow(x, 2)*lens_ipow(dy, 6) + -1.33001e-10 *lens_ipow(x, 3)*lens_ipow(y, 4)*dx + 1.0121 *lens_ipow(dx, 8)*lens_ipow(dy, 2) + -5.50343e-06 *lens_ipow(y, 4)*lens_ipow(dx, 6) + 8.15702e-14 *lens_ipow(y, 9)*dy+0.0f;
const float dx23 =  + 1.29215e-06 *x + -1.47969 *dx*dy + -0.0298412 *y*dx + -0.0323792 *x*dy + -0.000377757 *x*y + -0.275304 *lens_ipow(dx, 3)*dy + -5.91603e-06 *lens_ipow(y, 3)*dx + 0.000598336 *x*y*lens_ipow(dx, 2) + 2.45748e-07 *lens_ipow(x, 3)*y + 0.176363 *dx*lens_ipow(dy, 5) + -1.60768e-06 *lens_ipow(x, 4)*dx*dy + 0.045903 *x*lens_ipow(dy, 5)*lambda + 3.21302e-09 *x*lens_ipow(y, 5)*lens_ipow(dy, 2) + -0.000270254 *lens_ipow(x, 2)*dx*dy*lens_ipow(lambda, 4) + -0.00238077 *lens_ipow(x, 2)*dx*lens_ipow(dy, 5) + 0.224911 *lens_ipow(dx, 9)*dy + 8.15702e-14 *lens_ipow(y, 9)*dx + 0.00260034 *x*y*lens_ipow(dy, 2)*lens_ipow(lambda, 6) + 1.33906e-11 *lens_ipow(x, 7)*y*lens_ipow(dy, 2)+0.0f;
const float dx24 =  + 0.00733823 *x*lens_ipow(lambda, 2) + 0.254647 *dx*lens_ipow(lambda, 3) + -0.00233296 *x*lens_ipow(dx, 2)*lens_ipow(lambda, 3) + 2.08622e-08 *lens_ipow(x, 3)*lens_ipow(y, 2)*lambda + 0.00765051 *x*lens_ipow(dy, 6) + 2.85449e-08 *x*lens_ipow(y, 4)*lens_ipow(lambda, 3) + -0.000540508 *lens_ipow(x, 2)*dx*lens_ipow(dy, 2)*lens_ipow(lambda, 3) + 0.00520068 *x*y*lens_ipow(dy, 3)*lens_ipow(lambda, 5) + 6.89789e-15 *lens_ipow(x, 9)*lambda+0.0f;
const float dx30 =  + -0.0293766 *dx*dy + -0.000344246 *y*dx + -0.000545988 *x*dy + -1.36561e-05 *x*y + -0.000140643 *y*dx*lens_ipow(lambda, 2) + 0.00079135 *y*dx*lens_ipow(dy, 2) + -0.000361019 *x*lens_ipow(dy, 3) + -1.89706e-05 *x*y*lens_ipow(dx, 2) + -7.55128e-09 *x*lens_ipow(y, 3) + 5.21165e-09 *lens_ipow(x, 4)*dx*dy + -0.00369451 *lens_ipow(dx, 3)*dy*lens_ipow(lambda, 4) + -0.00848157 *lens_ipow(dx, 5)*lens_ipow(dy, 3) + 4.52973e-12 *lens_ipow(x, 3)*lens_ipow(y, 4)*dy*lambda + -0.000306188 *x*lens_ipow(dy, 9) + 6.3222e-12 *lens_ipow(x, 2)*lens_ipow(y, 5)*lens_ipow(dx, 3)+0.0f;
const float dx31 =  + -0.0151386  + -0.0440123 *lens_ipow(dy, 2) + -0.0160967 *lens_ipow(dx, 2) + -0.00126989 *y*dy + -2.28749e-05 *lens_ipow(y, 2) + -0.000344246 *x*dx + -6.82806e-06 *lens_ipow(x, 2) + 0.00201202 *lens_ipow(lambda, 3) + -2.54835e-05 *lens_ipow(y, 2)*lens_ipow(dy, 2) + -0.000140643 *x*dx*lens_ipow(lambda, 2) + 0.00079135 *x*dx*lens_ipow(dy, 2) + -9.48528e-06 *lens_ipow(x, 2)*lens_ipow(dx, 2) + -1.13269e-08 *lens_ipow(x, 2)*lens_ipow(y, 2) + 0.00552454 *lens_ipow(dx, 6) + -0.00077547 *y*lens_ipow(dx, 2)*lens_ipow(dy, 3) + -0.00104801 *y*lens_ipow(dx, 4)*dy + 4.00354e-08 *lens_ipow(y, 4)*lens_ipow(lambda, 3) + 4.52973e-12 *lens_ipow(x, 4)*lens_ipow(y, 3)*dy*lambda + -0.00975626 *lens_ipow(dx, 2)*lens_ipow(lambda, 8) + -2.95227e-05 *lens_ipow(y, 2)*lens_ipow(dx, 8) + 1.53849e-11 *lens_ipow(y, 7)*lens_ipow(dx, 2)*dy + -1.09179e-15 *lens_ipow(y, 9)*dy + 1.0537e-11 *lens_ipow(x, 3)*lens_ipow(y, 4)*lens_ipow(dx, 3)+0.0f;
const float dx32 =  + -1.45311 *dx*dy + -0.0321933 *y*dx + -0.0293766 *x*dy + -0.000344246 *x*y + -0.314084 *dx*lens_ipow(dy, 3) + -0.000140643 *x*y*lens_ipow(lambda, 2) + 0.00079135 *x*y*lens_ipow(dy, 2) + -1.89706e-05 *lens_ipow(x, 2)*y*dx + 0.0331472 *y*lens_ipow(dx, 5) + -0.00077547 *lens_ipow(y, 2)*dx*lens_ipow(dy, 3) + -0.00209602 *lens_ipow(y, 2)*lens_ipow(dx, 3)*dy + 1.04233e-09 *lens_ipow(x, 5)*dy + 0.279504 *lens_ipow(dx, 7)*dy + -0.0110835 *x*lens_ipow(dx, 2)*dy*lens_ipow(lambda, 4) + -0.0424079 *x*lens_ipow(dx, 4)*lens_ipow(dy, 3) + 0.364039 *dx*lens_ipow(dy, 9) + -0.0195125 *y*dx*lens_ipow(lambda, 8) + -7.87272e-05 *lens_ipow(y, 3)*lens_ipow(dx, 7) + 3.84622e-12 *lens_ipow(y, 8)*dx*dy + 6.3222e-12 *lens_ipow(x, 3)*lens_ipow(y, 5)*lens_ipow(dx, 2)+0.0f;
const float dx33 =  + 0.697549  + -2.09239 *lens_ipow(dy, 2) + -0.726557 *lens_ipow(dx, 2) + -0.0880246 *y*dy + -0.000634947 *lens_ipow(y, 2) + -0.0293766 *x*dx + -0.000272994 *lens_ipow(x, 2) + 0.0560699 *lens_ipow(lambda, 4) + -0.305548 *lens_ipow(dy, 4) + -0.471126 *lens_ipow(dx, 2)*lens_ipow(dy, 2) + -1.6989e-05 *lens_ipow(y, 3)*dy + 0.0015827 *x*y*dx*dy + -0.000541528 *lens_ipow(x, 2)*lens_ipow(dy, 2) + -0.0011632 *lens_ipow(y, 2)*lens_ipow(dx, 2)*lens_ipow(dy, 2) + -0.000524004 *lens_ipow(y, 2)*lens_ipow(dx, 4) + 1.04233e-09 *lens_ipow(x, 5)*dx + 0.112307 *lens_ipow(dy, 4)*lens_ipow(lambda, 3) + 0.034938 *lens_ipow(dx, 8) + -0.00369451 *x*lens_ipow(dx, 3)*lens_ipow(lambda, 4) + -0.0254447 *x*lens_ipow(dx, 5)*lens_ipow(dy, 2) + 1.13243e-12 *lens_ipow(x, 4)*lens_ipow(y, 4)*lambda + 1.63818 *lens_ipow(dx, 2)*lens_ipow(dy, 8) + 1.92311e-12 *lens_ipow(y, 8)*lens_ipow(dx, 2) + -1.09179e-16 *lens_ipow(y, 10) + -0.00137785 *lens_ipow(x, 2)*lens_ipow(dy, 8)+0.0f;
const float dx34 =  + 0.00603605 *y*lens_ipow(lambda, 2) + 0.224279 *dy*lens_ipow(lambda, 3) + -0.000281285 *x*y*dx*lambda + 0.0673845 *lens_ipow(dy, 5)*lens_ipow(lambda, 2) + 2.40213e-08 *lens_ipow(y, 5)*lens_ipow(lambda, 2) + -0.014778 *x*lens_ipow(dx, 3)*dy*lens_ipow(lambda, 3) + 1.13243e-12 *lens_ipow(x, 4)*lens_ipow(y, 4)*dy + -0.0780501 *y*lens_ipow(dx, 2)*lens_ipow(lambda, 7)+0.0f;
const float dx40 =  + -2.59491e-06 *dy + -0.000618751 *dx + -1.64866e-05 *x + 3.45364e-05 *x*lens_ipow(dy, 2) + -8.36572e-08 *x*lens_ipow(y, 2) + -0.00156874 *lens_ipow(dx, 3)*lambda + -5.70196e-06 *y*dx*dy*lens_ipow(lambda, 2) + 9.51957e-12 *lens_ipow(x, 5) + 0.000217906 *dx*lens_ipow(lambda, 5) + 0.000341378 *y*dx*lens_ipow(dy, 5) + 8.93567e-08 *lens_ipow(x, 2)*lens_ipow(y, 2)*dx*lens_ipow(dy, 2) + 0.000149333 *y*lens_ipow(dx, 7)*dy + -0.000161091 *x*lens_ipow(dy, 8) + -2.18186e-12 *x*lens_ipow(y, 6)*lens_ipow(dy, 2) + 1.90013e-06 *lens_ipow(x, 3)*lens_ipow(dx, 2)*lens_ipow(dy, 4) + 1.50798e-06 *lens_ipow(x, 3)*lens_ipow(dx, 6) + 2.06605e-10 *lens_ipow(x, 4)*lens_ipow(y, 2)*lens_ipow(dx, 3) + -3.62016e-13 *lens_ipow(x, 4)*lens_ipow(y, 4)*dx + -9.26215e-12 *lens_ipow(x, 7)*lens_ipow(dx, 2)+0.0f;
const float dx41 =  + 6.51812e-07  + -0.000596146 *dy + 3.72738e-06 *dx + -1.14904e-05 *y + -0.00139188 *lens_ipow(dx, 2)*dy + -8.36572e-08 *lens_ipow(x, 2)*y + 2.49768e-05 *y*lens_ipow(dy, 2)*lambda + -2.01756e-07 *lens_ipow(y, 3)*lens_ipow(dx, 2) + -3.20796e-10 *lens_ipow(y, 5) + -5.70196e-06 *x*dx*dy*lens_ipow(lambda, 2) + -0.00275713 *lens_ipow(dx, 4)*lens_ipow(dy, 3) + 0.000341378 *x*dx*lens_ipow(dy, 5) + 5.95711e-08 *lens_ipow(x, 3)*y*dx*lens_ipow(dy, 2) + 6.75466e-07 *lens_ipow(y, 3)*lens_ipow(dy, 6) + 0.000149333 *x*lens_ipow(dx, 7)*dy + -6.54559e-12 *lens_ipow(x, 2)*lens_ipow(y, 5)*lens_ipow(dy, 2) + 8.2642e-11 *lens_ipow(x, 5)*y*lens_ipow(dx, 3) + -2.89613e-13 *lens_ipow(x, 5)*lens_ipow(y, 3)*dx+0.0f;
const float dx42 =  + 1.65404e-05  + -0.0400833 *dx + 3.72738e-06 *y + -0.000618751 *x + -0.105891 *dx*lens_ipow(dy, 2) + -0.00278377 *y*dx*dy + -0.00470622 *x*lens_ipow(dx, 2)*lambda + -0.0979634 *dx*lens_ipow(dy, 4) + -0.266522 *lens_ipow(dx, 5) + -1.00878e-07 *lens_ipow(y, 4)*dx + -5.70196e-06 *x*y*dy*lens_ipow(lambda, 2) + 0.000217906 *x*lens_ipow(lambda, 5) + -0.131045 *lens_ipow(dx, 3)*lens_ipow(lambda, 4) + -0.0110285 *y*lens_ipow(dx, 3)*lens_ipow(dy, 3) + 0.000341378 *x*y*lens_ipow(dy, 5) + 2.97856e-08 *lens_ipow(x, 3)*lens_ipow(y, 2)*lens_ipow(dy, 2) + 0.00104533 *x*y*lens_ipow(dx, 6)*dy + 9.50064e-07 *lens_ipow(x, 4)*dx*lens_ipow(dy, 4) + 2.26197e-06 *lens_ipow(x, 4)*lens_ipow(dx, 5) + 1.23963e-10 *lens_ipow(x, 5)*lens_ipow(y, 2)*lens_ipow(dx, 2) + -7.24033e-14 *lens_ipow(x, 5)*lens_ipow(y, 4) + -2.31554e-12 *lens_ipow(x, 8)*dx+0.0f;
const float dx43 =  + -0.0388975 *dy + -0.000596146 *y + -2.59491e-06 *x + -0.105891 *lens_ipow(dx, 2)*dy + -0.00139188 *y*lens_ipow(dx, 2) + 3.45364e-05 *lens_ipow(x, 2)*dy + 2.49768e-05 *lens_ipow(y, 2)*dy*lambda + -0.157702 *lens_ipow(dy, 5) + -0.195927 *lens_ipow(dx, 2)*lens_ipow(dy, 3) + -5.70196e-06 *x*y*dx*lens_ipow(lambda, 2) + -0.0082714 *y*lens_ipow(dx, 4)*lens_ipow(dy, 2) + 0.00170689 *x*y*dx*lens_ipow(dy, 4) + 5.95711e-08 *lens_ipow(x, 3)*lens_ipow(y, 2)*dx*dy + 1.0132e-06 *lens_ipow(y, 4)*lens_ipow(dy, 5) + 0.000149333 *x*y*lens_ipow(dx, 7) + -0.000644366 *lens_ipow(x, 2)*lens_ipow(dy, 7) + -2.18186e-12 *lens_ipow(x, 2)*lens_ipow(y, 6)*dy + 1.90013e-06 *lens_ipow(x, 4)*lens_ipow(dx, 2)*lens_ipow(dy, 3)+0.0f;
const float dx44 =  + 0.0825872 *lens_ipow(lambda, 2) + 1.24884e-05 *lens_ipow(y, 2)*lens_ipow(dy, 2) + -0.00156874 *x*lens_ipow(dx, 3) + -1.14039e-05 *x*y*dx*dy*lambda + 0.00108953 *x*dx*lens_ipow(lambda, 4) + -0.131045 *lens_ipow(dx, 4)*lens_ipow(lambda, 3) + -2.15489 *lens_ipow(lambda, 10)+0.0f;