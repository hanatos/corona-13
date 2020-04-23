const float out_x =  + -5.21966e-05  + -0.00269054 *dy + 64.1735 *dx + -0.77313 *x + -0.0565598 *dx*dy + -0.205153 *lens_ipow(dx, 2) + -0.000445248 *x*dy + 72.7431 *dx*lens_ipow(dy, 2) + 30.7833 *lens_ipow(dx, 3) + 0.503941 *y*dx*dy + 0.000693242 *lens_ipow(y, 2)*dx + -1.33371 *x*lens_ipow(dy, 2) + -1.82635 *x*lens_ipow(dx, 2) + 0.002417 *x*y*dy + -0.000192708 *x*lens_ipow(y, 2) + 0.0172448 *lens_ipow(x, 2)*dx + -0.000263213 *lens_ipow(x, 3) + 38.5907 *dx*lens_ipow(lambda, 3) + -0.145518 *x*lens_ipow(dx, 3) + 0.125654 *x*lens_ipow(lambda, 4) + -19.7548 *x*lens_ipow(dx, 2)*lens_ipow(dy, 2) + -0.343184 *x*y*lens_ipow(dx, 2)*dy*lambda + 3338.13 *dx*lens_ipow(dy, 4)*lens_ipow(lambda, 2) + 2.0031 *lens_ipow(dx, 2)*lens_ipow(lambda, 5) + 5344.09 *lens_ipow(dx, 3)*lens_ipow(dy, 2)*lens_ipow(lambda, 2) + 465.556 *y*lens_ipow(dx, 3)*lens_ipow(dy, 3) + -241.703 *x*lens_ipow(dx, 6) + -0.0340574 *x*lens_ipow(y, 2)*lens_ipow(dx, 4) + -2.31564e-06 *x*lens_ipow(y, 4)*lens_ipow(dy, 2) + 8.12252e-08 *lens_ipow(x, 2)*lens_ipow(y, 4)*dx + 1.34481e-07 *lens_ipow(x, 3)*lens_ipow(y, 3)*dy + -1.00577e-05 *lens_ipow(x, 4)*y*dx*dy + 1.30979e-07 *lens_ipow(x, 6)*dx + 2501.24 *lens_ipow(dx, 5)*lens_ipow(lambda, 3) + 3.56501 *y*dx*dy*lens_ipow(lambda, 5) + 4.97816e-08 *lens_ipow(y, 6)*dx*lens_ipow(lambda, 2) + -2.96866e-08 *lens_ipow(x, 7)*lens_ipow(dy, 2) + -207.461 *dx*lens_ipow(lambda, 10) + -0.613647 *x*lens_ipow(lambda, 10) + 8.78977e-08 *lens_ipow(x, 6)*lens_ipow(y, 2)*lens_ipow(dx, 3);
const float out_y =  + 80.5245 *dy + -0.223757 *y + 1.76304e-05 *x + -0.183412 *dx*dy + -0.000847909 *y*dx + 0.00065769 *x*dy + 123.666 *lens_ipow(dy, 3) + 65.1156 *lens_ipow(dx, 2)*dy + 0.671081 *y*lens_ipow(dy, 2) + -0.0323779 *y*lens_ipow(dx, 2) + 0.00819248 *lens_ipow(y, 2)*dy + -0.000147184 *lens_ipow(y, 3) + -0.0196473 *x*lens_ipow(dy, 2) + -1.29881 *x*dx*dy + -0.000258926 *x*y*dy + 0.00207194 *x*y*dx + 0.00412981 *lens_ipow(x, 2)*dy + -0.000252213 *lens_ipow(x, 2)*y + 47.4524 *dy*lens_ipow(lambda, 3) + -0.111497 *x*dx*lens_ipow(dy, 2) + 0.503787 *y*lens_ipow(lambda, 4) + 0.102108 *lens_ipow(x, 2)*lens_ipow(dy, 3) + 0.842141 *y*lens_ipow(dx, 2)*lens_ipow(lambda, 3) + 2423.9 *lens_ipow(dy, 5)*lens_ipow(lambda, 2) + 3576.51 *lens_ipow(dx, 2)*lens_ipow(dy, 3)*lens_ipow(lambda, 2) + 3.16618 *y*lens_ipow(dy, 2)*lens_ipow(lambda, 4) + 0.750642 *x*dx*dy*lens_ipow(lambda, 4) + -6.14108e-06 *x*lens_ipow(y, 4)*dx*dy + -6.25724e-06 *lens_ipow(x, 2)*lens_ipow(y, 3)*lens_ipow(dx, 2) + 9.86536e-08 *lens_ipow(x, 2)*lens_ipow(y, 4)*dy + -4.96344e-06 *lens_ipow(x, 5)*dx*dy + 6.78151e-08 *lens_ipow(x, 5)*y*dx + 2713.77 *lens_ipow(dx, 4)*dy*lens_ipow(lambda, 3) + -7.51022e-09 *lens_ipow(y, 7)*lens_ipow(dy, 2) + -28069.5 *x*lens_ipow(dx, 5)*lens_ipow(dy, 3) + -0.0107871 *lens_ipow(x, 2)*dy*lens_ipow(lambda, 6) + 5.4134e-10 *lens_ipow(x, 6)*lens_ipow(y, 2)*dy + -248.74 *dy*lens_ipow(lambda, 10) + -2.49029 *y*lens_ipow(lambda, 10) + 8.84606e-05 *lens_ipow(y, 4)*dy*lens_ipow(lambda, 6);
const float out_dx =  + 3.16697e-06  + -3.97472e-05 *dy + -0.250874 *dx + -0.0124777 *x + -0.000844746 *dx*dy + -0.00180622 *lens_ipow(dx, 2) + 6.51006e-06 *y*dy + 7.3116e-08 *x*y + 2.93183 *dx*lens_ipow(dy, 2) + 1.83648 *lens_ipow(dx, 3) + 0.0214165 *y*dx*dy + -1.01791e-05 *lens_ipow(y, 2)*dx + -0.00867881 *x*lens_ipow(dy, 2) + -0.0236576 *x*lens_ipow(dx, 2) + -0.000146358 *x*y*dy + 1.30571e-06 *x*y*dx + 9.05968e-07 *x*lens_ipow(y, 2) + -0.000124999 *lens_ipow(x, 2)*dx + 2.2047e-06 *lens_ipow(x, 3) + 0.508911 *dx*lens_ipow(lambda, 3) + -0.00154828 *x*lens_ipow(dx, 3) + 0.0757805 *y*lens_ipow(dx, 3)*dy + 0.00111915 *x*lens_ipow(lambda, 4) + -2.29544e-05 *lens_ipow(x, 2)*y*dx*dy + -2.56549e-05 *lens_ipow(x, 3)*lens_ipow(dx, 2) + 0.0316159 *x*lens_ipow(dx, 2)*lens_ipow(lambda, 3) + 40.8083 *lens_ipow(dx, 3)*lens_ipow(dy, 2)*lens_ipow(lambda, 2) + 17.1954 *lens_ipow(dx, 5)*lens_ipow(lambda, 2) + 4.19667e-06 *lens_ipow(y, 4)*dx*lens_ipow(dy, 2) + -2.69981e-05 *x*lens_ipow(y, 2)*lens_ipow(dy, 2)*lens_ipow(lambda, 2) + -0.0536859 *lens_ipow(x, 2)*lens_ipow(dx, 5) + 0.0396166 *y*dx*dy*lens_ipow(lambda, 5) + 2.47959e-07 *lens_ipow(y, 4)*lens_ipow(dx, 2)*lens_ipow(lambda, 2) + 0.435418 *x*lens_ipow(dx, 2)*lens_ipow(dy, 2)*lens_ipow(lambda, 3) + -0.0237313 *x*lens_ipow(y, 2)*lens_ipow(dx, 6) + -0.000260487 *lens_ipow(x, 2)*dx*lens_ipow(lambda, 6) + 0.850931 *x*lens_ipow(dy, 4)*lens_ipow(lambda, 5) + -2.70524 *dx*lens_ipow(lambda, 10) + -0.00626974 *x*lens_ipow(lambda, 10) + -0.000322547 *x*lens_ipow(y, 2)*lens_ipow(dx, 2)*lens_ipow(lambda, 6);
const float out_dy =  + -0.228998 *dy + -0.0116081 *y + 2.00034e-07 *x + -0.002396 *dx*dy + -1.09544e-05 *y*dx + 1.3016e-05 *x*dy + 2.22304e-07 *x*y + 3.75246 *lens_ipow(dy, 3) + 2.26751 *lens_ipow(dx, 2)*dy + 0.0410277 *y*lens_ipow(dy, 2) + 0.00816191 *y*lens_ipow(dx, 2) + 6.29182e-06 *lens_ipow(y, 2)*dy + 5.40986e-07 *lens_ipow(y, 3) + -0.000260123 *x*lens_ipow(dy, 2) + -0.0177934 *x*dx*dy + -2.73159e-06 *x*y*dy + -0.000124414 *x*y*dx + -0.000113656 *lens_ipow(x, 2)*dy + 1.32418e-06 *lens_ipow(x, 2)*y + 0.610365 *dy*lens_ipow(lambda, 3) + -0.0013763 *x*dx*lens_ipow(dy, 2) + -2.19979e-05 *x*y*lens_ipow(dx, 2) + 0.830329 *lens_ipow(dx, 2)*dy*lens_ipow(lambda, 2) + 0.00612239 *y*lens_ipow(lambda, 4) + 1.21493e-07 *x*lens_ipow(y, 3)*dx + -1.52484e-05 *lens_ipow(x, 2)*y*lens_ipow(dy, 2) + 0.0165547 *y*lens_ipow(dx, 2)*lens_ipow(lambda, 4) + -0.000202563 *lens_ipow(y, 3)*lens_ipow(dy, 4) + 0.0269228 *x*dx*dy*lens_ipow(lambda, 4) + 0.00431654 *lens_ipow(x, 2)*lens_ipow(dy, 3)*lens_ipow(lambda, 2) + 7.67607e-12 *lens_ipow(x, 2)*lens_ipow(y, 5) + -1.20998e-07 *lens_ipow(x, 5)*dx*dy + 2.6574e-12 *lens_ipow(x, 6)*y + 0.0689428 *y*lens_ipow(dy, 2)*lens_ipow(lambda, 5) + -3.37741e-06 *lens_ipow(x, 2)*y*lens_ipow(lambda, 5) + 7.45551e-15 *lens_ipow(y, 9) + -0.000445349 *lens_ipow(x, 2)*dy*lens_ipow(lambda, 6) + 1.90644e-06 *lens_ipow(y, 4)*dy*lens_ipow(lambda, 5) + -3.16057 *dy*lens_ipow(lambda, 10) + -0.0301143 *y*lens_ipow(lambda, 10);
const float out_transmittance =  + 0.0776915  + 0.225284 *lambda + 2.76322e-06 *dx + 2.71748e-07 *y + -0.0118818 *lens_ipow(dy, 2) + -0.000513865 *dx*dy + -0.00666743 *lens_ipow(dx, 2) + -3.9058e-06 *y*dy + -3.90868e-06 *y*dx + -4.43335e-07 *lens_ipow(y, 2) + -2.01219e-06 *x*dy + -1.54203e-05 *x*dx + -1.91457e-06 *lens_ipow(x, 2) + -0.182881 *lens_ipow(lambda, 3) + -0.00515078 *dx*lens_ipow(dy, 2) + -4.26785e-05 *x*dx*dy + 1.04773e-05 *x*lens_ipow(dx, 2) + 1.91714e-07 *lens_ipow(x, 2)*dy + -2.16923 *lens_ipow(dx, 2)*lens_ipow(dy, 2) + -0.0369953 *y*lens_ipow(dy, 3) + -0.0286968 *y*lens_ipow(dx, 2)*dy + -0.000319611 *lens_ipow(y, 2)*lens_ipow(dy, 2) + -7.99536e-05 *lens_ipow(y, 2)*lens_ipow(dx, 2) + -6.76217e-07 *lens_ipow(y, 3)*dy + -0.0021188 *x*dx*lens_ipow(dy, 2) + -3.78524e-05 *x*y*dx*dy + 4.72996e-07 *lens_ipow(x, 2)*y*dy + -1.08086e-06 *lens_ipow(x, 3)*dx + -8.54732e-07 *lens_ipow(x, 2)*y*dx*dy + -50.6803 *lens_ipow(dy, 6) + -18.7803 *lens_ipow(dx, 6) + -0.0020069 *lens_ipow(x, 2)*lens_ipow(dx, 4) + -1.29687e-11 *lens_ipow(x, 2)*lens_ipow(y, 4) + 2.73684e-09 *lens_ipow(x, 3)*lens_ipow(y, 2)*dx + -1.77131e-10 *lens_ipow(x, 6) + -1304.07 *lens_ipow(dx, 4)*lens_ipow(dy, 4) + -1.73455e-05 *lens_ipow(x, 4)*lens_ipow(dy, 4) + -163.444 *y*lens_ipow(dy, 9) + -7.04828e-14 *lens_ipow(x, 9)*dx + 0.318525 *lens_ipow(lambda, 11);
