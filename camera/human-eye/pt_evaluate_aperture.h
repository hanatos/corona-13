const float out_x =  + 3.28602e-05  + 20.6514 *dx + -5.63678e-06 *y + 0.958992 *x + -0.000264769 *dx*dy + -6.10268e-07 *lens_ipow(y, 2) + 1.78183e-07 *lens_ipow(x, 2) + -2.00887 *dx*lens_ipow(dy, 2) + -0.805632 *lens_ipow(dx, 3) + 1.37149e-05 *y*lens_ipow(dx, 2) + 0.00388807 *lens_ipow(y, 2)*dx + 0.000996169 *x*lens_ipow(lambda, 2) + -0.0796915 *x*lens_ipow(dy, 2) + 0.000146804 *x*lens_ipow(y, 2) + 0.00331119 *lens_ipow(x, 2)*dx + 7.79241e-05 *lens_ipow(x, 3) + 0.00187821 *lens_ipow(dx, 3)*dy*lambda + 0.0496628 *lens_ipow(dx, 5) + 0.0415684 *y*lens_ipow(dx, 3)*dy + -0.00471672 *x*lens_ipow(dx, 2)*lens_ipow(lambda, 2) + 1.51626e-05 *lens_ipow(x, 2)*lens_ipow(y, 2)*dx + -0.000368344 *lens_ipow(x, 3)*lens_ipow(dy, 2) + -2.83985e-05 *lens_ipow(x, 3)*y*dy + 7.98496e-06 *lens_ipow(y, 2)*dx*lens_ipow(lambda, 3) + 8.78111e-06 *lens_ipow(y, 2)*lens_ipow(dx, 4) + -0.006344 *x*lens_ipow(dy, 6) + 1.79008e-08 *x*lens_ipow(y, 4)*lens_ipow(dx, 2) + -2.34358e-09 *x*lens_ipow(y, 5)*dy + 9.60394e-05 *lens_ipow(y, 2)*lens_ipow(dx, 7) + -3.11898e-05 *lens_ipow(y, 3)*lens_ipow(dx, 3)*lens_ipow(dy, 3) + -1.22963e-07 *lens_ipow(x, 4)*lens_ipow(dy, 2)*lens_ipow(lambda, 3) + 3.43725e-05 *lens_ipow(x, 3)*lens_ipow(lambda, 7) + 0.264863 *dx*lens_ipow(lambda, 10) + -0.0526424 *lens_ipow(dx, 7)*lens_ipow(lambda, 4) + 0.00836218 *lens_ipow(dx, 11) + 4.95146e-15 *lens_ipow(y, 10)*dx + -0.00054028 *lens_ipow(x, 2)*dx*lens_ipow(dy, 8) + 1.0679e-12 *lens_ipow(x, 8)*lens_ipow(dx, 3) + 4.20136e-17 *lens_ipow(x, 9)*lens_ipow(y, 2) + 3.44444e-16 *lens_ipow(x, 11);
const float out_y =  + -1.37254e-05  + 20.6714 *dy + 0.00014689 *dx + 0.960287 *y + -0.000605106 *lens_ipow(dy, 2) + 0.000288029 *dx*dy + -0.256365 *lens_ipow(dy, 3) + -2.05504 *lens_ipow(dx, 2)*dy + -0.0920141 *y*lens_ipow(dx, 2) + -1.12147e-05 *lens_ipow(y, 3) + 2.80653e-05 *x*lens_ipow(dy, 2) + 0.00295942 *lens_ipow(x, 2)*dy + 0.000131986 *lens_ipow(x, 2)*y + -0.022321 *dy*lens_ipow(lambda, 4) + -0.624661 *lens_ipow(dy, 5) + -0.000149128 *lens_ipow(y, 3)*lens_ipow(dy, 2) + 6.90913e-05 *x*lens_ipow(dy, 3)*lambda + -4.04562e-06 *x*lens_ipow(y, 3)*dx + 0.0012934 *lens_ipow(x, 2)*dy*lens_ipow(lambda, 2) + 0.00104629 *lens_ipow(x, 2)*lens_ipow(dy, 3) + 2.47413e-06 *lens_ipow(x, 4)*dy + -7.51521e-06 *lens_ipow(x, 3)*y*dx*lambda + -7.54881e-05 *lens_ipow(y, 3)*lens_ipow(dx, 4) + 5.4747e-10 *lens_ipow(y, 7) + -9.85075e-05 *lens_ipow(y, 3)*lens_ipow(dx, 2)*lens_ipow(lambda, 3) + -0.00162525 *x*y*dx*lens_ipow(lambda, 5) + 2.03814e-13 *lens_ipow(x, 2)*lens_ipow(y, 6) + -0.00377389 *y*lens_ipow(dx, 8) + 8.00912e-05 *lens_ipow(x, 3)*lens_ipow(dx, 3)*lens_ipow(dy, 3) + 0.000879528 *lens_ipow(y, 2)*dy*lens_ipow(lambda, 7) + 2.3683e-10 *lens_ipow(y, 7)*lens_ipow(lambda, 3) + -0.0266561 *y*lens_ipow(dx, 6)*lens_ipow(lambda, 4) + 0.000407786 *lens_ipow(y, 2)*lens_ipow(dy, 7)*lens_ipow(lambda, 2) + -6.63749e-13 *lens_ipow(y, 9)*lens_ipow(dy, 2) + 0.0882062 *x*dx*lens_ipow(dy, 3)*lens_ipow(lambda, 6) + -6.73945e-05 *lens_ipow(x, 2)*lens_ipow(dy, 9) + 8.06044e-18 *lens_ipow(x, 2)*lens_ipow(y, 9) + -1.21235e-14 *lens_ipow(x, 5)*lens_ipow(y, 4)*lens_ipow(lambda, 2) + 2.40636e-09 *lens_ipow(x, 6)*lens_ipow(dy, 5) + -8.72095e-10 *lens_ipow(x, 6)*y*lens_ipow(lambda, 4);
const float out_dx =  + -1.3018e-05  + 0.000143682 *dy + 0.652959 *dx + -0.0178346 *x + 0.000227036 *lens_ipow(dy, 2) + -0.000146325 *dx*dy + -0.000127855 *lens_ipow(dx, 2) + 0.0158749 *dx*lens_ipow(lambda, 2) + -0.0783214 *dx*lens_ipow(dy, 2) + -0.00119969 *lens_ipow(dx, 2)*dy + -0.113371 *lens_ipow(dx, 3) + -0.000325238 *lens_ipow(y, 2)*dx + 1.63558e-08 *lens_ipow(y, 3) + -2.94335e-05 *x*lens_ipow(y, 2) + -1.2444e-07 *lens_ipow(x, 2)*y + -1.78927e-05 *lens_ipow(x, 3) + 2.19682e-05 *y*lens_ipow(lambda, 3) + -2.552e-09 *lens_ipow(y, 4) + -0.0182644 *lens_ipow(dx, 5) + -1.92814e-05 *lens_ipow(y, 2)*lens_ipow(dx, 3) + 0.00194052 *x*lens_ipow(lambda, 4) + 1.69424e-06 *lens_ipow(y, 2)*lens_ipow(dx, 4) + 5.2939e-12 *lens_ipow(x, 6) + -6.40093e-05 *lens_ipow(x, 2)*dx*lens_ipow(dy, 4) + 4.91689e-06 *lens_ipow(x, 3)*lens_ipow(dx, 2)*lens_ipow(lambda, 2) + -4.47452e-11 *lens_ipow(x, 5)*lens_ipow(y, 2) + -3.2386e-11 *lens_ipow(x, 7) + 0.00345605 *lens_ipow(dx, 3)*lens_ipow(dy, 3)*lens_ipow(lambda, 2) + 3.16853e-13 *lens_ipow(y, 8)*dx + 0.000146273 *x*lens_ipow(dy, 8) + 5.84563e-11 *lens_ipow(x, 7)*lens_ipow(dy, 2) + 0.0754515 *dx*lens_ipow(lambda, 9) + 0.0740975 *lens_ipow(dx, 5)*lens_ipow(lambda, 5) + -9.91617e-09 *lens_ipow(x, 4)*lens_ipow(dy, 6) + 0.112051 *lens_ipow(dx, 9)*lens_ipow(dy, 2) + 0.0341507 *lens_ipow(dx, 11) + -4.71868e-05 *y*lens_ipow(dx, 10) + -0.000252111 *lens_ipow(y, 2)*lens_ipow(dx, 9) + 5.09132e-11 *lens_ipow(x, 7)*lens_ipow(dx, 4) + -1.76295e-16 *lens_ipow(x, 7)*lens_ipow(y, 4);
const float out_dy =  + -4.41351e-05  + 0.651747 *dy + -0.0180723 *y + 1.58486e-07 *x + -9.54883e-06 *lens_ipow(dy, 2) + 3.82465e-05 *lens_ipow(dx, 2) + -2.25784e-07 *x*y + -0.112733 *lens_ipow(dy, 3) + -0.0737849 *lens_ipow(dx, 2)*dy + 0.00151504 *y*lens_ipow(lambda, 2) + -1.81833e-05 *lens_ipow(y, 3) + 1.50527e-05 *x*lens_ipow(dy, 2) + -0.000325749 *lens_ipow(x, 2)*dy + -3.01079e-05 *lens_ipow(x, 2)*y + -1.25466e-08 *lens_ipow(x, 3) + 0.000233115 *lens_ipow(lambda, 4) + 0.0396612 *dy*lens_ipow(lambda, 3) + -3.59882e-05 *lens_ipow(x, 2)*dy*lambda + -0.0163436 *lens_ipow(dy, 5) + 0.000796254 *y*lens_ipow(dx, 4) + -6.42563e-05 *lens_ipow(x, 2)*lens_ipow(dy, 3) + 7.83266e-08 *lens_ipow(x, 4)*dy + -0.000858782 *y*lens_ipow(dx, 2)*lens_ipow(lambda, 4) + 1.17285e-06 *lens_ipow(y, 3)*lens_ipow(lambda, 4) + -6.9994e-13 *x*lens_ipow(y, 6) + -4.46182e-11 *lens_ipow(x, 2)*lens_ipow(y, 5) + 0.0294655 *lens_ipow(dx, 4)*lens_ipow(dy, 5) + 4.5837e-11 *lens_ipow(y, 7)*lens_ipow(dx, 2) + -9.50418e-14 *lens_ipow(y, 9) + 4.98666e-06 *lens_ipow(x, 2)*lens_ipow(dy, 3)*lens_ipow(lambda, 4) + 1.11132e-14 *lens_ipow(x, 6)*lens_ipow(y, 3) + 0.00159688 *y*lens_ipow(lambda, 9) + 6.58317e-11 *lens_ipow(y, 6)*lens_ipow(dx, 4) + -2.45029e-08 *lens_ipow(x, 4)*lens_ipow(dy, 6) + 1.12633e-10 *lens_ipow(x, 6)*dy*lens_ipow(lambda, 3) + 0.0171727 *lens_ipow(dy, 3)*lens_ipow(lambda, 8) + 0.0379702 *lens_ipow(dy, 11) + 3.64315e-06 *lens_ipow(y, 3)*lens_ipow(dx, 6)*lens_ipow(lambda, 2) + 6.70795e-11 *lens_ipow(y, 7)*lens_ipow(dy, 4) + 2.14337e-09 *lens_ipow(x, 5)*lens_ipow(dy, 4)*lens_ipow(lambda, 2);
const float out_transmittance =  + 0.99806  + -5.59974e-06 *dy + -2.59033e-07 *x + 5.36132e-05 *lens_ipow(dy, 2) + 0.000842483 *lens_ipow(dx, 2) + -4.43169e-08 *lens_ipow(y, 2) + -4.85853e-08 *x*y + -4.21071e-06 *lens_ipow(x, 2) + 0.000601819 *lens_ipow(lambda, 3) + -1.02365e-06 *y*lens_ipow(dx, 2) + -6.1525e-08 *lens_ipow(y, 2)*dx + -0.0451489 *lens_ipow(dx, 2)*lens_ipow(dy, 2) + -0.0140775 *lens_ipow(dx, 4) + 0.000234318 *x*y*dx*dy + 4.77574e-05 *lens_ipow(x, 2)*lens_ipow(dx, 2) + -3.1512e-07 *lens_ipow(x, 2)*lens_ipow(y, 2) + -4.08491e-09 *lens_ipow(x, 3)*dy + 5.6236e-05 *lens_ipow(dx, 2)*lens_ipow(dy, 3) + -4.15987e-08 *lens_ipow(x, 2)*lens_ipow(dx, 3) + -0.0173334 *lens_ipow(dy, 6) + -4.47952e-06 *lens_ipow(y, 3)*lens_ipow(dy, 3) + -2.89511e-10 *lens_ipow(y, 6) + -5.23784e-06 *x*lens_ipow(dy, 5) + -0.000481397 *x*lens_ipow(dx, 5) + 3.62145e-06 *lens_ipow(x, 2)*lens_ipow(dy, 4) + -2.77125e-10 *lens_ipow(x, 6) + 4.98755e-10 *x*lens_ipow(y, 3)*lens_ipow(dy, 3) + -1.82512e-09 *lens_ipow(x, 3)*y*dx*dy*lambda + 2.08385e-09 *lens_ipow(y, 3)*lens_ipow(dx, 5) + 1.02364e-11 *lens_ipow(y, 6)*lens_ipow(dx, 2) + -6.40171e-12 *x*lens_ipow(y, 5)*lens_ipow(lambda, 3) + 1.50043e-08 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(lambda, 5) + -0.00205336 *lens_ipow(dx, 4)*lens_ipow(dy, 6) + -0.00860845 *lens_ipow(dx, 10) + 4.22825e-06 *lens_ipow(y, 2)*lens_ipow(dx, 8) + -1.10233e-11 *lens_ipow(x, 7)*lens_ipow(dx, 3) + 2.02955e-14 *lens_ipow(x, 8)*lens_ipow(dy, 2) + -0.00512579 *lens_ipow(lambda, 11) + -0.00230046 *lens_ipow(dx, 6)*lens_ipow(dy, 2)*lens_ipow(lambda, 3) + 2.06484e-10 *lens_ipow(y, 5)*lens_ipow(dx, 4)*lens_ipow(lambda, 2);
