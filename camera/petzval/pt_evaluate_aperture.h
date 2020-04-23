const float out_x =  + -2.39596e-05  + 0.000170721 *dy + 50.3715 *dx + 1.20573e-06 *y + 0.792054 *x + 0.003814 *lens_ipow(dx, 2) + -1.02565e-06 *lens_ipow(x, 2) + 0.101002 *dx*lens_ipow(dy, 2) + 0.240926 *lens_ipow(dx, 3) + 0.534569 *y*dx*dy + 0.00916713 *lens_ipow(y, 2)*dx + 0.123821 *x*lens_ipow(dy, 2) + 0.658175 *x*lens_ipow(dx, 2) + 0.013701 *x*y*dy + 0.000166393 *x*lens_ipow(y, 2) + -7.99133e-06 *lens_ipow(x, 2)*dy + 0.022954 *lens_ipow(x, 2)*dx + 0.000168299 *lens_ipow(x, 3) + -0.0537504 *dx*lens_ipow(dy, 3) + 0.0538104 *lens_ipow(dx, 3)*dy + 0.0553584 *x*lens_ipow(lambda, 3) + 3.57819 *dx*lens_ipow(lambda, 4) + 0.00194404 *lens_ipow(y, 2)*dx*lens_ipow(dy, 2) + 6.06162e-06 *x*lens_ipow(y, 2)*lens_ipow(lambda, 3) + 0.00031874 *lens_ipow(x, 2)*lens_ipow(dy, 4) + 0.00750986 *x*y*lens_ipow(dy, 5) + 0.0268014 *lens_ipow(x, 2)*dx*lens_ipow(dy, 4) + 4.90502e-05 *lens_ipow(x, 3)*y*lens_ipow(dx, 2)*dy + 1.80118e-09 *lens_ipow(x, 6)*dx + 0.0583032 *lens_ipow(y, 2)*lens_ipow(dx, 6)*dy + -0.144778 *lens_ipow(x, 2)*lens_ipow(dx, 5)*lens_ipow(dy, 2) + 1.88903e-10 *lens_ipow(x, 5)*lens_ipow(y, 2)*lens_ipow(lambda, 3) + 9.66942e-11 *lens_ipow(x, 7)*lens_ipow(lambda, 3) + -17.3519 *dx*lens_ipow(lambda, 10) + -2.25406 *lens_ipow(y, 2)*lens_ipow(dx, 9) + -0.298715 *x*lens_ipow(lambda, 10) + 0.691107 *lens_ipow(x, 2)*lens_ipow(dx, 9) + 1.13005e-15 *lens_ipow(x, 3)*lens_ipow(y, 8) + 2.36409e-08 *lens_ipow(x, 7)*lens_ipow(dy, 4) + 1.25173e-15 *lens_ipow(x, 9)*lens_ipow(y, 2);
const float out_y =  + 5.36935e-05  + 50.4606 *dy + -7.98426e-05 *dx + 0.793472 *y + -7.02206e-07 *x + 2.57101e-05 *y*dy + -6.73559e-06 *x*dy + -11.1018 *lens_ipow(dy, 3) + -1.28628 *lens_ipow(dx, 2)*dy + 0.104314 *y*lens_ipow(dx, 2) + 0.0185838 *lens_ipow(y, 2)*dy + 0.00014928 *lens_ipow(y, 3) + 0.504188 *x*dx*dy + 0.013111 *x*y*dx + 0.00900447 *lens_ipow(x, 2)*dy + 0.000161239 *lens_ipow(x, 2)*y + 0.0490546 *y*lens_ipow(lambda, 3) + 1.01948 *y*lens_ipow(dy, 2)*lambda + 2.62835 *dy*lens_ipow(lambda, 4) + 35.2614 *lens_ipow(dy, 3)*lens_ipow(lambda, 2) + -0.00169852 *lens_ipow(x, 2)*lens_ipow(dy, 3) + 3.99942 *lens_ipow(dx, 2)*dy*lens_ipow(lambda, 3) + 0.0127049 *lens_ipow(y, 2)*dy*lens_ipow(lambda, 3) + 0.0430536 *lens_ipow(y, 2)*lens_ipow(dx, 2)*lens_ipow(dy, 3) + 0.00229183 *x*y*lens_ipow(dx, 5) + 4.37134e-05 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(dx, 2)*dy + 0.0655604 *y*lens_ipow(dx, 4)*lens_ipow(lambda, 4) + -0.00308816 *lens_ipow(y, 3)*lens_ipow(dy, 6) + -0.00228635 *lens_ipow(x, 2)*dy*lens_ipow(lambda, 6) + -17.2515 *lens_ipow(dy, 5)*lens_ipow(lambda, 5) + 0.00147988 *lens_ipow(y, 3)*lens_ipow(dx, 2)*lens_ipow(lambda, 5) + -3.60928e-14 *lens_ipow(y, 9)*lambda + -9.44344e-14 *x*lens_ipow(y, 8)*dx + -12.1588 *dy*lens_ipow(lambda, 10) + -77.059 *lens_ipow(dy, 3)*lens_ipow(lambda, 8) + -0.24787 *y*lens_ipow(lambda, 10) + 2.87159e-07 *lens_ipow(y, 5)*lens_ipow(lambda, 6) + 0.000274665 *x*lens_ipow(y, 3)*dx*lens_ipow(dy, 2)*lens_ipow(lambda, 4) + 1.22456e-12 *lens_ipow(x, 2)*lens_ipow(y, 7)*lens_ipow(lambda, 2) + 4.54444e-15 *lens_ipow(x, 6)*lens_ipow(y, 5);
const float out_dx =  + -1.85825e-07  + 0.608182 *dx + -3.83536e-08 *y + -0.0102713 *x + 7.47082e-05 *lens_ipow(dx, 2) + 7.60148e-07 *y*dy + 3.7942e-08 *x*dy + -0.534525 *dx*lens_ipow(dy, 2) + -0.513478 *lens_ipow(dx, 3) + -2.51485e-05 *y*lens_ipow(dx, 2) + 0.000139479 *lens_ipow(y, 2)*dx + -0.00216452 *x*lens_ipow(dy, 2) + 0.000337277 *x*y*dy + 3.87528e-06 *x*lens_ipow(y, 2) + 0.000416054 *lens_ipow(x, 2)*dx + 1.17821e-08 *lens_ipow(x, 2)*y + 0.00225279 *x*lens_ipow(lambda, 3) + 6.0162e-06 *lens_ipow(x, 3)*lambda + -0.0112627 *x*lens_ipow(dx, 2)*lens_ipow(lambda, 2) + 9.42778e-10 *x*lens_ipow(y, 4) + 1.44427e-05 *lens_ipow(x, 3)*lens_ipow(dy, 2) + 4.22515e-09 *lens_ipow(x, 3)*lens_ipow(y, 2) + 3.14996e-09 *lens_ipow(x, 5) + 0.242838 *dx*lens_ipow(lambda, 5) + 1.89557e-06 *x*lens_ipow(y, 2)*lens_ipow(lambda, 3) + -4.01037e-05 *lens_ipow(x, 2)*y*dx*dy*lambda + 0.000885214 *x*y*lens_ipow(dy, 5) + -0.000136779 *lens_ipow(x, 3)*lens_ipow(dx, 2)*lens_ipow(dy, 2) + 0.00801133 *x*y*lens_ipow(dx, 2)*lens_ipow(dy, 5) + 0.000654922 *lens_ipow(x, 2)*y*lens_ipow(dx, 5)*dy + -0.000435043 *lens_ipow(x, 3)*lens_ipow(dx, 6) + -9.03874e-05 *lens_ipow(x, 3)*lens_ipow(dx, 2)*lens_ipow(lambda, 5) + -0.901541 *dx*lens_ipow(lambda, 10) + -1.65352 *lens_ipow(dx, 3)*lens_ipow(lambda, 8) + -19.0025 *lens_ipow(dx, 3)*lens_ipow(dy, 2)*lens_ipow(lambda, 6) + -49.1993 *lens_ipow(dx, 7)*lens_ipow(lambda, 4) + 0.00168554 *lens_ipow(y, 2)*dx*lens_ipow(dy, 2)*lens_ipow(lambda, 6) + -0.011703 *x*lens_ipow(lambda, 10) + -6.49905e-05 *lens_ipow(x, 3)*lens_ipow(dy, 2)*lens_ipow(lambda, 6) + -5.26643e-07 *lens_ipow(x, 5)*lens_ipow(dy, 6);
const float out_dy =  + 2.84649e-06  + 0.607813 *dy + -1.87294e-05 *dx + -0.0103448 *y + -4.74824e-07 *x + 1.78648e-06 *y*dy + 3.2244e-06 *y*dx + -0.52495 *lens_ipow(dy, 3) + -0.483809 *lens_ipow(dx, 2)*dy + 0.000545522 *lens_ipow(y, 2)*dy + 4.27867e-06 *lens_ipow(y, 3) + 0.000431408 *x*y*dx + 4.40996e-06 *lens_ipow(x, 2)*y + 0.0024682 *y*lens_ipow(lambda, 3) + 1.34402e-06 *lens_ipow(y, 3)*lambda + -8.36458e-05 *x*lens_ipow(dx, 2)*dy + 0.00023829 *lens_ipow(x, 2)*dy*lambda + 6.31698e-07 *lens_ipow(y, 3)*lens_ipow(dx, 2) + 4.5719e-10 *lens_ipow(y, 5) + -0.000164039 *x*y*dx*lens_ipow(lambda, 2) + 3.20538e-11 *x*lens_ipow(y, 4) + 9.79683e-05 *lens_ipow(x, 2)*lens_ipow(dy, 3) + 1.33789e-09 *lens_ipow(x, 2)*lens_ipow(y, 3) + 1.23125e-07 *lens_ipow(x, 4)*dy + 0.243968 *dy*lens_ipow(lambda, 5) + -2.52177 *lens_ipow(dx, 4)*dy*lens_ipow(lambda, 2) + 1.02018e-08 *lens_ipow(x, 4)*y*lens_ipow(lambda, 2) + 0.00239545 *lens_ipow(y, 2)*lens_ipow(dx, 2)*dy*lens_ipow(lambda, 4) + -9.81465e-05 *lens_ipow(y, 3)*lens_ipow(dy, 6) + 0.0207307 *x*y*lens_ipow(dx, 5)*lens_ipow(dy, 2) + 0.00782527 *x*y*lens_ipow(dx, 7) + -0.000942381 *x*lens_ipow(y, 2)*dx*lens_ipow(dy, 5) + 3.16602e-07 *lens_ipow(y, 4)*dy*lens_ipow(lambda, 5) + -0.955327 *dy*lens_ipow(lambda, 10) + -69.6604 *lens_ipow(dx, 2)*lens_ipow(dy, 5)*lens_ipow(lambda, 4) + -0.0128775 *y*lens_ipow(lambda, 10) + -1.48939e-09 *lens_ipow(y, 5)*dx*lens_ipow(lambda, 5) + 0.00783278 *x*y*lens_ipow(dx, 3)*lens_ipow(lambda, 6) + -0.000645319 *lens_ipow(x, 2)*dy*lens_ipow(lambda, 8) + 7.24717e-17 *lens_ipow(x, 2)*lens_ipow(y, 9);
const float out_transmittance =  + 0.786202  + 0.0996122 *lambda + 9.86266e-06 *dx + 3.70175e-07 *y + -0.00910754 *lens_ipow(dy, 2) + -0.00568612 *lens_ipow(dx, 2) + -1.59751e-05 *lens_ipow(y, 2) + 7.32371e-09 *x*y + -1.56246e-05 *lens_ipow(x, 2) + -0.083723 *lens_ipow(lambda, 3) + -0.000981915 *y*dy*lambda + 1.47383e-05 *y*dx*dy + -0.000779615 *x*dx*lambda + 0.000259471 *lens_ipow(y, 2)*lens_ipow(dy, 2) + 9.2514e-05 *lens_ipow(y, 2)*lens_ipow(dx, 2) + -0.0127392 *x*dx*lens_ipow(dy, 2) + 0.00026543 *lens_ipow(x, 2)*lens_ipow(dx, 2) + -8.78338e-06 *lens_ipow(y, 2)*lens_ipow(lambda, 3) + -8.29356 *lens_ipow(dx, 2)*lens_ipow(dy, 4) + -12.7569 *lens_ipow(dx, 4)*lens_ipow(dy, 2) + -2.42497 *lens_ipow(dx, 6) + -4.74632e-10 *lens_ipow(y, 6) + 0.00746722 *x*y*lens_ipow(dx, 3)*dy + -1.26025e-09 *lens_ipow(x, 2)*lens_ipow(y, 4) + -7.53024e-08 *lens_ipow(x, 4)*lens_ipow(lambda, 2) + -1.20692e-09 *lens_ipow(x, 4)*lens_ipow(y, 2) + -3.8521e-10 *lens_ipow(x, 6) + -0.0533115 *lens_ipow(dy, 2)*lens_ipow(lambda, 5) + -1.56198e-09 *lens_ipow(x, 2)*lens_ipow(y, 3)*lens_ipow(dy, 2) + -20.866 *lens_ipow(dy, 8) + -0.0641845 *lens_ipow(dx, 2)*lens_ipow(lambda, 6) + 2.73869e-05 *x*lens_ipow(y, 3)*dx*lens_ipow(dy, 3) + -0.0713101 *lens_ipow(y, 2)*lens_ipow(dx, 8) + 1.82671e-08 *lens_ipow(y, 6)*lens_ipow(dy, 4) + 3.85362e-05 *lens_ipow(x, 4)*lens_ipow(dx, 6) + 1.07602e-07 *lens_ipow(x, 4)*lens_ipow(y, 2)*lens_ipow(dx, 2)*lens_ipow(dy, 2) + 1.53421e-08 *lens_ipow(x, 6)*lens_ipow(dx, 2)*lens_ipow(dy, 2) + 0.160036 *lens_ipow(lambda, 11) + 1.68215e-06 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(dy, 2)*lens_ipow(lambda, 5) + -2.90093e-08 *lens_ipow(x, 5)*dx*lens_ipow(lambda, 5);