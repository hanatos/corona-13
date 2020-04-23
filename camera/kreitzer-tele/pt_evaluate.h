const float out_x =  + -7.53719e-05  + -0.0040198 *dy + 399.583 *dx + -2.28631e-05 *y + 2.58883 *x + 0.0282439 *lens_ipow(dx, 2) + 2.41907e-06 *x*y + 5.21373e-06 *lens_ipow(x, 2) + 1.70625 *dx*lens_ipow(lambda, 2) + -241.438 *dx*lens_ipow(dy, 2) + -245.154 *lens_ipow(dx, 3) + -1.18733 *y*dx*dy + 0.0114028 *lens_ipow(y, 2)*dx + 1.26004 *x*lens_ipow(dy, 2) + 0.0414114 *x*y*dy + 0.000350368 *x*lens_ipow(y, 2) + 0.0576284 *lens_ipow(x, 2)*dx + 0.000350798 *lens_ipow(x, 3) + -0.0874267 *x*lens_ipow(lambda, 3) + 0.0227785 *x*lens_ipow(dx, 2)*lambda + 0.0131072 *x*y*dy*lambda + 2.20274e-06 *x*lens_ipow(y, 2)*dy + 343.293 *lens_ipow(dx, 5) + -3.40182e-05 *lens_ipow(x, 2)*lens_ipow(lambda, 3) + 0.00275759 *lens_ipow(x, 2)*lens_ipow(dy, 3)*lambda + -106.895 *lens_ipow(dx, 3)*lens_ipow(lambda, 4) + 29.8208 *y*lens_ipow(dx, 3)*dy*lens_ipow(lambda, 2) + 1.00066 *x*lens_ipow(dy, 2)*lens_ipow(lambda, 4) + -2.10546e-06 *lens_ipow(x, 2)*lens_ipow(y, 3)*dx*dy + 0.0177254 *lens_ipow(x, 3)*lens_ipow(dx, 4) + 0.000103079 *x*lens_ipow(y, 2)*lens_ipow(lambda, 5) + 3.66702e-05 *lens_ipow(x, 4)*dx*lens_ipow(lambda, 3) + -2.0833 *lens_ipow(y, 2)*dx*lens_ipow(dy, 4)*lens_ipow(lambda, 2) + -6.11315 *dx*lens_ipow(lambda, 9) + 0.277279 *x*lens_ipow(lambda, 9) + 2943.8 *lens_ipow(dx, 2)*lens_ipow(dy, 4)*lens_ipow(lambda, 5) + 9.99168e+06 *lens_ipow(dx, 9)*lens_ipow(dy, 2) + 0.719103 *lens_ipow(y, 2)*lens_ipow(dx, 3)*lens_ipow(lambda, 6) + 4.37915e-07 *lens_ipow(x, 5)*lens_ipow(lambda, 6) + 6.93607e-16 *lens_ipow(x, 7)*lens_ipow(y, 4);
const float out_y =  + -6.14701e-05  + 398.883 *dy + 2.58684 *y + 1.48627e-05 *x + 0.000427641 *y*dy + -2.87632e-06 *x*y + 1.79306e-06 *lens_ipow(x, 2) + 5.16056 *dy*lens_ipow(lambda, 2) + -241.665 *lens_ipow(dy, 3) + -290.011 *lens_ipow(dx, 2)*dy + 0.413533 *y*lens_ipow(dx, 2) + 0.057092 *lens_ipow(y, 2)*dy + 0.000349051 *lens_ipow(y, 3) + -3.05602 *x*dx*dy + 5.3719e-05 *x*y*dy + 0.00877699 *lens_ipow(x, 2)*dy + 0.000301458 *lens_ipow(x, 2)*y + -0.0599654 *y*lens_ipow(lambda, 3) + 0.0330832 *y*lens_ipow(dy, 2)*lambda + -0.019729 *y*lens_ipow(dx, 3) + 0.0821105 *x*y*dx*lambda + 199.497 *lens_ipow(dy, 5) + 5.50017 *x*dx*dy*lens_ipow(lambda, 2) + -0.00051614 *x*lens_ipow(y, 2)*dx*dy + 5.21794e-06 *lens_ipow(x, 3)*y*dx + -128.211 *lens_ipow(dy, 3)*lens_ipow(lambda, 4) + 7.10664 *y*lens_ipow(dx, 2)*lens_ipow(lambda, 4) + 160.905 *y*lens_ipow(dx, 2)*lens_ipow(dy, 4) + 52.3091 *y*lens_ipow(dx, 6) + -0.000133896 *lens_ipow(y, 4)*lens_ipow(dy, 3) + -0.0104561 *x*y*lens_ipow(dy, 3)*lens_ipow(lambda, 2) + 0.000427871 *lens_ipow(x, 2)*y*lens_ipow(lambda, 4) + 6.65948e-11 *lens_ipow(x, 6)*y + 4.79792e-05 *lens_ipow(y, 4)*dy*lens_ipow(lambda, 3) + 2.88671e-09 *lens_ipow(x, 2)*lens_ipow(y, 3)*lens_ipow(lambda, 4) + -40.335 *dy*lens_ipow(lambda, 9) + 1958.2 *lens_ipow(dx, 2)*dy*lens_ipow(lambda, 8) + 5.94329e-07 *lens_ipow(y, 5)*lens_ipow(lambda, 6) + 0.136651 *lens_ipow(x, 2)*dy*lens_ipow(lambda, 8) + 4.95212e-16 *lens_ipow(x, 2)*lens_ipow(y, 9);
const float out_dx =  + -1.36984e-06  + -3.02386 *dx + -0.0220831 *x + -12.6295 *dx*lens_ipow(dy, 2) + 1.85207 *lens_ipow(dx, 3) + -0.181349 *y*dx*dy + -0.00070382 *lens_ipow(y, 2)*dx + -0.105079 *x*lens_ipow(dy, 2) + 0.00846378 *x*lens_ipow(dx, 2) + -0.0016196 *x*y*dy + -6.67979e-06 *x*lens_ipow(y, 2) + -0.000318603 *lens_ipow(x, 2)*dx + -2.13191e-06 *lens_ipow(x, 3) + 0.000778602 *x*lens_ipow(lambda, 3) + -0.298228 *lens_ipow(dx, 5) + -0.000632288 *y*lens_ipow(dx, 3)*dy + 8.03208e-09 *lens_ipow(x, 3)*lens_ipow(y, 2)*lens_ipow(dx, 4) + -0.00411936 *x*lens_ipow(lambda, 10);
const float out_dy =  + -1.39664e-06  + -3.02497 *dy + -5.30828e-06 *dx + -0.0220656 *y + 0.000282626 *lens_ipow(dy, 2) + 2.79523e-06 *y*dx + -2.11973e-06 *x*dy + 1.89904 *lens_ipow(dy, 3) + 16.361 *lens_ipow(dx, 2)*dy + 0.00984536 *y*lens_ipow(dy, 2) + 0.0845495 *y*lens_ipow(dx, 2) + -0.000303261 *lens_ipow(y, 2)*dy + -2.07305e-06 *lens_ipow(y, 3) + 1.20882e-06 *x*lens_ipow(dy, 2) + 0.223632 *x*dx*dy + 0.00106874 *x*y*dx + 0.000690967 *lens_ipow(x, 2)*dy + 2.47702e-06 *lens_ipow(x, 2)*y + 5.98681e-07 *lens_ipow(y, 2)*lens_ipow(dx, 2) + -0.00870921 *x*dx*dy*lambda + -0.000127605 *x*y*dx*lambda + 0.00111463 *y*lens_ipow(lambda, 4) + 0.0171457 *x*dx*lens_ipow(dy, 3) + 2.16915e-09 *x*lens_ipow(y, 3)*dx + -0.840573 *lens_ipow(dx, 2)*dy*lens_ipow(lambda, 3) + 3.8148 *lens_ipow(dx, 4)*lens_ipow(dy, 3) + -0.0128798 *y*lens_ipow(dx, 2)*lens_ipow(lambda, 4) + 0.0016008 *y*lens_ipow(dx, 6) + -6.72972e-05 *lens_ipow(x, 2)*dy*lens_ipow(lambda, 4) + -2.33936 *lens_ipow(dy, 5)*lens_ipow(lambda, 4) + 20.5479 *lens_ipow(dx, 2)*lens_ipow(dy, 3)*lens_ipow(lambda, 4) + -1.4105e-06 *lens_ipow(x, 2)*y*lens_ipow(lambda, 6) + -0.00165138 *dy*lens_ipow(lambda, 9) + 12.7052 *lens_ipow(dx, 4)*dy*lens_ipow(lambda, 5) + -0.00501469 *y*lens_ipow(lambda, 10) + 0.000271308 *lens_ipow(y, 3)*lens_ipow(dx, 4)*lens_ipow(lambda, 4);
const float out_transmittance =  + 0.416426  + 0.271139 *lambda + 1.25411e-05 *dy + -9.02847e-05 *dx + 2.92475e-07 *y + -4.35994e-07 *x + -0.143453 *lens_ipow(dy, 2) + -0.134492 *lens_ipow(dx, 2) + -0.00255901 *y*dy + -7.55823e-06 *y*dx + -1.1374e-05 *lens_ipow(y, 2) + -0.00243245 *x*dx + -1.16659e-05 *lens_ipow(x, 2) + -0.225852 *lens_ipow(lambda, 3) + 0.00454149 *lens_ipow(dx, 3) + 6.74301e-05 *y*lens_ipow(dx, 2) + -4.28026e-09 *lens_ipow(x, 2)*y + 0.808027 *lens_ipow(dx, 2)*lens_ipow(dy, 2) + 0.0666357 *y*lens_ipow(dy, 3) + -0.000179446 *lens_ipow(y, 2)*lens_ipow(dx, 2) + -4.81253e-06 *lens_ipow(y, 3)*dy + 0.000739319 *x*lens_ipow(dx, 2)*dy + 0.0388081 *x*lens_ipow(dx, 3) + -5.01914e-06 *x*lens_ipow(y, 2)*dx + 0.000312941 *lens_ipow(x, 2)*lens_ipow(dy, 2) + -2.94753e-08 *lens_ipow(x, 2)*lens_ipow(y, 2) + -3.96148e-06 *lens_ipow(x, 3)*dx + 0.105056 *x*dx*lens_ipow(dy, 2)*lambda + 134.378 *lens_ipow(dy, 6) + 38.0779 *lens_ipow(dx, 6) + -6.74165e-11 *lens_ipow(y, 6) + -4.71779e-11 *lens_ipow(x, 6) + -8.04499e-06 *lens_ipow(x, 2)*y*dy*lens_ipow(lambda, 5) + 514.727 *y*lens_ipow(dy, 9) + -3.58428e-14 *lens_ipow(y, 9)*dy + -0.000707711 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(dx, 4)*lens_ipow(dy, 2) + -0.000128013 *lens_ipow(x, 4)*lens_ipow(dy, 6) + -2.22617e-14 *lens_ipow(x, 8)*y*dy + 0.404952 *lens_ipow(lambda, 11) + 15659.6 *lens_ipow(dx, 4)*lens_ipow(dy, 4)*lens_ipow(lambda, 3);