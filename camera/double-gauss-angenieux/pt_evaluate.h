const float out_x =  + 49.7564 *dx + -0.45111 *x + -20.6769 *dx*lens_ipow(dy, 2) + -21.6629 *lens_ipow(dx, 3) + 0.158343 *y*dx*dy + 0.00540493 *lens_ipow(y, 2)*dx + 0.222106 *x*lens_ipow(dy, 2) + 0.325749 *x*lens_ipow(dx, 2) + 0.00588602 *x*y*dy + -0.000925557 *x*lens_ipow(y, 2) + -5.36279e-05 *lens_ipow(x, 2)*dy + 0.00983509 *lens_ipow(x, 2)*dx + -0.00105935 *lens_ipow(x, 3) + 0.327477 *x*lens_ipow(lambda, 3) + 1.70099 *x*lens_ipow(dx, 2)*lens_ipow(dy, 2) + 0.121221 *x*y*lens_ipow(dx, 2)*dy + 0.0633018 *lens_ipow(x, 2)*lens_ipow(dx, 3) + -0.000739957 *lens_ipow(x, 2)*y*dx*dy + -2.39731e-06 *lens_ipow(x, 3)*lens_ipow(y, 2) + 1.18041 *dx*lens_ipow(lambda, 5) + 110.698 *lens_ipow(dx, 3)*lens_ipow(dy, 4) + 96.7421 *lens_ipow(dx, 5)*lens_ipow(dy, 2) + 9.06137 *y*lens_ipow(dx, 3)*lens_ipow(dy, 3) + 0.000325899 *lens_ipow(y, 4)*dx*lens_ipow(dy, 2) + 4.44504 *x*lens_ipow(dx, 6) + -2.32777e-06 *x*lens_ipow(y, 4)*lens_ipow(dx, 2) + -4.70582e-09 *x*lens_ipow(y, 6) + 0.128425 *lens_ipow(x, 2)*dx*lens_ipow(dy, 4) + -1.07899e-07 *lens_ipow(x, 3)*lens_ipow(y, 3)*dy + 135.503 *lens_ipow(dx, 9) + 0.182123 *x*y*lens_ipow(dy, 7) + -1.3638 *x*y*lens_ipow(dx, 6)*dy + 2.81701e-08 *lens_ipow(x, 6)*y*dx*dy + -8.89742e-12 *lens_ipow(x, 9) + 5.31624e-06 *lens_ipow(x, 3)*lens_ipow(y, 2)*lens_ipow(lambda, 5) + 2.12111e-06 *lens_ipow(x, 5)*lens_ipow(lambda, 5) + -2.09036 *x*lens_ipow(lambda, 10) + 7.36594e-09 *x*lens_ipow(y, 6)*lens_ipow(lambda, 4) + 0.0427966 *lens_ipow(x, 3)*lens_ipow(dy, 8) + -1.69755e-13 *lens_ipow(x, 7)*lens_ipow(y, 4);
const float out_y =  + -0.000255802  + 49.7076 *dy + -0.450998 *y + -21.7135 *lens_ipow(dy, 3) + -20.4721 *lens_ipow(dx, 2)*dy + 0.302858 *y*lens_ipow(dy, 2) + 0.213691 *y*lens_ipow(dx, 2) + 0.0108702 *lens_ipow(y, 2)*dy + -0.00104189 *lens_ipow(y, 3) + 0.0158971 *x*dx*dy + 0.00801017 *x*y*dx + 0.0057819 *lens_ipow(x, 2)*dy + -0.000902568 *lens_ipow(x, 2)*y + 0.324103 *y*lens_ipow(lambda, 3) + 2.21888 *y*lens_ipow(dx, 2)*lens_ipow(dy, 2) + 0.069591 *lens_ipow(y, 2)*lens_ipow(dy, 3) + 0.0258352 *lens_ipow(y, 2)*lens_ipow(dx, 2)*dy + 1.55084 *x*dx*lens_ipow(dy, 3) + 0.908836 *x*lens_ipow(dx, 3)*dy + 0.120407 *x*y*dx*lens_ipow(dy, 2) + -2.19662e-05 *x*lens_ipow(y, 3)*dx + -1.73416e-05 *lens_ipow(x, 2)*lens_ipow(y, 2)*dy + -2.89843e-06 *lens_ipow(x, 2)*lens_ipow(y, 3) + 1.24163 *dy*lens_ipow(lambda, 5) + 107.925 *lens_ipow(dx, 2)*lens_ipow(dy, 5) + 6.00745 *y*lens_ipow(dy, 6) + 0.259704 *lens_ipow(x, 2)*lens_ipow(dx, 2)*lens_ipow(dy, 3) + 0.000217323 *lens_ipow(x, 4)*lens_ipow(dx, 2)*dy + -2.49484e-06 *lens_ipow(x, 4)*y*lens_ipow(dy, 2) + -4.30423e-09 *lens_ipow(x, 6)*y + -7.24966e-08 *lens_ipow(y, 6)*dy*lambda + -1.60675e-07 *lens_ipow(x, 5)*y*dx*lambda + 180.736 *lens_ipow(dy, 9) + 559.35 *lens_ipow(dx, 6)*lens_ipow(dy, 3) + 0.00652112 *lens_ipow(y, 3)*lens_ipow(dx, 6) + -9.80956e-12 *lens_ipow(y, 9) + 1.55747e-06 *lens_ipow(y, 5)*lens_ipow(lambda, 5) + 7.243e-06 *lens_ipow(x, 2)*lens_ipow(y, 3)*lens_ipow(lambda, 5) + -2.0189 *y*lens_ipow(lambda, 10) + -1.44756e-13 *lens_ipow(x, 4)*lens_ipow(y, 7);
const float out_dx =  + -0.606975 *dx + -0.0145462 *x + 0.139927 *dx*lens_ipow(dy, 2) + 0.230118 *lens_ipow(dx, 3) + -2.07784e-05 *lens_ipow(y, 2)*dx + -0.0020595 *x*lens_ipow(dy, 2) + -7.26888e-05 *x*y*dy + 1.41104e-05 *x*lens_ipow(y, 2) + 8.55316e-07 *lens_ipow(x, 2)*dy + -0.000173851 *lens_ipow(x, 2)*dx + 1.48755e-05 *lens_ipow(x, 3) + 0.0083638 *dx*lens_ipow(lambda, 3) + -0.00396336 *x*lens_ipow(lambda, 3) + -0.00373201 *y*dx*lens_ipow(dy, 3) + 0.00050211 *lens_ipow(y, 2)*lens_ipow(dx, 3) + 9.26117e-06 *lens_ipow(y, 3)*dx*dy + 9.97076e-06 *x*lens_ipow(y, 2)*lens_ipow(dy, 2) + 0.000797341 *lens_ipow(x, 2)*dx*lens_ipow(dy, 2) + 0.0010059 *lens_ipow(x, 2)*lens_ipow(dx, 3) + 4.51109e-05 *lens_ipow(x, 2)*y*dx*dy + 0.000397735 *lens_ipow(y, 2)*dx*lens_ipow(dy, 2)*lambda + 0.00300323 *x*y*lens_ipow(dx, 4)*dy + 5.05886e-11 *x*lens_ipow(y, 6) + -1.96346e-09 *lens_ipow(x, 3)*lens_ipow(y, 3)*dy + 1.1502e-08 *lens_ipow(x, 5)*lens_ipow(dy, 2) + 1.40362e-10 *lens_ipow(x, 5)*lens_ipow(y, 2) + 0.0197079 *x*y*lens_ipow(dx, 2)*lens_ipow(dy, 5) + -3.53467e-08 *x*lens_ipow(y, 4)*lens_ipow(lambda, 4) + 1.18601e-05 *lens_ipow(x, 3)*y*lens_ipow(dy, 5) + -3.43286e-08 *lens_ipow(x, 3)*lens_ipow(y, 3)*lens_ipow(dx, 2)*dy + -1.33799e-09 *lens_ipow(x, 4)*lens_ipow(y, 3)*dx*dy + -1.23413e-11 *lens_ipow(x, 4)*lens_ipow(y, 4)*dx + -2.34649e-08 *lens_ipow(x, 5)*lens_ipow(lambda, 4) + 6.38443e-07 *lens_ipow(x, 5)*lens_ipow(dx, 4) + -8.8726e-14 *lens_ipow(x, 8)*dx + 1.06211e-13 *lens_ipow(x, 9) + -0.0100528 *lens_ipow(y, 2)*lens_ipow(dx, 9) + 0.0267101 *x*lens_ipow(lambda, 10) + -2.41145e-10 *lens_ipow(x, 5)*lens_ipow(y, 2)*lens_ipow(lambda, 4) + 1.90932e-15 *lens_ipow(x, 5)*lens_ipow(y, 6);
const float out_dy =  + -0.60801 *dy + -0.0145595 *y + 0.245273 *lens_ipow(dy, 3) + 0.34091 *lens_ipow(dx, 2)*dy + -0.00316615 *y*lens_ipow(dx, 2) + -0.000146289 *lens_ipow(y, 2)*dy + 1.46753e-05 *lens_ipow(y, 3) + 0.00366859 *x*dx*dy + 1.2712e-06 *x*y*dy + -0.000104641 *x*y*dx + -7.15017e-05 *lens_ipow(x, 2)*dy + 1.44672e-05 *lens_ipow(x, 2)*y + 0.009071 *dy*lens_ipow(lambda, 3) + -0.00396411 *y*lens_ipow(lambda, 3) + -0.00494903 *y*lens_ipow(dy, 4) + 0.000699537 *lens_ipow(y, 2)*lens_ipow(dy, 3) + 0.000640116 *lens_ipow(y, 2)*lens_ipow(dx, 2)*dy + 1.61631e-05 *lens_ipow(y, 3)*lens_ipow(dy, 2) + -5.75467e-06 *lens_ipow(y, 3)*lens_ipow(dx, 2) + 1.98884e-05 *x*lens_ipow(y, 2)*dx*dy + 0.000390292 *lens_ipow(x, 2)*lens_ipow(dy, 3) + 0.000513268 *lens_ipow(x, 2)*lens_ipow(dx, 2)*dy + -1.26958 *lens_ipow(dx, 2)*lens_ipow(dy, 5) + -0.552503 *lens_ipow(dx, 6)*dy + -5.22359e-05 *lens_ipow(y, 3)*lens_ipow(dx, 2)*lens_ipow(dy, 2) + 0.000773699 *x*y*lens_ipow(dx, 5) + -4.78816e-07 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(dx, 2)*dy + 1.30194e-10 *lens_ipow(x, 2)*lens_ipow(y, 5) + -5.98885e-10 *lens_ipow(x, 4)*lens_ipow(y, 2)*dy + 1.05778e-10 *lens_ipow(x, 4)*lens_ipow(y, 3) + 4.31431e-11 *lens_ipow(x, 6)*y + -0.943972 *lens_ipow(dy, 9) + 1.3409e-13 *lens_ipow(y, 8)*dy + 1.12021e-13 *lens_ipow(y, 9) + 0.0370831 *x*y*lens_ipow(dx, 3)*lens_ipow(dy, 4) + 6.01526e-08 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(dy, 3)*lens_ipow(lambda, 2) + -2.28646e-08 *lens_ipow(x, 4)*y*lens_ipow(lambda, 4) + -3.25636e-08 *lens_ipow(y, 5)*lens_ipow(lambda, 5) + -7.83459e-08 *lens_ipow(x, 2)*lens_ipow(y, 3)*lens_ipow(lambda, 5) + 0.0273723 *y*lens_ipow(lambda, 10);
const float out_transmittance =  + 0.333512  + 0.280781 *lambda + -2.84557e-06 *dy + -1.76181e-06 *dx + -0.027121 *lens_ipow(dy, 2) + -0.0321791 *lens_ipow(dx, 2) + -0.000337426 *y*dy + -6.6444e-06 *lens_ipow(y, 2) + -0.000388471 *x*dx + -1.35534e-05 *lens_ipow(x, 2) + -0.232101 *lens_ipow(lambda, 3) + -0.00688131 *y*lens_ipow(dy, 3) + -0.00592774 *y*lens_ipow(dx, 2)*dy + -0.000324578 *lens_ipow(y, 2)*lens_ipow(dy, 2) + -7.56015e-07 *lens_ipow(y, 3)*dy + -0.00595614 *x*dx*lens_ipow(dy, 2) + -0.00731767 *x*lens_ipow(dx, 3) + -0.000428209 *x*y*dx*dy + 1.60803e-06 *lens_ipow(x, 3)*dx + -0.554329 *lens_ipow(dy, 6) + -1.60078 *lens_ipow(dx, 2)*lens_ipow(dy, 4) + -1.40849 *lens_ipow(dx, 4)*lens_ipow(dy, 2) + -0.435074 *lens_ipow(dx, 6) + -0.000536006 *lens_ipow(y, 2)*lens_ipow(dx, 4) + -1.00769e-10 *lens_ipow(y, 6) + -1.10345e-09 *x*lens_ipow(y, 4)*dx + -0.000422883 *lens_ipow(x, 2)*lens_ipow(dy, 4) + -0.00174843 *lens_ipow(x, 2)*lens_ipow(dx, 4) + 1.33591e-05 *lens_ipow(x, 2)*y*lens_ipow(dy, 3) + -1.19873e-08 *lens_ipow(x, 2)*lens_ipow(y, 3)*dy + -6.00716e-09 *lens_ipow(x, 3)*lens_ipow(y, 2)*dx + -6.69677e-10 *lens_ipow(y, 6)*lens_ipow(dx, 2) + -1.9219e-05 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(dx, 2)*lens_ipow(dy, 2) + -2.16786e-12 *lens_ipow(x, 2)*lens_ipow(y, 6) + 7.41756e-06 *lens_ipow(x, 3)*y*dx*lens_ipow(dy, 3) + 4.20375e-12 *lens_ipow(x, 6)*y*dy + -2.30144e-12 *lens_ipow(x, 6)*lens_ipow(y, 2) + -0.749089 *lens_ipow(dx, 2)*lens_ipow(dy, 2)*lens_ipow(lambda, 5) + 4.60505e-10 *lens_ipow(x, 7)*dx*lens_ipow(dy, 2) + 0.420019 *lens_ipow(lambda, 11);
