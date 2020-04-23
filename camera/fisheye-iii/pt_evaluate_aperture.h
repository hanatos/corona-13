const float out_x =  + 0.000902259  + 37.6368 *dx + 0.116556 *x + -9.58971e-06 *lens_ipow(y, 2) + -5.62534e-06 *lens_ipow(x, 2) + -0.153206 *dx*lens_ipow(dy, 2) + 1.17228 *lens_ipow(dx, 3) + 0.838948 *y*dx*dy + 0.00613825 *lens_ipow(y, 2)*dx + -0.0167978 *x*lens_ipow(dy, 2) + 0.889901 *x*lens_ipow(dx, 2) + 0.0112782 *x*y*dy + 0.0171549 *lens_ipow(x, 2)*dx + -5.31713e-06 *x*lens_ipow(y, 2)*dx + -7.4088 *dx*lens_ipow(lambda, 4) + -0.0160046 *lens_ipow(y, 2)*lens_ipow(dx, 3) + -0.234674 *x*lens_ipow(lambda, 4) + 2.9518e-07 *x*lens_ipow(y, 4) + 0.000797167 *lens_ipow(x, 2)*y*dx*dy + 3.07656e-05 *lens_ipow(x, 2)*lens_ipow(y, 2)*dx + 6.9492e-07 *lens_ipow(x, 3)*lens_ipow(y, 2) + 2.97897e-07 *lens_ipow(x, 5) + 386.299 *lens_ipow(dx, 5)*lens_ipow(dy, 2) + 8.74808e-07 *lens_ipow(y, 5)*dx*dy + 1.5698e-08 *lens_ipow(y, 6)*dx + -0.0827476 *x*y*lens_ipow(dy, 5) + 0.00203023 *lens_ipow(x, 3)*lens_ipow(dx, 2)*lens_ipow(dy, 2) + 8.74489e-08 *lens_ipow(x, 3)*lens_ipow(y, 3)*dy + -0.00495765 *x*y*dy*lens_ipow(lambda, 5) + -0.00023624 *lens_ipow(x, 2)*lens_ipow(y, 2)*dx*lens_ipow(dy, 2)*lambda + -0.000116432 *lens_ipow(x, 4)*lens_ipow(dx, 3)*lambda + 538.864 *dx*lens_ipow(dy, 8) + 2377.24 *lens_ipow(dx, 3)*lens_ipow(dy, 6) + 427.056 *lens_ipow(dx, 9) + 7.6351e-11 *x*lens_ipow(y, 7)*dy + 1.0104e-10 *lens_ipow(x, 7)*y*dy + 1.41572e-10 *lens_ipow(x, 8)*dx + 35.8328 *dx*lens_ipow(lambda, 10) + 1.14029 *x*lens_ipow(lambda, 10) + -3.68189e-09 *lens_ipow(x, 5)*lens_ipow(y, 2)*lens_ipow(lambda, 4);
const float out_y =  + -0.000532841  + 37.6385 *dy + 0.114645 *y + 0.772044 *lens_ipow(dy, 3) + 1.72184 *lens_ipow(dx, 2)*dy + 0.885817 *y*lens_ipow(dy, 2) + 0.0125457 *y*lens_ipow(dx, 2) + 0.0164211 *lens_ipow(y, 2)*dy + 9.08096e-05 *lens_ipow(y, 3) + 0.918219 *x*dx*dy + 0.0108813 *x*y*dx + 0.00652653 *lens_ipow(x, 2)*dy + 0.000100077 *lens_ipow(x, 2)*y + -0.152055 *y*lens_ipow(lambda, 3) + -7.44186 *dy*lens_ipow(lambda, 4) + -0.000281488 *lens_ipow(y, 3)*lens_ipow(dx, 2) + 2.47932e-05 *lens_ipow(x, 2)*lens_ipow(y, 2)*dy + 1.57224e-05 *lens_ipow(x, 3)*y*dx + -0.0836668 *lens_ipow(y, 2)*lens_ipow(dx, 2)*dy*lens_ipow(lambda, 2) + 5.31678e-08 *lens_ipow(y, 6)*dy + 0.00387435 *x*lens_ipow(y, 2)*dx*lens_ipow(dy, 3) + 0.00554742 *x*lens_ipow(y, 2)*lens_ipow(dx, 3)*dy + -0.00015151 *lens_ipow(x, 2)*y*lens_ipow(lambda, 4) + 0.00213603 *lens_ipow(x, 2)*y*lens_ipow(dx, 4) + 1.40558e-09 *lens_ipow(x, 2)*lens_ipow(y, 5) + -3.93158e-05 *lens_ipow(x, 4)*lens_ipow(dx, 2)*dy + 5.45554e-10 *lens_ipow(x, 6)*y + -0.474885 *x*dx*dy*lens_ipow(lambda, 5) + -0.0152121 *x*y*dx*lens_ipow(lambda, 5) + -0.0482852 *lens_ipow(x, 2)*lens_ipow(dy, 3)*lens_ipow(lambda, 3) + 590.433 *lens_ipow(dy, 9) + 8013.18 *lens_ipow(dx, 4)*lens_ipow(dy, 5) + 786.742 *lens_ipow(dx, 8)*dy + 1.14028e-05 *lens_ipow(y, 5)*lens_ipow(dy, 4) + 1.25035e-07 *x*lens_ipow(y, 5)*dx*lens_ipow(lambda, 2) + 5.36908e-06 *lens_ipow(x, 4)*y*lens_ipow(dy, 4) + 3.66262e-11 *lens_ipow(x, 8)*dy + 35.8743 *dy*lens_ipow(lambda, 10) + 0.810168 *y*lens_ipow(lambda, 10) + 2.93762e-15 *lens_ipow(y, 11);
const float out_dx =  + -9.58194e-06  + -0.667314 *dx + -0.0285989 *x + -1.44854 *dx*lens_ipow(dy, 2) + -1.10292 *lens_ipow(dx, 3) + -0.0169563 *y*dx*dy + -0.00050775 *lens_ipow(y, 2)*dx + 0.000398916 *x*lens_ipow(dy, 2) + -0.0109949 *x*lens_ipow(dx, 2) + -0.000213514 *x*y*dy + -9.9664e-06 *x*lens_ipow(y, 2) + -0.000699554 *lens_ipow(x, 2)*dx + -9.96278e-06 *lens_ipow(x, 3) + -0.236851 *dx*lens_ipow(lambda, 3) + -1.3038 *lens_ipow(dx, 3)*lens_ipow(lambda, 2) + 0.00112127 *x*y*lens_ipow(dx, 2)*dy + 7.95606e-06 *lens_ipow(x, 3)*lens_ipow(dy, 2) + -0.000508211 *x*lens_ipow(lambda, 5) + 0.00444486 *lens_ipow(y, 2)*dx*lens_ipow(dy, 4) + -0.0870952 *x*lens_ipow(dx, 2)*lens_ipow(lambda, 4) + 5.53917e-05 *x*lens_ipow(y, 2)*lens_ipow(dx, 4) + -3.15492e-06 *x*lens_ipow(y, 3)*lens_ipow(dy, 3) + -6.08308e-12 *x*lens_ipow(y, 6) + 1.00772e-07 *lens_ipow(x, 2)*lens_ipow(y, 3)*dx*dy + 4.73016e-05 *lens_ipow(x, 3)*lens_ipow(dx, 2)*lens_ipow(dy, 2) + -0.0566703 *y*dx*dy*lens_ipow(lambda, 5) + 6.75305e-07 *lens_ipow(x, 4)*y*dx*lens_ipow(dy, 3) + 1.10954e-11 *lens_ipow(x, 7)*lens_ipow(dx, 2) + -6.00966e-14 *lens_ipow(x, 7)*lens_ipow(y, 2) + -16.5736 *dx*lens_ipow(dy, 4)*lens_ipow(lambda, 5) + -1.00218e-05 *lens_ipow(y, 4)*lens_ipow(dx, 3)*lens_ipow(lambda, 3) + -0.195445 *x*lens_ipow(dy, 4)*lens_ipow(lambda, 5) + -8.96731e-11 *lens_ipow(x, 3)*lens_ipow(y, 4)*lens_ipow(lambda, 3) + -6.67883e-09 *lens_ipow(x, 4)*lens_ipow(y, 2)*dx*lens_ipow(lambda, 3) + 1.11761 *dx*lens_ipow(lambda, 10) + -53.5004 *lens_ipow(dx, 5)*lens_ipow(dy, 2)*lens_ipow(lambda, 4) + -0.00247517 *x*y*dy*lens_ipow(lambda, 8) + 2.40364e-05 *x*lens_ipow(y, 3)*lens_ipow(dx, 2)*dy*lens_ipow(lambda, 4) + -7.79088e-06 *lens_ipow(x, 4)*dx*lens_ipow(lambda, 6) + -6.7343e-11 *lens_ipow(x, 7)*lens_ipow(lambda, 4);
const float out_dy =  + 1.25512e-06  + -0.667848 *dy + -0.0285768 *y + -0.001593 *lens_ipow(dy, 2) + -3.91201e-05 *y*dy + -1.96027e-05 *y*dx + -1.08648 *lens_ipow(dy, 3) + -1.47453 *lens_ipow(dx, 2)*dy + -0.0108042 *y*lens_ipow(dy, 2) + -0.000295841 *y*lens_ipow(dx, 2) + -0.0006781 *lens_ipow(y, 2)*dy + -9.73318e-06 *lens_ipow(y, 3) + -0.0183763 *x*dx*dy + -0.000251102 *x*y*dx + -0.000507707 *lens_ipow(x, 2)*dy + -9.87308e-06 *lens_ipow(x, 2)*y + -0.240532 *dy*lens_ipow(lambda, 3) + -1.07008e-05 *x*y*lens_ipow(dx, 2) + -1.38878 *lens_ipow(dy, 3)*lens_ipow(lambda, 2) + -0.000665779 *y*lens_ipow(lambda, 4) + 0.000927268 *x*y*dx*lens_ipow(dy, 2) + -0.0902573 *y*lens_ipow(dy, 2)*lens_ipow(lambda, 4) + -3.88614e-06 *lens_ipow(y, 4)*lens_ipow(dx, 2)*dy + 8.22724e-05 *lens_ipow(x, 2)*y*lens_ipow(dx, 4) + -6.33144e-07 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(dy, 3) + -1.92613e-06 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(dx, 2)*dy + -3.20931e-11 *lens_ipow(x, 2)*lens_ipow(y, 5) + -6.89377e-05 *lens_ipow(x, 3)*lens_ipow(dx, 3)*dy + -7.86867e-10 *lens_ipow(x, 3)*lens_ipow(y, 3)*dx + -8.46515e-12 *lens_ipow(x, 6)*y + -0.0313921 *x*dx*dy*lens_ipow(lambda, 5) + -3.43879e-06 *lens_ipow(x, 4)*lens_ipow(dy, 3)*lambda + -6.24534e-09 *lens_ipow(x, 2)*lens_ipow(y, 4)*dy*lens_ipow(lambda, 2) + -7.71813 *lens_ipow(dx, 4)*dy*lens_ipow(lambda, 5) + -4.00997e-05 *lens_ipow(x, 2)*y*lens_ipow(dy, 2)*lens_ipow(lambda, 5) + 1.50561 *dy*lens_ipow(lambda, 10) + -0.00461612 *lens_ipow(y, 2)*dy*lens_ipow(lambda, 8) + -2.06316e-08 *lens_ipow(y, 6)*lens_ipow(dy, 5) + -4.58061e-11 *lens_ipow(y, 7)*lens_ipow(lambda, 4) + -1.37717e-06 *x*lens_ipow(y, 3)*dx*lens_ipow(lambda, 6);
const float out_transmittance =  + 0.51491  + 0.27353 *lambda + 2.19793e-05 *dy + 5.47844e-07 *y + -0.0241163 *lens_ipow(dy, 2) + -0.0277282 *lens_ipow(dx, 2) + -0.000376043 *y*dy + -7.06492e-07 *y*dx + -3.45023e-06 *lens_ipow(y, 2) + -0.000396702 *x*dx + -3.35017e-06 *lens_ipow(x, 2) + -0.228683 *lens_ipow(lambda, 3) + -5.25735e-06 *lens_ipow(y, 2)*lens_ipow(dx, 2) + -1.18833e-06 *lens_ipow(y, 3)*dy + -0.00215847 *x*dx*lens_ipow(dy, 2) + -3.39926e-05 *x*y*dx*dy + -1.93753e-06 *x*lens_ipow(y, 2)*dx + -2.70406e-05 *lens_ipow(x, 2)*lens_ipow(dy, 2) + 2.15712e-05 *lens_ipow(x, 2)*lens_ipow(dx, 2) + -1.87889e-06 *lens_ipow(x, 2)*y*dy + -6.50916e-08 *lens_ipow(x, 2)*lens_ipow(y, 2) + -1.63496 *lens_ipow(dy, 6) + -4.79483 *lens_ipow(dx, 2)*lens_ipow(dy, 4) + -5.09638 *lens_ipow(dx, 4)*lens_ipow(dy, 2) + -1.2765 *lens_ipow(dx, 6) + -3.46002e-08 *lens_ipow(y, 4)*lens_ipow(dy, 2) + -7.45789e-11 *lens_ipow(y, 6) + 1.26184e-08 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(dy, 2) + -3.54509e-09 *lens_ipow(x, 5)*dx + -7.39851e-11 *lens_ipow(x, 6) + -0.00192761 *lens_ipow(y, 2)*lens_ipow(dx, 6) + 0.00726618 *x*y*dx*lens_ipow(dy, 5) + 0.0146502 *x*y*lens_ipow(dx, 3)*lens_ipow(dy, 3) + 0.00888838 *x*y*lens_ipow(dx, 5)*dy + 8.19731e-08 *x*lens_ipow(y, 4)*dx*lens_ipow(dy, 2) + 6.08274e-08 *lens_ipow(x, 4)*y*lens_ipow(dx, 2)*dy + -0.000251479 *x*y*dx*dy*lens_ipow(lambda, 5) + -0.000285921 *lens_ipow(y, 3)*lens_ipow(dy, 7) + -1.30705e-10 *lens_ipow(y, 7)*lens_ipow(dx, 2)*dy + 0.414686 *lens_ipow(lambda, 11);
