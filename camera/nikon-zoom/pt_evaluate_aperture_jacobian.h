const float dx00 =  + 0.995072  + 1.89009e-06 *y + 0.652066 *lens_ipow(dy, 2) + 1.85633 *lens_ipow(dx, 2) + -0.000339802 *lens_ipow(y, 2) + -0.000959562 *lens_ipow(x, 2) + 0.0079601 *lens_ipow(lambda, 3) + 3.44221e-05 *y*dx*dy + -0.990878 *lens_ipow(dx, 2)*lens_ipow(dy, 2) + 6.59443e-05 *lens_ipow(y, 2)*lens_ipow(lambda, 3) + -0.000645977 *lens_ipow(x, 3)*dx*lens_ipow(dy, 2)*lambda + 5.46095e-07 *lens_ipow(x, 4)*lens_ipow(lambda, 3) + 8.9391e-08 *lens_ipow(y, 4)*lens_ipow(dy, 2)*lens_ipow(lambda, 2) + -0.00228802 *x*lens_ipow(y, 2)*dx*lens_ipow(dy, 4) + 0.000299782 *x*lens_ipow(y, 2)*lens_ipow(dx, 5) + -3.96224e-07 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(lambda, 4) + 0.00279309 *lens_ipow(x, 3)*lens_ipow(dx, 5) + -6.33779e-05 *lens_ipow(x, 4)*lens_ipow(dy, 4) + -5.26933e-12 *lens_ipow(x, 4)*lens_ipow(y, 4) + -5.52884e-12 *lens_ipow(x, 8) + -1.49913 *lens_ipow(dx, 4)*lens_ipow(lambda, 5) + -0.622995 *lens_ipow(dy, 2)*lens_ipow(lambda, 8) + -0.786273 *lens_ipow(x, 2)*lens_ipow(dx, 4)*lens_ipow(dy, 2)*lens_ipow(lambda, 2) + -2.40621e-05 *lens_ipow(x, 3)*lens_ipow(y, 2)*lens_ipow(dx, 3)*lens_ipow(dy, 2) + 1.02867e-10 *lens_ipow(x, 8)*lens_ipow(dx, 2)+0.0f;
const float dx01 =  + 1.89009e-06 *x + 1.11534 *dx*dy + -0.00319777 *y*dx + -0.000679603 *x*y + 3.44221e-05 *x*dx*dy + 0.000131889 *x*y*lens_ipow(lambda, 3) + -0.0958787 *lens_ipow(y, 2)*lens_ipow(dx, 5)*dy + 2.52426e-05 *lens_ipow(y, 3)*dx*lens_ipow(lambda, 4) + 3.39712e-11 *lens_ipow(y, 7)*dx + 3.57564e-07 *x*lens_ipow(y, 3)*lens_ipow(dy, 2)*lens_ipow(lambda, 2) + -0.00228802 *lens_ipow(x, 2)*y*dx*lens_ipow(dy, 4) + 0.000299782 *lens_ipow(x, 2)*y*lens_ipow(dx, 5) + -2.64149e-07 *lens_ipow(x, 3)*y*lens_ipow(lambda, 4) + -4.21546e-12 *lens_ipow(x, 5)*lens_ipow(y, 3) + 212.242 *lens_ipow(dx, 5)*lens_ipow(dy, 3)*lens_ipow(lambda, 2) + -1.41367 *lens_ipow(y, 2)*lens_ipow(dx, 3)*lens_ipow(dy, 5) + -1.2031e-05 *lens_ipow(x, 4)*y*lens_ipow(dx, 3)*lens_ipow(dy, 2)+0.0f;
const float dx02 =  + 61.7449  + -0.00620699 *dy + 53.4254 *lens_ipow(dy, 2) + 166.667 *lens_ipow(dx, 2) + 1.11534 *y*dy + -0.00159888 *lens_ipow(y, 2) + 3.71266 *x*dx + 3.44221e-05 *x*y*dy + -1.98176 *x*dx*lens_ipow(dy, 2) + 0.707743 *lens_ipow(lambda, 5) + -0.000161494 *lens_ipow(x, 4)*lens_ipow(dy, 2)*lambda + -90221.9 *lens_ipow(dx, 6)*lens_ipow(dy, 2) + -0.159798 *lens_ipow(y, 3)*lens_ipow(dx, 4)*dy + 6.31065e-06 *lens_ipow(y, 4)*lens_ipow(lambda, 4) + 4.2464e-12 *lens_ipow(y, 8) + -0.00114401 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(dy, 4) + 0.000749455 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(dx, 4) + 0.00349136 *lens_ipow(x, 4)*lens_ipow(dx, 4) + -5.99653 *x*lens_ipow(dx, 3)*lens_ipow(lambda, 5) + -327837 *lens_ipow(dx, 10) + 1061.21 *y*lens_ipow(dx, 4)*lens_ipow(dy, 3)*lens_ipow(lambda, 2) + -1.41367 *lens_ipow(y, 3)*lens_ipow(dx, 2)*lens_ipow(dy, 5) + -1.04836 *lens_ipow(x, 3)*lens_ipow(dx, 3)*lens_ipow(dy, 2)*lens_ipow(lambda, 2) + -1.80466e-05 *lens_ipow(x, 4)*lens_ipow(y, 2)*lens_ipow(dx, 2)*lens_ipow(dy, 2) + 2.28594e-11 *lens_ipow(x, 9)*dx+0.0f;
const float dx03 =  + -0.00620699 *dx + 106.851 *dx*dy + 1.11534 *y*dx + 1.30413 *x*dy + 3.44221e-05 *x*y*dx + -1.98176 *x*lens_ipow(dx, 2)*dy + -0.000322988 *lens_ipow(x, 4)*dx*dy*lambda + -25777.7 *lens_ipow(dx, 7)*dy + -0.0319596 *lens_ipow(y, 3)*lens_ipow(dx, 5) + 1.78782e-07 *x*lens_ipow(y, 4)*dy*lens_ipow(lambda, 2) + -0.00457604 *lens_ipow(x, 2)*lens_ipow(y, 2)*dx*lens_ipow(dy, 3) + -5.07023e-05 *lens_ipow(x, 5)*lens_ipow(dy, 3) + 636.725 *y*lens_ipow(dx, 5)*lens_ipow(dy, 2)*lens_ipow(lambda, 2) + -2.35612 *lens_ipow(y, 3)*lens_ipow(dx, 3)*lens_ipow(dy, 4) + -1.24599 *x*dy*lens_ipow(lambda, 8) + -0.524182 *lens_ipow(x, 3)*lens_ipow(dx, 4)*dy*lens_ipow(lambda, 2) + -1.2031e-05 *lens_ipow(x, 4)*lens_ipow(y, 2)*lens_ipow(dx, 3)*dy+0.0f;
const float dx04 =  + 0.0238803 *x*lens_ipow(lambda, 2) + 3.53872 *dx*lens_ipow(lambda, 4) + 0.000197833 *x*lens_ipow(y, 2)*lens_ipow(lambda, 2) + -0.000161494 *lens_ipow(x, 4)*dx*lens_ipow(dy, 2) + 3.27657e-07 *lens_ipow(x, 5)*lens_ipow(lambda, 2) + 2.52426e-05 *lens_ipow(y, 4)*dx*lens_ipow(lambda, 3) + 1.78782e-07 *x*lens_ipow(y, 4)*lens_ipow(dy, 2)*lambda + -5.28298e-07 *lens_ipow(x, 3)*lens_ipow(y, 2)*lens_ipow(lambda, 3) + -7.49567 *x*lens_ipow(dx, 4)*lens_ipow(lambda, 4) + 424.483 *y*lens_ipow(dx, 5)*lens_ipow(dy, 3)*lambda + -4.98396 *x*lens_ipow(dy, 2)*lens_ipow(lambda, 7) + -0.524182 *lens_ipow(x, 3)*lens_ipow(dx, 4)*lens_ipow(dy, 2)*lambda+0.0f;
const float dx10 =  + -1.60575e-06 *y + 1.05814 *dx*dy + 8.86919e-08 *lens_ipow(y, 2) + -0.00360782 *x*dy + -0.000644842 *x*y + 2.0696e-08 *x*lens_ipow(y, 2) + 0.000447567 *lens_ipow(y, 2)*dx*dy + -1.82717e-07 *x*lens_ipow(y, 3) + 0.302942 *x*lens_ipow(dx, 2)*lens_ipow(dy, 3) + -0.00316799 *x*y*lens_ipow(dy, 4) + -7.02994e-10 *lens_ipow(x, 5)*y + -1.96863e-05 *lens_ipow(x, 2)*y*dx*lens_ipow(lambda, 3) + 3.87288e-07 *x*lens_ipow(y, 3)*lens_ipow(lambda, 4) + -2.39433e-12 *lens_ipow(x, 3)*lens_ipow(y, 5) + -3.69649e-07 *lens_ipow(x, 4)*y*lens_ipow(dx, 3) + 9.10723e-12 *lens_ipow(x, 7)*dy + 21.6478 *y*lens_ipow(dx, 3)*lens_ipow(dy, 4)*lens_ipow(lambda, 2) + -1.43734e-05 *x*lens_ipow(y, 4)*lens_ipow(dx, 2)*lens_ipow(dy, 3)+0.0f;
const float dx11 =  + 0.994438  + -1.60575e-06 *x + 1.82013 *lens_ipow(dy, 2) + 0.636496 *lens_ipow(dx, 2) + -0.000952833 *lens_ipow(y, 2) + 1.77384e-07 *x*y + -0.000322421 *lens_ipow(x, 2) + 0.0125377 *lens_ipow(lambda, 3) + 2.0696e-08 *lens_ipow(x, 2)*y + 0.000895134 *x*y*dx*dy + -2.74076e-07 *lens_ipow(x, 2)*lens_ipow(y, 2) + -1.77789e-06 *lens_ipow(y, 4)*lens_ipow(dx, 2) + -0.00158399 *lens_ipow(x, 2)*lens_ipow(dy, 4) + -1.17166e-10 *lens_ipow(x, 6) + 4.05113e-07 *lens_ipow(y, 4)*lens_ipow(lambda, 3) + -6.56211e-06 *lens_ipow(x, 3)*dx*lens_ipow(lambda, 3) + -18.8099 *lens_ipow(dx, 8) + 6.76467e-09 *lens_ipow(y, 6)*lens_ipow(dy, 2) + -3.06314e-12 *lens_ipow(y, 8) + 5.80932e-07 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(lambda, 4) + -2.99291e-12 *lens_ipow(x, 4)*lens_ipow(y, 4) + -7.39297e-08 *lens_ipow(x, 5)*lens_ipow(dx, 3) + -1.95615 *lens_ipow(dy, 4)*lens_ipow(lambda, 5) + -0.0346746 *lens_ipow(lambda, 10) + -0.394193 *lens_ipow(dx, 2)*lens_ipow(lambda, 8) + -1.84907 *lens_ipow(y, 2)*lens_ipow(dx, 2)*lens_ipow(dy, 6) + 21.6478 *x*lens_ipow(dx, 3)*lens_ipow(dy, 4)*lens_ipow(lambda, 2) + -2.87468e-05 *lens_ipow(x, 2)*lens_ipow(y, 3)*lens_ipow(dx, 2)*lens_ipow(dy, 3)+0.0f;
const float dx12 =  + 103.681 *dx*dy + 1.27299 *y*dx + 1.05814 *x*dy + 0.291685 *lens_ipow(dx, 2)*dy + 0.000447567 *x*lens_ipow(y, 2)*dy + -7.11156e-07 *lens_ipow(y, 5)*dx + 0.302942 *lens_ipow(x, 2)*dx*lens_ipow(dy, 3) + -6.56211e-06 *lens_ipow(x, 3)*y*lens_ipow(lambda, 3) + -150.479 *y*lens_ipow(dx, 7) + -2.21789e-07 *lens_ipow(x, 5)*y*lens_ipow(dx, 2) + -1.16137e+06 *lens_ipow(dx, 3)*lens_ipow(dy, 7) + -0.788386 *y*dx*lens_ipow(lambda, 8) + -1.23271 *lens_ipow(y, 3)*dx*lens_ipow(dy, 6) + 64.9434 *x*y*lens_ipow(dx, 2)*lens_ipow(dy, 4)*lens_ipow(lambda, 2) + -1.43734e-05 *lens_ipow(x, 2)*lens_ipow(y, 4)*dx*lens_ipow(dy, 3)+0.0f;
const float dx13 =  + 61.7517  + 162.884 *lens_ipow(dy, 2) + 51.8407 *lens_ipow(dx, 2) + 3.64026 *y*dy + 1.05814 *x*dx + -0.00180391 *lens_ipow(x, 2) + 0.0972283 *lens_ipow(dx, 3) + 0.000447567 *x*lens_ipow(y, 2)*dx + 0.698269 *lens_ipow(lambda, 5) + 0.454413 *lens_ipow(x, 2)*lens_ipow(dx, 2)*lens_ipow(dy, 2) + -0.00633598 *lens_ipow(x, 2)*y*lens_ipow(dy, 3) + 1.93276e-09 *lens_ipow(y, 7)*dy + 1.1384e-12 *lens_ipow(x, 8) + -7.82459 *y*lens_ipow(dy, 3)*lens_ipow(lambda, 5) + 6574.73 *lens_ipow(dy, 10) + -2.03239e+06 *lens_ipow(dx, 4)*lens_ipow(dy, 6) + -3.69814 *lens_ipow(y, 3)*lens_ipow(dx, 2)*lens_ipow(dy, 5) + 86.5912 *x*y*lens_ipow(dx, 3)*lens_ipow(dy, 3)*lens_ipow(lambda, 2) + -2.15601e-05 *lens_ipow(x, 2)*lens_ipow(y, 4)*lens_ipow(dx, 2)*lens_ipow(dy, 2)+0.0f;
const float dx14 =  + 0.0376132 *y*lens_ipow(lambda, 2) + 3.49134 *dy*lens_ipow(lambda, 4) + 2.43068e-07 *lens_ipow(y, 5)*lens_ipow(lambda, 2) + -1.96863e-05 *lens_ipow(x, 3)*y*dx*lens_ipow(lambda, 2) + 7.74576e-07 *lens_ipow(x, 2)*lens_ipow(y, 3)*lens_ipow(lambda, 3) + -9.78073 *y*lens_ipow(dy, 4)*lens_ipow(lambda, 4) + -0.346746 *y*lens_ipow(lambda, 9) + -3.15354 *y*lens_ipow(dx, 2)*lens_ipow(lambda, 7) + 43.2956 *x*y*lens_ipow(dx, 3)*lens_ipow(dy, 4)*lambda+0.0f;
const float dx20 =  + 0.0022763  + 0.107079 *lens_ipow(dy, 2) + 0.318294 *lens_ipow(dx, 2) + 0.00318295 *y*dy + 1.74886e-05 *lens_ipow(y, 2) + 0.00947827 *x*dx + 5.15404e-05 *lens_ipow(x, 2) + 0.000463255 *lens_ipow(lambda, 3) + 3.60722e-10 *lens_ipow(y, 3) + 7.40822e-05 *x*y*dx*dy + 1.03738e-05 *lens_ipow(x, 2)*lens_ipow(dy, 2) + 0.012579 *y*lens_ipow(dx, 4)*dy + 0.000193739 *lens_ipow(y, 2)*lens_ipow(dy, 4) + 7.91896e-05 *lens_ipow(x, 2)*lens_ipow(dy, 4) + -4.15932e-11 *lens_ipow(x, 2)*lens_ipow(y, 4) + -1.75815e-10 *lens_ipow(x, 4)*lens_ipow(y, 2) + -6.23185e-09 *lens_ipow(x, 5)*dx + 1.40472e-05 *lens_ipow(y, 2)*lens_ipow(dy, 2)*lens_ipow(lambda, 4) + -1.64294e-08 *lens_ipow(y, 5)*lens_ipow(dx, 2)*dy + -2.49063e-14 *lens_ipow(y, 8) + 0.0907677 *x*lens_ipow(dx, 7) + -1.22662e-08 *lens_ipow(x, 3)*lens_ipow(y, 2)*dx*lens_ipow(lambda, 3) + -0.00231114 *lens_ipow(lambda, 10) + -31.5911 *lens_ipow(dy, 10) + -0.390253 *lens_ipow(dx, 2)*lens_ipow(dy, 2)*lens_ipow(lambda, 6) + 6.35522e-12 *lens_ipow(x, 6)*lens_ipow(y, 2)*lens_ipow(dx, 2) + -1.24009e-15 *lens_ipow(x, 10)+0.0f;
const float dx21 =  + -3.63781e-07  + 0.208823 *dx*dy + 0.0030898 *y*dx + 0.00318295 *x*dy + 3.49772e-05 *x*y + 1.08217e-09 *x*lens_ipow(y, 2) + -8.57862e-08 *lens_ipow(y, 3)*dx + 3.70411e-05 *lens_ipow(x, 2)*dx*dy + 0.012579 *x*lens_ipow(dx, 4)*dy + 0.000387477 *x*y*lens_ipow(dy, 4) + -5.54576e-11 *lens_ipow(x, 3)*lens_ipow(y, 3) + -7.03261e-11 *lens_ipow(x, 5)*y + 2.80944e-05 *x*y*lens_ipow(dy, 2)*lens_ipow(lambda, 4) + -8.21471e-08 *x*lens_ipow(y, 4)*lens_ipow(dx, 2)*dy + -1.9925e-13 *x*lens_ipow(y, 7) + -6.1331e-09 *lens_ipow(x, 4)*y*dx*lens_ipow(lambda, 3) + 0.000364403 *lens_ipow(y, 3)*lens_ipow(dx, 7) + 1.81578e-12 *lens_ipow(x, 7)*y*lens_ipow(dx, 2)+0.0f;
const float dx22 =  + 1.14424  + 6.36787 *lens_ipow(dy, 2) + 19.2649 *lens_ipow(dx, 2) + 0.208823 *y*dy + 0.0015449 *lens_ipow(y, 2) + 0.636588 *x*dx + 0.00473914 *lens_ipow(x, 2) + -2.14465e-08 *lens_ipow(y, 4) + 3.70411e-05 *lens_ipow(x, 2)*y*dy + 12.41 *lens_ipow(dy, 6) + 219.144 *lens_ipow(dx, 2)*lens_ipow(dy, 4) + 69.4182 *lens_ipow(dx, 4)*lens_ipow(dy, 2) + 0.050316 *x*y*lens_ipow(dx, 3)*dy + -1.03864e-09 *lens_ipow(x, 6) + -3.28588e-08 *x*lens_ipow(y, 5)*dx*dy + 0.317687 *lens_ipow(x, 2)*lens_ipow(dx, 6) + -3.06655e-09 *lens_ipow(x, 4)*lens_ipow(y, 2)*lens_ipow(lambda, 3) + 16.1368 *lens_ipow(dx, 4)*lens_ipow(lambda, 6) + 0.000637705 *lens_ipow(y, 4)*lens_ipow(dx, 6) + -0.780506 *x*dx*lens_ipow(dy, 2)*lens_ipow(lambda, 6) + 1.81578e-12 *lens_ipow(x, 7)*lens_ipow(y, 2)*dx+0.0f;
const float dx23 =  + 12.7357 *dx*dy + 0.208823 *y*dx + 0.214159 *x*dy + 0.00318295 *x*y + 3.70411e-05 *lens_ipow(x, 2)*y*dx + 6.91589e-06 *lens_ipow(x, 3)*dy + 74.4601 *dx*lens_ipow(dy, 5) + 292.193 *lens_ipow(dx, 3)*lens_ipow(dy, 3) + 27.7673 *lens_ipow(dx, 5)*dy + 0.012579 *x*y*lens_ipow(dx, 4) + 0.000774955 *x*lens_ipow(y, 2)*lens_ipow(dy, 3) + 0.000105586 *lens_ipow(x, 3)*lens_ipow(dy, 3) + 2.80944e-05 *x*lens_ipow(y, 2)*dy*lens_ipow(lambda, 4) + -1.64294e-08 *x*lens_ipow(y, 5)*lens_ipow(dx, 2) + -315.911 *x*lens_ipow(dy, 9) + -0.780506 *x*lens_ipow(dx, 2)*dy*lens_ipow(lambda, 6)+0.0f;
const float dx24 =  + 0.00138977 *x*lens_ipow(lambda, 2) + 5.61889e-05 *x*lens_ipow(y, 2)*lens_ipow(dy, 2)*lens_ipow(lambda, 3) + -9.19965e-09 *lens_ipow(x, 4)*lens_ipow(y, 2)*dx*lens_ipow(lambda, 2) + 19.3642 *lens_ipow(dx, 5)*lens_ipow(lambda, 5) + -0.0231114 *x*lens_ipow(lambda, 9) + -2.34152 *x*lens_ipow(dx, 2)*lens_ipow(dy, 2)*lens_ipow(lambda, 5)+0.0f;
const float dx30 =  + 8.25174e-08  + -8.67855e-08 *y + 0.207388 *dx*dy + 0.00313476 *y*dx + 0.00303303 *x*dy + 3.32729e-05 *x*y + 2.7554e-05 *lens_ipow(y, 2)*dx*dy + 0.00524085 *x*lens_ipow(dy, 3) + 0.000115327 *x*y*lens_ipow(dy, 2) + -1.39037e-07 *lens_ipow(x, 3)*dy + 0.0181226 *y*dx*lens_ipow(dy, 4) + 0.000193938 *x*y*lens_ipow(dx, 4) + -7.76599e-11 *x*lens_ipow(y, 5) + -0.185364 *x*lens_ipow(dy, 7) + 1.99567e-08 *x*lens_ipow(y, 3)*lens_ipow(lambda, 4) + 1.56513e-09 *lens_ipow(x, 2)*lens_ipow(y, 4)*dx*dy + -3.07516e-09 *lens_ipow(x, 5)*lens_ipow(dy, 3) + -5.35475e-13 *lens_ipow(x, 5)*lens_ipow(y, 3) + 5.1205e-09 *lens_ipow(x, 4)*y*lens_ipow(dx, 3)*lens_ipow(dy, 2)+0.0f;
const float dx31 =  + 0.0023323  + -8.67855e-08 *x + 0.314181 *lens_ipow(dy, 2) + 0.106684 *lens_ipow(dx, 2) + 0.00941157 *y*dy + 5.0703e-05 *lens_ipow(y, 2) + 0.00313476 *x*dx + 1.66364e-05 *lens_ipow(x, 2) + 0.000259097 *lens_ipow(lambda, 3) + 0.00282923 *y*lens_ipow(dx, 2)*dy + 5.51079e-05 *x*y*dx*dy + 5.76635e-05 *lens_ipow(x, 2)*lens_ipow(dy, 2) + 0.0155056 *y*lens_ipow(dy, 5) + 0.000525634 *lens_ipow(y, 2)*lens_ipow(dx, 4) + -1.03913e-08 *lens_ipow(y, 5)*dy + 0.0181226 *x*dx*lens_ipow(dy, 4) + 9.69688e-05 *lens_ipow(x, 2)*lens_ipow(dx, 4) + -1.9415e-10 *lens_ipow(x, 2)*lens_ipow(y, 4) + 7.58444e-10 *lens_ipow(y, 6)*lens_ipow(dx, 2) + 2.99351e-08 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(lambda, 4) + 2.08684e-09 *lens_ipow(x, 3)*lens_ipow(y, 3)*dx*dy + -2.67737e-13 *lens_ipow(x, 6)*lens_ipow(y, 2) + -0.0290378 *lens_ipow(y, 2)*lens_ipow(dx, 2)*lens_ipow(dy, 6) + -0.271565 *lens_ipow(y, 2)*lens_ipow(dx, 6)*lens_ipow(dy, 2) + -1.73523e-15 *lens_ipow(y, 10) + 1.0241e-09 *lens_ipow(x, 5)*lens_ipow(dx, 3)*lens_ipow(dy, 2)+0.0f;
const float dx32 =  + 12.6122 *dx*dy + 0.213368 *y*dx + 0.207388 *x*dy + 0.00313476 *x*y + 0.00282923 *lens_ipow(y, 2)*dx*dy + 2.7554e-05 *x*lens_ipow(y, 2)*dy + -38.4536 *dx*lens_ipow(dy, 5) + 140.219 *lens_ipow(dx, 5)*dy + 0.000700845 *lens_ipow(y, 3)*lens_ipow(dx, 3) + 0.0181226 *x*y*lens_ipow(dy, 4) + 0.000387875 *lens_ipow(x, 2)*y*lens_ipow(dx, 3) + 2.16698e-10 *lens_ipow(y, 7)*dx + 5.2171e-10 *lens_ipow(x, 3)*lens_ipow(y, 4)*dy + -0.0193585 *lens_ipow(y, 3)*dx*lens_ipow(dy, 6) + -0.54313 *lens_ipow(y, 3)*lens_ipow(dx, 5)*lens_ipow(dy, 2) + 3.0723e-09 *lens_ipow(x, 5)*y*lens_ipow(dx, 2)*lens_ipow(dy, 2)+0.0f;
const float dx33 =  + 1.146  + 18.7821 *lens_ipow(dy, 2) + 6.30612 *lens_ipow(dx, 2) + 0.628362 *y*dy + 0.00470579 *lens_ipow(y, 2) + 0.207388 *x*dx + 0.00151652 *lens_ipow(x, 2) + 0.00141462 *lens_ipow(y, 2)*lens_ipow(dx, 2) + 2.7554e-05 *x*lens_ipow(y, 2)*dx + 0.00786128 *lens_ipow(x, 2)*lens_ipow(dy, 2) + 0.000115327 *lens_ipow(x, 2)*y*dy + -3.47593e-08 *lens_ipow(x, 4) + -96.134 *lens_ipow(dx, 2)*lens_ipow(dy, 4) + 23.3698 *lens_ipow(dx, 6) + 0.0387639 *lens_ipow(y, 2)*lens_ipow(dy, 4) + -1.73188e-09 *lens_ipow(y, 6) + 0.0724905 *x*y*dx*lens_ipow(dy, 3) + -0.648775 *lens_ipow(x, 2)*lens_ipow(dy, 6) + 5.2171e-10 *lens_ipow(x, 3)*lens_ipow(y, 4)*dx + -1.53758e-09 *lens_ipow(x, 6)*lens_ipow(dy, 2) + 2.28935 *lens_ipow(dy, 2)*lens_ipow(lambda, 8) + -0.0580756 *lens_ipow(y, 3)*lens_ipow(dx, 2)*lens_ipow(dy, 5) + -0.181043 *lens_ipow(y, 3)*lens_ipow(dx, 6)*dy + 2.0482e-09 *lens_ipow(x, 5)*y*lens_ipow(dx, 3)*dy+0.0f;
const float dx34 =  + 0.000777292 *y*lens_ipow(lambda, 2) + 3.99135e-08 *lens_ipow(x, 2)*lens_ipow(y, 3)*lens_ipow(lambda, 3) + 6.10494 *lens_ipow(dy, 3)*lens_ipow(lambda, 7)+0.0f;
const float dx40 =  + -5.86881e-07  + -0.00280849 *dx + 1.27446e-07 *y + -4.68889e-05 *x + -1.83583e-05 *dx*dy + -5.44112e-05 *y*dx*dy + 2.92596e-06 *y*lens_ipow(dx, 2) + -7.26883e-05 *x*lens_ipow(dx, 2) + -1.99624e-08 *x*lens_ipow(y, 2) + 1.29851e-07 *lens_ipow(x, 2)*dy + 4.36522e-09 *lens_ipow(x, 4)*dx + -1.193e-06 *x*y*dy*lens_ipow(lambda, 3) + -0.0535531 *y*lens_ipow(dx, 3)*lens_ipow(dy, 3) + -0.00182858 *x*lens_ipow(dy, 6) + -4.61249e-09 *lens_ipow(x, 5)*lens_ipow(dy, 2) + 0.00340751 *x*y*lens_ipow(dx, 6)*dy + -4.07246e-12 *x*lens_ipow(y, 6)*lens_ipow(dx, 2) + 0.000785302 *lens_ipow(x, 3)*lens_ipow(dx, 2)*lens_ipow(dy, 4) + 1.5028e-13 *lens_ipow(x, 5)*lens_ipow(y, 3)*dy+0.0f;
const float dx41 =  + -1.37711e-07  + -0.00273422 *dy + 1.6236e-05 *dx + -4.57379e-05 *y + 1.27446e-07 *x + 2.47658e-05 *lens_ipow(dx, 2) + -5.44112e-05 *x*dx*dy + 2.92596e-06 *x*lens_ipow(dx, 2) + -1.99624e-08 *lens_ipow(x, 2)*y + 2.44486e-05 *lens_ipow(y, 2)*lens_ipow(dy, 3) + -5.54469e-07 *lens_ipow(y, 3)*lens_ipow(dx, 2) + -8.36886e-11 *lens_ipow(y, 5) + -5.96499e-07 *lens_ipow(x, 2)*dy*lens_ipow(lambda, 3) + -0.0535531 *x*lens_ipow(dx, 3)*lens_ipow(dy, 3) + -0.0203234 *lens_ipow(y, 2)*lens_ipow(dx, 6)*dy + -3.45628e-06 *lens_ipow(y, 4)*lens_ipow(dx, 2)*lens_ipow(dy, 3) + 0.00170375 *lens_ipow(x, 2)*lens_ipow(dx, 6)*dy + -1.22174e-11 *lens_ipow(x, 2)*lens_ipow(y, 5)*lens_ipow(dx, 2) + 7.514e-14 *lens_ipow(x, 6)*lens_ipow(y, 2)*dy+0.0f;
const float dx42 =  + -3.13117e-05  + 0.000600637 *dy + -0.169893 *dx + 1.6236e-05 *y + -0.00280849 *x + 4.95316e-05 *y*dx + -1.83583e-05 *x*dy + 0.398808 *dx*lens_ipow(dy, 2) + -5.44112e-05 *x*y*dy + 5.85192e-06 *x*y*dx + -7.26883e-05 *lens_ipow(x, 2)*dx + -2.77235e-07 *lens_ipow(y, 4)*dx + 8.73043e-10 *lens_ipow(x, 5) + 58.5532 *lens_ipow(dx, 7) + -0.160659 *x*y*lens_ipow(dx, 2)*lens_ipow(dy, 3) + -0.0406468 *lens_ipow(y, 3)*lens_ipow(dx, 5)*dy + -1.38251e-06 *lens_ipow(y, 5)*dx*lens_ipow(dy, 3) + 0.0102225 *lens_ipow(x, 2)*y*lens_ipow(dx, 5)*dy + -4.07246e-12 *lens_ipow(x, 2)*lens_ipow(y, 6)*dx + 0.000392651 *lens_ipow(x, 4)*dx*lens_ipow(dy, 4)+0.0f;
const float dx43 =  + -0.169913 *dy + 0.000600637 *dx + -0.00273422 *y + -1.83583e-05 *x*dx + 0.398808 *lens_ipow(dx, 2)*dy + -5.44112e-05 *x*y*dx + 4.32838e-08 *lens_ipow(x, 3) + 2.44486e-05 *lens_ipow(y, 3)*lens_ipow(dy, 2) + -5.96499e-07 *lens_ipow(x, 2)*y*lens_ipow(lambda, 3) + 76.2705 *lens_ipow(dy, 7) + -0.160659 *x*y*lens_ipow(dx, 3)*lens_ipow(dy, 2) + -0.00548573 *lens_ipow(x, 2)*lens_ipow(dy, 5) + -1.5375e-09 *lens_ipow(x, 6)*dy + -0.00677447 *lens_ipow(y, 3)*lens_ipow(dx, 6) + -2.07377e-06 *lens_ipow(y, 5)*lens_ipow(dx, 2)*lens_ipow(dy, 2) + 0.00170375 *lens_ipow(x, 2)*y*lens_ipow(dx, 6) + 0.000785302 *lens_ipow(x, 4)*lens_ipow(dx, 2)*lens_ipow(dy, 3) + 2.50467e-14 *lens_ipow(x, 6)*lens_ipow(y, 3)+0.0f;
const float dx44 =  + 0.15019  + -0.377619 *lens_ipow(lambda, 2) + -1.7895e-06 *lens_ipow(x, 2)*y*dy*lens_ipow(lambda, 2) + 2.50774 *lens_ipow(lambda, 10)+0.0f;
