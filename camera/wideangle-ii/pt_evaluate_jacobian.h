const float dx00 =  + -1.79206  + 0.702177 *lens_ipow(dy, 2) + 1.96079 *lens_ipow(dx, 2) + 0.0634463 *y*dy + 0.00273604 *lens_ipow(y, 2) + 0.228484 *x*dx + 0.00823362 *lens_ipow(x, 2) + -0.126033 *lens_ipow(lambda, 3) + 0.00344132 *x*dx*lens_ipow(lambda, 2) + -0.0161354 *x*y*dx*dy + 0.00941357 *lens_ipow(x, 2)*lens_ipow(dy, 2) + -15.3263 *lens_ipow(dx, 2)*lens_ipow(dy, 2)*lens_ipow(lambda, 2) + -0.00904287 *lens_ipow(y, 3)*lens_ipow(dx, 2)*dy + 1.49982 *x*lens_ipow(dx, 5) + 0.000188357 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(dx, 2) + -0.0273947 *y*dy*lens_ipow(lambda, 6) + -4.36807e-11 *lens_ipow(y, 8) + -0.0618198 *lens_ipow(x, 2)*y*lens_ipow(dy, 5) + -1.42539e-09 *lens_ipow(x, 4)*lens_ipow(y, 4) + -6.48277e-05 *lens_ipow(x, 5)*lens_ipow(dx, 3) + -6.45333e-10 *lens_ipow(x, 6)*lens_ipow(y, 2) + -2.65837e-10 *lens_ipow(x, 8) + 0.526558 *lens_ipow(lambda, 10) + 8.03647 *lens_ipow(dx, 2)*lens_ipow(lambda, 8) + 0.000147137 *lens_ipow(x, 5)*y*lens_ipow(dx, 3)*dy+0.0f;
const float dx01 =  + 0.000508213 *dx + 1.40919 *dx*dy + 0.102822 *y*dx + 0.0634463 *x*dy + 0.00547208 *x*y + -0.0747238 *lens_ipow(dx, 2)*dy + -0.00806772 *lens_ipow(x, 2)*dx*dy + 8.14437 *y*lens_ipow(dx, 3)*lens_ipow(dy, 2) + -0.0271286 *x*lens_ipow(y, 2)*lens_ipow(dx, 2)*dy + 0.000125571 *lens_ipow(x, 3)*y*lens_ipow(dx, 2) + -12.738 *dx*lens_ipow(dy, 3)*lens_ipow(lambda, 3) + 0.183187 *lens_ipow(y, 3)*lens_ipow(dx, 5) + 4.70736e-07 *lens_ipow(y, 6)*dx*dy + -0.0273947 *x*dy*lens_ipow(lambda, 6) + -3.49446e-10 *x*lens_ipow(y, 7) + -0.0206066 *lens_ipow(x, 3)*lens_ipow(dy, 5) + -1.14031e-09 *lens_ipow(x, 5)*lens_ipow(y, 3) + -1.84381e-10 *lens_ipow(x, 7)*y + -15.785 *y*lens_ipow(dx, 5)*lens_ipow(lambda, 4) + 2.45229e-05 *lens_ipow(x, 6)*lens_ipow(dx, 3)*dy+0.0f;
const float dx02 =  + 11.3942  + 0.000508213 *y + 18.3599 *lens_ipow(dy, 2) + 1.40919 *y*dy + 0.051411 *lens_ipow(y, 2) + 3.92159 *x*dx + 0.114242 *lens_ipow(x, 2) + -0.149448 *y*dx*dy + 180.882 *lens_ipow(dx, 2)*lens_ipow(lambda, 2) + 0.00172066 *lens_ipow(x, 2)*lens_ipow(lambda, 2) + -0.00806772 *lens_ipow(x, 2)*y*dy + 12.2166 *lens_ipow(y, 2)*lens_ipow(dx, 2)*lens_ipow(dy, 2) + -30.6525 *x*dx*lens_ipow(dy, 2)*lens_ipow(lambda, 2) + -0.0180857 *x*lens_ipow(y, 3)*dx*dy + 3.74955 *lens_ipow(x, 2)*lens_ipow(dx, 4) + 0.000125571 *lens_ipow(x, 3)*lens_ipow(y, 2)*dx + -12.738 *y*lens_ipow(dy, 3)*lens_ipow(lambda, 3) + 0.228983 *lens_ipow(y, 4)*lens_ipow(dx, 4) + 6.7248e-08 *lens_ipow(y, 7)*dy + -3.24138e-05 *lens_ipow(x, 6)*lens_ipow(dx, 2) + -1.17839e+07 *lens_ipow(dx, 4)*lens_ipow(dy, 6) + -39.4624 *lens_ipow(y, 2)*lens_ipow(dx, 4)*lens_ipow(lambda, 4) + 16.0729 *x*dx*lens_ipow(lambda, 8) + 7.35687e-05 *lens_ipow(x, 6)*y*lens_ipow(dx, 2)*dy+0.0f;
const float dx03 =  + -0.052806 *dy + 36.7197 *dx*dy + 1.40919 *y*dx + 1.40435 *x*dy + 0.0634463 *x*y + -0.0747238 *y*lens_ipow(dx, 2) + -0.00806772 *lens_ipow(x, 2)*y*dx + 0.00627571 *lens_ipow(x, 3)*dy + 8.14437 *lens_ipow(y, 2)*lens_ipow(dx, 3)*dy + -30.6525 *x*lens_ipow(dx, 2)*dy*lens_ipow(lambda, 2) + -0.00904287 *x*lens_ipow(y, 3)*lens_ipow(dx, 2) + -38.2141 *y*dx*lens_ipow(dy, 2)*lens_ipow(lambda, 3) + 6.7248e-08 *lens_ipow(y, 7)*dx + -0.0273947 *x*y*lens_ipow(lambda, 6) + -0.103033 *lens_ipow(x, 3)*y*lens_ipow(dy, 4) + -1.41407e+07 *lens_ipow(dx, 5)*lens_ipow(dy, 5) + 2.45229e-05 *lens_ipow(x, 6)*y*lens_ipow(dx, 3)+0.0f;
const float dx04 =  + -0.378099 *x*lens_ipow(lambda, 2) + 120.588 *lens_ipow(dx, 3)*lambda + 0.00344132 *lens_ipow(x, 2)*dx*lambda + -30.6525 *x*lens_ipow(dx, 2)*lens_ipow(dy, 2)*lambda + -38.2141 *y*dx*lens_ipow(dy, 3)*lens_ipow(lambda, 2) + -0.164368 *x*y*dy*lens_ipow(lambda, 5) + -31.5699 *lens_ipow(y, 2)*lens_ipow(dx, 5)*lens_ipow(lambda, 3) + 5.26558 *x*lens_ipow(lambda, 9) + 64.2917 *x*lens_ipow(dx, 2)*lens_ipow(lambda, 7)+0.0f;
const float dx10 =  + -0.000642049 *dy + 1.42244 *dx*dy + 0.070137 *y*dx + 0.10091 *x*dy + 0.00567694 *x*y + -0.012487 *y*dx*lens_ipow(lambda, 2) + -0.349981 *y*dx*lens_ipow(dy, 2) + -0.0642196 *y*lens_ipow(dx, 3) + -0.00822533 *lens_ipow(y, 2)*dx*dy + 0.379005 *x*lens_ipow(dy, 3) + -10.9848 *lens_ipow(dx, 3)*dy*lens_ipow(lambda, 2) + 0.000231732 *x*lens_ipow(y, 3)*lens_ipow(dy, 2) + -5.92849e-08 *x*lens_ipow(y, 5) + -9.56589e-08 *lens_ipow(x, 3)*lens_ipow(y, 3) + -4.20249e-08 *lens_ipow(x, 2)*lens_ipow(y, 5)*dx + -1.28995e-05 *lens_ipow(x, 3)*y*lens_ipow(lambda, 4) + -0.000193596 *lens_ipow(x, 4)*y*dx*lens_ipow(dy, 2) + -2.5727e-10 *lens_ipow(x, 7)*y + -2.17193e-07 *lens_ipow(y, 7)*dx*lens_ipow(dy, 2) + -0.122343 *x*dy*lens_ipow(lambda, 8)+0.0f;
const float dx11 =  + -1.78474  + -0.103219 *lens_ipow(lambda, 2) + 2.16252 *lens_ipow(dy, 2) + 0.687603 *lens_ipow(dx, 2) + 0.24108 *y*dy + 0.000160172 *y*dx + 0.00833531 *lens_ipow(y, 2) + 0.070137 *x*dx + 0.00283847 *lens_ipow(x, 2) + -0.0217767 *y*dy*lens_ipow(lambda, 2) + 0.0141417 *lens_ipow(y, 2)*lens_ipow(dx, 2) + -0.012487 *x*dx*lens_ipow(lambda, 2) + -0.349981 *x*dx*lens_ipow(dy, 2) + -0.0642196 *x*lens_ipow(dx, 3) + -0.0164507 *x*y*dx*dy + -19.79 *lens_ipow(dx, 2)*lens_ipow(dy, 2)*lens_ipow(lambda, 2) + 0.000347598 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(dy, 2) + -1.48212e-07 *lens_ipow(x, 2)*lens_ipow(y, 4) + -7.17442e-08 *lens_ipow(x, 4)*lens_ipow(y, 2) + -5.14162e-05 *lens_ipow(y, 5)*lens_ipow(dy, 3) + -8.45686e-07 *lens_ipow(y, 6)*lens_ipow(dx, 2) + -2.48962e-10 *lens_ipow(y, 8) + -7.00415e-08 *lens_ipow(x, 3)*lens_ipow(y, 4)*dx + -3.22486e-06 *lens_ipow(x, 4)*lens_ipow(lambda, 4) + -3.87192e-05 *lens_ipow(x, 5)*dx*lens_ipow(dy, 2) + -3.21587e-11 *lens_ipow(x, 8) + -101.879 *lens_ipow(dy, 6)*lens_ipow(lambda, 3) + 0.582681 *lens_ipow(lambda, 10) + -1.52035e-06 *x*lens_ipow(y, 6)*dx*lens_ipow(dy, 2)+0.0f;
const float dx12 =  + 36.6005 *dx*dy + 1.37521 *y*dx + 8.00862e-05 *lens_ipow(y, 2) + 1.42244 *x*dy + 0.070137 *x*y + 0.00942781 *lens_ipow(y, 3)*dx + -0.012487 *x*y*lens_ipow(lambda, 2) + -0.349981 *x*y*lens_ipow(dy, 2) + -0.192659 *x*y*lens_ipow(dx, 2) + -0.00822533 *x*lens_ipow(y, 2)*dy + -39.5799 *y*dx*lens_ipow(dy, 2)*lens_ipow(lambda, 2) + -32.9545 *x*lens_ipow(dx, 2)*dy*lens_ipow(lambda, 2) + -2.41625e-07 *lens_ipow(y, 7)*dx + -1.40083e-08 *lens_ipow(x, 3)*lens_ipow(y, 5) + -3.87192e-05 *lens_ipow(x, 5)*y*lens_ipow(dy, 2) + -2.17193e-07 *x*lens_ipow(y, 7)*lens_ipow(dy, 2)+0.0f;
const float dx13 =  + 11.0908  + -0.000642049 *x + 56.2843 *lens_ipow(dy, 2) + 18.3002 *lens_ipow(dx, 2) + 4.32505 *y*dy + 0.12054 *lens_ipow(y, 2) + 1.42244 *x*dx + 0.0504552 *lens_ipow(x, 2) + 1.57168 *lens_ipow(lambda, 3) + -0.0108884 *lens_ipow(y, 2)*lens_ipow(lambda, 2) + -0.699962 *x*y*dx*dy + -0.00822533 *x*lens_ipow(y, 2)*dx + 0.568507 *lens_ipow(x, 2)*lens_ipow(dy, 2) + -39.5799 *y*lens_ipow(dx, 2)*dy*lens_ipow(lambda, 2) + -10.9848 *x*lens_ipow(dx, 3)*lens_ipow(lambda, 2) + 0.000231732 *lens_ipow(x, 2)*lens_ipow(y, 3)*dy + -2.57081e-05 *lens_ipow(y, 6)*lens_ipow(dy, 2) + -7.74384e-05 *lens_ipow(x, 5)*y*dx*dy + -611.276 *y*lens_ipow(dy, 5)*lens_ipow(lambda, 3) + -4.34386e-07 *x*lens_ipow(y, 7)*dx*dy + -0.0611717 *lens_ipow(x, 2)*lens_ipow(lambda, 8)+0.0f;
const float dx14 =  + -0.206439 *y*lambda + 4.71503 *dy*lens_ipow(lambda, 2) + -0.0217767 *lens_ipow(y, 2)*dy*lambda + -0.024974 *x*y*dx*lambda + -39.5799 *y*lens_ipow(dx, 2)*lens_ipow(dy, 2)*lambda + -21.9697 *x*lens_ipow(dx, 3)*dy*lambda + -1.28995e-05 *lens_ipow(x, 4)*y*lens_ipow(lambda, 3) + -305.638 *y*lens_ipow(dy, 6)*lens_ipow(lambda, 2) + 5.82681 *y*lens_ipow(lambda, 9) + -0.489374 *lens_ipow(x, 2)*dy*lens_ipow(lambda, 7)+0.0f;
const float dx20 =  + -0.0269329  + 4.96973e-08 *y + -4.64787e-07 *x + -0.000182299 *y*dy + -3.7398e-06 *lens_ipow(y, 2) + -0.000285772 *x*dx + 5.13997e-05 *lens_ipow(x, 2) + 0.00518719 *lens_ipow(lambda, 3) + -0.0166865 *lens_ipow(dy, 4) + -0.000100659 *x*y*dx*dy + 0.000101545 *lens_ipow(x, 2)*lens_ipow(dy, 2) + -2.77013e-06 *lens_ipow(x, 2)*y*dy + -2.75211e-08 *lens_ipow(x, 3)*y*dy + -9.26389e-05 *lens_ipow(y, 3)*lens_ipow(dx, 2)*dy + 7.68238e-05 *x*lens_ipow(y, 2)*dx*lens_ipow(dy, 2) + -4.05792e-05 *lens_ipow(x, 2)*y*lens_ipow(dy, 3) + 2.20974e-06 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(dx, 2) + -3.10889e-05 *lens_ipow(x, 3)*lens_ipow(dx, 3) + 0.00050185 *lens_ipow(x, 2)*lens_ipow(dx, 2)*lens_ipow(lambda, 3) + -1.62812e-08 *lens_ipow(x, 6)*lens_ipow(dy, 2) + -9.12446e-12 *lens_ipow(x, 6)*lens_ipow(y, 2) + -3.54006e-12 *lens_ipow(x, 8) + -0.593565 *lens_ipow(dx, 2)*lens_ipow(dy, 2)*lens_ipow(lambda, 5) + -3.46583e-07 *lens_ipow(x, 4)*lens_ipow(lambda, 5) + -0.0255423 *lens_ipow(lambda, 10) + -4.11955 *lens_ipow(dx, 6)*lens_ipow(lambda, 4) + 3.20863e-06 *lens_ipow(y, 3)*dy*lens_ipow(lambda, 6) + 1.87392e-05 *lens_ipow(x, 3)*y*dx*dy*lens_ipow(lambda, 4)+0.0f;
const float dx21 =  + 4.96973e-08 *x + 0.00861408 *dx*dy + 0.000399222 *y*dx + -0.000182299 *x*dy + -7.47959e-06 *x*y + -5.03297e-05 *lens_ipow(x, 2)*dx*dy + -9.23377e-07 *lens_ipow(x, 3)*dy + -6.88028e-09 *lens_ipow(x, 4)*dy + 0.0330432 *y*lens_ipow(dx, 5) + -0.000277917 *x*lens_ipow(y, 2)*lens_ipow(dx, 2)*dy + 7.68238e-05 *lens_ipow(x, 2)*y*dx*lens_ipow(dy, 2) + -1.35264e-05 *lens_ipow(x, 3)*lens_ipow(dy, 3) + 1.47316e-06 *lens_ipow(x, 3)*y*lens_ipow(dx, 2) + 2.03052e-05 *lens_ipow(y, 4)*dx*lens_ipow(dy, 3) + -8.14927e-05 *lens_ipow(y, 4)*lens_ipow(dx, 3)*dy + -2.60699e-12 *lens_ipow(x, 7)*y + 9.62589e-06 *x*lens_ipow(y, 2)*dy*lens_ipow(lambda, 6) + 4.6848e-06 *lens_ipow(x, 4)*dx*dy*lens_ipow(lambda, 4)+0.0f;
const float dx22 =  + -0.38828  + 0.443308 *lens_ipow(dy, 2) + 1.21386 *lens_ipow(dx, 2) + 0.00861408 *y*dy + 0.000199611 *lens_ipow(y, 2) + -0.000142886 *lens_ipow(x, 2) + 0.0153536 *lens_ipow(lambda, 4) + -5.03297e-05 *lens_ipow(x, 2)*y*dy + 0.0826081 *lens_ipow(y, 2)*lens_ipow(dx, 4) + -0.000185278 *x*lens_ipow(y, 3)*dx*dy + 3.84119e-05 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(dy, 2) + 1.47316e-06 *lens_ipow(x, 3)*lens_ipow(y, 2)*dx + -2.33167e-05 *lens_ipow(x, 4)*lens_ipow(dx, 2) + 0.000334567 *lens_ipow(x, 3)*dx*lens_ipow(lambda, 3) + -5272.53 *lens_ipow(dx, 4)*lens_ipow(dy, 4) + 4.06104e-06 *lens_ipow(y, 5)*lens_ipow(dy, 3) + -4.88956e-05 *lens_ipow(y, 5)*lens_ipow(dx, 2)*dy + -1.18713 *x*dx*lens_ipow(dy, 2)*lens_ipow(lambda, 5) + -19568.6 *lens_ipow(dx, 10) + -24.7173 *x*lens_ipow(dx, 5)*lens_ipow(lambda, 4) + 4.6848e-06 *lens_ipow(x, 4)*y*dy*lens_ipow(lambda, 4)+0.0f;
const float dx23 =  + 0.886617 *dx*dy + 0.00861408 *y*dx + -0.000182299 *x*y + -0.0667459 *x*lens_ipow(dy, 3) + -5.03297e-05 *lens_ipow(x, 2)*y*dx + 6.76969e-05 *lens_ipow(x, 3)*dy + -9.23377e-07 *lens_ipow(x, 3)*y + -6.88028e-09 *lens_ipow(x, 4)*y + -9.26389e-05 *x*lens_ipow(y, 3)*lens_ipow(dx, 2) + 7.68238e-05 *lens_ipow(x, 2)*lens_ipow(y, 2)*dx*dy + -4.05792e-05 *lens_ipow(x, 3)*y*lens_ipow(dy, 2) + -4218.03 *lens_ipow(dx, 5)*lens_ipow(dy, 3) + 1.21831e-05 *lens_ipow(y, 5)*dx*lens_ipow(dy, 2) + -1.62985e-05 *lens_ipow(y, 5)*lens_ipow(dx, 3) + -4.65176e-09 *lens_ipow(x, 7)*dy + -1.18713 *x*lens_ipow(dx, 2)*dy*lens_ipow(lambda, 5) + 3.20863e-06 *x*lens_ipow(y, 3)*lens_ipow(lambda, 6) + 4.6848e-06 *lens_ipow(x, 4)*y*dx*lens_ipow(lambda, 4)+0.0f;
const float dx24 =  + 0.0155616 *x*lens_ipow(lambda, 2) + 0.0614146 *dx*lens_ipow(lambda, 3) + 0.00050185 *lens_ipow(x, 3)*lens_ipow(dx, 2)*lens_ipow(lambda, 2) + -2.96782 *x*lens_ipow(dx, 2)*lens_ipow(dy, 2)*lens_ipow(lambda, 4) + -3.46583e-07 *lens_ipow(x, 5)*lens_ipow(lambda, 4) + -0.255423 *x*lens_ipow(lambda, 9) + -16.4782 *x*lens_ipow(dx, 6)*lens_ipow(lambda, 3) + 1.92518e-05 *x*lens_ipow(y, 3)*dy*lens_ipow(lambda, 5) + 1.87392e-05 *lens_ipow(x, 4)*y*dx*dy*lens_ipow(lambda, 3)+0.0f;
const float dx30 =  + -1.91681e-07 *y + -3.06053e-05 *y*dx + -3.09542e-08 *lens_ipow(y, 2) + -0.000421085 *x*dy + 7.83319e-05 *x*y + -0.00126511 *y*dx*lens_ipow(dy, 2) + 0.000131163 *lens_ipow(x, 2)*dx*dy + -0.00265214 *lens_ipow(y, 2)*lens_ipow(dx, 3)*dy + -4.80643e-05 *lens_ipow(y, 3)*lens_ipow(dx, 3) + 2.11365e-06 *x*lens_ipow(y, 3)*lens_ipow(dy, 2) + -1.00206e-09 *x*lens_ipow(y, 5) + -0.00016648 *lens_ipow(x, 2)*y*dx*lens_ipow(dy, 2) + -2.51928e-09 *lens_ipow(x, 3)*lens_ipow(y, 3) + -2.30874e-09 *lens_ipow(x, 5)*y + -3.58501e-07 *lens_ipow(y, 5)*dx*lens_ipow(dy, 2) + 1.18033e-06 *lens_ipow(x, 5)*lens_ipow(dy, 3) + 1.22149e-07 *lens_ipow(x, 2)*lens_ipow(y, 4)*dx*dy*lambda + -3.23564e-07 *lens_ipow(x, 3)*y*lens_ipow(lambda, 5) + -0.00118255 *x*dy*lens_ipow(lambda, 8) + 8.52964e-10 *lens_ipow(x, 6)*lens_ipow(y, 2)*dx*dy + 3.89081e-10 *lens_ipow(x, 7)*y*lens_ipow(dx, 2)+0.0f;
const float dx31 =  + -0.0268871  + -1.91681e-07 *x + -0.00767475 *lens_ipow(dx, 2) + -0.000370774 *y*dy + 4.9643e-05 *lens_ipow(y, 2) + -3.06053e-05 *x*dx + -6.19085e-08 *x*y + 3.9166e-05 *lens_ipow(x, 2) + 0.0051492 *lens_ipow(lambda, 3) + -0.00126511 *x*dx*lens_ipow(dy, 2) + -0.00530428 *x*y*lens_ipow(dx, 3)*dy + -0.000144193 *x*lens_ipow(y, 2)*lens_ipow(dx, 3) + 3.17047e-06 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(dy, 2) + -2.50516e-09 *lens_ipow(x, 2)*lens_ipow(y, 4) + -5.54933e-05 *lens_ipow(x, 3)*dx*lens_ipow(dy, 2) + -1.88946e-09 *lens_ipow(x, 4)*lens_ipow(y, 2) + -3.8479e-10 *lens_ipow(x, 6) + 4.30022e-06 *lens_ipow(y, 4)*lens_ipow(dy, 4) + -4.16173e-12 *lens_ipow(y, 8) + -1.79251e-06 *x*lens_ipow(y, 4)*dx*lens_ipow(dy, 2) + -0.154681 *lens_ipow(dy, 4)*lens_ipow(lambda, 5) + 7.05762e-06 *lens_ipow(y, 4)*lens_ipow(dy, 2)*lens_ipow(lambda, 3) + 1.62866e-07 *lens_ipow(x, 3)*lens_ipow(y, 3)*dx*dy*lambda + -8.0891e-08 *lens_ipow(x, 4)*lens_ipow(lambda, 5) + -0.0247506 *lens_ipow(lambda, 10) + -4.70624e-07 *lens_ipow(y, 4)*lens_ipow(lambda, 6) + 2.43704e-10 *lens_ipow(x, 7)*y*dx*dy + 4.86351e-11 *lens_ipow(x, 8)*lens_ipow(dx, 2)+0.0f;
const float dx32 =  + -1.97645e-05  + 0.759587 *dx*dy + -0.0153495 *y*dx + -3.06053e-05 *x*y + -0.00126511 *x*y*lens_ipow(dy, 2) + 4.3721e-05 *lens_ipow(x, 3)*dy + -0.00795641 *x*lens_ipow(y, 2)*lens_ipow(dx, 2)*dy + -0.000144193 *x*lens_ipow(y, 3)*lens_ipow(dx, 2) + -5.54933e-05 *lens_ipow(x, 3)*y*lens_ipow(dy, 2) + 25.0232 *dx*lens_ipow(dy, 3)*lens_ipow(lambda, 4) + -3.58501e-07 *x*lens_ipow(y, 5)*lens_ipow(dy, 2) + 4.07165e-08 *lens_ipow(x, 3)*lens_ipow(y, 4)*dy*lambda + -160515 *lens_ipow(dx, 5)*lens_ipow(dy, 5) + 1.21852e-10 *lens_ipow(x, 7)*lens_ipow(y, 2)*dy + 9.72702e-11 *lens_ipow(x, 8)*y*dx+0.0f;
const float dx33 =  + -0.387565  + 1.17963 *lens_ipow(dy, 2) + 0.379793 *lens_ipow(dx, 2) + -0.000185387 *lens_ipow(y, 2) + -0.000210542 *lens_ipow(x, 2) + 0.0173095 *lens_ipow(lambda, 4) + -0.00253023 *x*y*dx*dy + 4.3721e-05 *lens_ipow(x, 3)*dx + -0.00265214 *x*lens_ipow(y, 2)*lens_ipow(dx, 3) + 2.11365e-06 *lens_ipow(x, 2)*lens_ipow(y, 3)*dy + -0.000110987 *lens_ipow(x, 3)*y*dx*dy + 37.5348 *lens_ipow(dx, 2)*lens_ipow(dy, 2)*lens_ipow(lambda, 4) + 3.44018e-06 *lens_ipow(y, 5)*lens_ipow(dy, 3) + -7.17003e-07 *x*lens_ipow(y, 5)*dx*dy + 5.90163e-07 *lens_ipow(x, 6)*lens_ipow(dy, 2) + -0.618722 *y*lens_ipow(dy, 3)*lens_ipow(lambda, 5) + 2.82305e-06 *lens_ipow(y, 5)*dy*lens_ipow(lambda, 3) + 4.07165e-08 *lens_ipow(x, 3)*lens_ipow(y, 4)*dx*lambda + -133762 *lens_ipow(dx, 6)*lens_ipow(dy, 4) + -0.000591273 *lens_ipow(x, 2)*lens_ipow(lambda, 8) + 1.21852e-10 *lens_ipow(x, 7)*lens_ipow(y, 2)*dx+0.0f;
const float dx34 =  + 0.0154476 *y*lens_ipow(lambda, 2) + 0.0692381 *dy*lens_ipow(lambda, 3) + 50.0464 *lens_ipow(dx, 2)*lens_ipow(dy, 3)*lens_ipow(lambda, 3) + -0.773403 *y*lens_ipow(dy, 4)*lens_ipow(lambda, 4) + 4.23457e-06 *lens_ipow(y, 5)*lens_ipow(dy, 2)*lens_ipow(lambda, 2) + 4.07165e-08 *lens_ipow(x, 3)*lens_ipow(y, 4)*dx*dy + -4.04455e-07 *lens_ipow(x, 4)*y*lens_ipow(lambda, 4) + -0.247506 *y*lens_ipow(lambda, 9) + -5.64749e-07 *lens_ipow(y, 5)*lens_ipow(lambda, 5) + -0.00473018 *lens_ipow(x, 2)*dy*lens_ipow(lambda, 7)+0.0f;
const float dx40 =  + 1.81291e-07  + -0.00104348 *dx + -1.1659e-07 *y + -2.25321e-05 *x + -0.0469432 *dx*lens_ipow(dy, 2) + -0.024014 *lens_ipow(dx, 3) + -0.000889007 *y*dx*dy + -2.5788e-07 *x*lens_ipow(y, 2) + -2.67111e-07 *lens_ipow(x, 3) + -0.00028129 *lens_ipow(y, 2)*lens_ipow(dx, 3) + -0.0341913 *x*lens_ipow(dy, 4) + -0.00136933 *x*y*lens_ipow(dy, 3) + -8.3488e-05 *lens_ipow(x, 2)*lens_ipow(dx, 3) + 3.12386e-06 *lens_ipow(x, 3)*lens_ipow(dy, 2) + 1.22364 *y*lens_ipow(dx, 3)*lens_ipow(dy, 3) + -0.00168603 *x*lens_ipow(y, 2)*lens_ipow(dx, 2)*lens_ipow(dy, 2) + -1.38023e-07 *x*lens_ipow(y, 4)*lens_ipow(dy, 2) + -3.98577e-08 *lens_ipow(x, 3)*lens_ipow(y, 2)*lens_ipow(dx, 2) + -286.151 *lens_ipow(dx, 7)*lens_ipow(dy, 2) + -0.0658487 *lens_ipow(y, 2)*dx*lens_ipow(dy, 6) + 0.00531327 *lens_ipow(y, 3)*lens_ipow(dx, 5)*dy + 0.00978524 *x*y*lens_ipow(dx, 6)*dy + 5.05141e-09 *lens_ipow(x, 2)*lens_ipow(y, 5)*dx*dy+0.0f;
const float dx41 =  + -0.00124334 *dy + -2.18219e-05 *y + -1.1659e-07 *x + -0.0470381 *lens_ipow(dx, 2)*dy + -3.13596e-07 *lens_ipow(y, 3) + -0.000889007 *x*dx*dy + -2.5788e-07 *lens_ipow(x, 2)*y + -0.00059429 *lens_ipow(y, 2)*lens_ipow(dy, 3) + -4.93415e-06 *lens_ipow(y, 3)*lens_ipow(dx, 2) + -0.000562581 *x*y*lens_ipow(dx, 3) + -0.000684665 *lens_ipow(x, 2)*lens_ipow(dy, 3) + -0.32826 *y*lens_ipow(dx, 6) + 1.22364 *x*lens_ipow(dx, 3)*lens_ipow(dy, 3) + -0.00168603 *lens_ipow(x, 2)*y*lens_ipow(dx, 2)*lens_ipow(dy, 2) + -2.76045e-07 *lens_ipow(x, 2)*lens_ipow(y, 3)*lens_ipow(dy, 2) + -1.99288e-08 *lens_ipow(x, 4)*y*lens_ipow(dx, 2) + -0.131697 *x*y*dx*lens_ipow(dy, 6) + 0.0159398 *x*lens_ipow(y, 2)*lens_ipow(dx, 5)*dy + 0.00489262 *lens_ipow(x, 2)*lens_ipow(dx, 6)*dy + 8.41902e-09 *lens_ipow(x, 3)*lens_ipow(y, 4)*dx*dy+0.0f;
const float dx42 =  + -0.115435 *dx + -0.00104348 *x + -0.0940763 *y*dx*dy + -0.0469432 *x*lens_ipow(dy, 2) + -0.0720421 *x*lens_ipow(dx, 2) + -0.000889007 *x*y*dy + -119.94 *dx*lens_ipow(dy, 4) + -382.521 *lens_ipow(dx, 3)*lens_ipow(dy, 2) + -116.504 *lens_ipow(dx, 5) + -2.46707e-06 *lens_ipow(y, 4)*dx + -0.000843871 *x*lens_ipow(y, 2)*lens_ipow(dx, 2) + -8.3488e-05 *lens_ipow(x, 3)*lens_ipow(dx, 2) + -0.98478 *lens_ipow(y, 2)*lens_ipow(dx, 5) + 3.67092 *x*y*lens_ipow(dx, 2)*lens_ipow(dy, 3) + -0.00168603 *lens_ipow(x, 2)*lens_ipow(y, 2)*dx*lens_ipow(dy, 2) + -1.99288e-08 *lens_ipow(x, 4)*lens_ipow(y, 2)*dx + -2003.06 *x*lens_ipow(dx, 6)*lens_ipow(dy, 2) + -0.0658487 *x*lens_ipow(y, 2)*lens_ipow(dy, 6) + 0.0265664 *x*lens_ipow(y, 3)*lens_ipow(dx, 4)*dy + 0.0293557 *lens_ipow(x, 2)*y*lens_ipow(dx, 5)*dy + 1.6838e-09 *lens_ipow(x, 3)*lens_ipow(y, 5)*dy+0.0f;
const float dx43 =  + -0.109123 *dy + -0.00124334 *y + -0.0470381 *y*lens_ipow(dx, 2) + -0.0938864 *x*dx*dy + -0.000889007 *x*y*dx + -91.42 *lens_ipow(dy, 5) + -239.88 *lens_ipow(dx, 2)*lens_ipow(dy, 3) + -191.261 *lens_ipow(dx, 4)*dy + -0.00059429 *lens_ipow(y, 3)*lens_ipow(dy, 2) + -0.0683825 *lens_ipow(x, 2)*lens_ipow(dy, 3) + -0.002054 *lens_ipow(x, 2)*y*lens_ipow(dy, 2) + 1.56193e-06 *lens_ipow(x, 4)*dy + 3.67092 *x*y*lens_ipow(dx, 3)*lens_ipow(dy, 2) + -0.00168603 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(dx, 2)*dy + -1.38023e-07 *lens_ipow(x, 2)*lens_ipow(y, 4)*dy + -572.302 *x*lens_ipow(dx, 7)*dy + -0.395092 *x*lens_ipow(y, 2)*dx*lens_ipow(dy, 5) + 0.00531327 *x*lens_ipow(y, 3)*lens_ipow(dx, 5) + 0.00489262 *lens_ipow(x, 2)*y*lens_ipow(dx, 6) + 1.6838e-09 *lens_ipow(x, 3)*lens_ipow(y, 5)*dx+0.0f;
const float dx44 =  + 0.246012  + -0.604002 *lens_ipow(lambda, 2) + 3.9107 *lens_ipow(lambda, 10)+0.0f;
