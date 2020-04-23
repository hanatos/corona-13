const float dx00 =  + 0.747087  + 5.0693e-05 *dy + 1.0533 *lens_ipow(dy, 2) + 2.60023 *lens_ipow(dx, 2) + 0.0246052 *y*dy + 0.000119221 *lens_ipow(y, 2) + 0.0695846 *x*dx + 0.000364638 *lens_ipow(x, 2) + 0.0690261 *lens_ipow(lambda, 4) + -3.47823e-06 *x*lens_ipow(y, 2)*dx + 6.22609 *lens_ipow(dx, 4)*lens_ipow(dy, 2) + 1.04057e-07 *lens_ipow(y, 4)*lens_ipow(dy, 2) + 3.13269e-07 *lens_ipow(y, 4)*lens_ipow(dx, 2) + 0.286426 *x*dx*lens_ipow(dy, 4) + -0.00777629 *x*y*lens_ipow(dx, 3)*dy + -0.00634012 *lens_ipow(x, 2)*lens_ipow(dx, 4) + 1.60793e-06 *lens_ipow(x, 4)*lens_ipow(dy, 2) + -5.15699e-08 *lens_ipow(x, 5)*dx + -0.0764266 *lens_ipow(dy, 2)*lens_ipow(lambda, 5) + -3.28726e-08 *lens_ipow(x, 2)*lens_ipow(y, 3)*dy*lens_ipow(lambda, 2) + -1.2732e-12 *lens_ipow(x, 4)*lens_ipow(y, 4) + -0.33329 *lens_ipow(lambda, 10) + -5.03113e-15 *lens_ipow(x, 10)+0.0f;
const float dx01 =  + -1.33622e-05 *lambda + 0.000150834 *dy + 1.53717 *dx*dy + 0.0195216 *y*dx + 0.0246052 *x*dy + 0.000238442 *x*y + -3.47823e-06 *lens_ipow(x, 2)*y*dx + -0.00674175 *lens_ipow(y, 2)*lens_ipow(dx, 3)*dy + 4.16227e-07 *x*lens_ipow(y, 3)*lens_ipow(dy, 2) + 1.25308e-06 *x*lens_ipow(y, 3)*lens_ipow(dx, 2) + -0.00388815 *lens_ipow(x, 2)*lens_ipow(dx, 3)*dy + -0.995897 *y*lens_ipow(dx, 7) + -3.28726e-08 *lens_ipow(x, 3)*lens_ipow(y, 2)*dy*lens_ipow(lambda, 2) + -1.01856e-12 *lens_ipow(x, 5)*lens_ipow(y, 3)+0.0f;
const float dx02 =  + 84.6305  + 0.00705907 *dx + 46.1915 *lens_ipow(dy, 2) + 138.576 *lens_ipow(dx, 2) + 1.53717 *y*dy + 0.00976078 *lens_ipow(y, 2) + 5.20046 *x*dx + 0.0347923 *lens_ipow(x, 2) + 3.61538 *lens_ipow(lambda, 3) + -1.73912e-06 *lens_ipow(x, 2)*lens_ipow(y, 2) + -0.00674175 *lens_ipow(y, 3)*lens_ipow(dx, 2)*dy + 24.9043 *x*lens_ipow(dx, 3)*lens_ipow(dy, 2) + 6.26538e-07 *x*lens_ipow(y, 4)*dx + 0.143213 *lens_ipow(x, 2)*lens_ipow(dy, 4) + -0.0116644 *lens_ipow(x, 2)*y*lens_ipow(dx, 2)*dy + -0.00845349 *lens_ipow(x, 3)*lens_ipow(dx, 3) + -8.59498e-09 *lens_ipow(x, 6) + -3.48564 *lens_ipow(y, 2)*lens_ipow(dx, 6) + -19.3704 *lens_ipow(lambda, 10) + -671774 *lens_ipow(dx, 4)*lens_ipow(dy, 6) + -209570 *lens_ipow(dx, 10)+0.0f;
const float dx03 =  + 1.85033e-05  + 0.000150834 *y + 5.0693e-05 *x + 0.0393372 *dy*lambda + 92.383 *dx*dy + 1.53717 *y*dx + 2.10661 *x*dy + 0.0246052 *x*y + -0.00224725 *lens_ipow(y, 3)*lens_ipow(dx, 3) + 12.4522 *x*lens_ipow(dx, 4)*dy + 2.08113e-07 *x*lens_ipow(y, 4)*dy + 0.572851 *lens_ipow(x, 2)*dx*lens_ipow(dy, 3) + -0.00388815 *lens_ipow(x, 2)*y*lens_ipow(dx, 3) + 6.43172e-07 *lens_ipow(x, 5)*dy + -0.152853 *x*dy*lens_ipow(lambda, 5) + -1.09575e-08 *lens_ipow(x, 3)*lens_ipow(y, 3)*lens_ipow(lambda, 2) + -806128 *lens_ipow(dx, 5)*lens_ipow(dy, 5)+0.0f;
const float dx04 =  + -1.33622e-05 *y + 0.0196686 *lens_ipow(dy, 2) + 10.8462 *dx*lens_ipow(lambda, 2) + 0.276104 *x*lens_ipow(lambda, 3) + -0.382133 *x*lens_ipow(dy, 2)*lens_ipow(lambda, 4) + -2.1915e-08 *lens_ipow(x, 3)*lens_ipow(y, 3)*dy*lambda + -193.704 *dx*lens_ipow(lambda, 9) + -3.3329 *x*lens_ipow(lambda, 9)+0.0f;
const float dx10 =  + -2.20108e-06  + 1.51862 *dx*dy + -0.000320125 *lens_ipow(dx, 2) + 1.75895e-05 *y*dy + 0.0242141 *y*dx + 2.68061e-07 *lens_ipow(y, 2) + 0.0192815 *x*dy + 0.000237408 *x*y + -0.000194699 *lens_ipow(y, 2)*dx*dy + 0.000281407 *x*y*lens_ipow(dy, 2)*lens_ipow(lambda, 2) + -1.16441e-08 *x*lens_ipow(y, 4)*dy + 5.21073e-05 *lens_ipow(x, 3)*lens_ipow(dy, 5) + -9.58537e-13 *x*lens_ipow(y, 7)*lens_ipow(lambda, 2)+0.0f;
const float dx11 =  + 0.746666  + 0.000107669 *dx + -1.10446e-06 *y + 2.57977 *lens_ipow(dy, 2) + 1.03579 *lens_ipow(dx, 2) + 0.0667967 *y*dy + 0.000293398 *lens_ipow(y, 2) + 1.75895e-05 *x*dy + 0.0242141 *x*dx + 5.36121e-07 *x*y + 0.000118704 *lens_ipow(x, 2) + 0.0397129 *lens_ipow(lambda, 3) + 0.000134021 *lens_ipow(y, 2)*lambda + 0.00650403 *y*dy*lens_ipow(lambda, 2) + -0.000389397 *x*y*dx*dy + 0.205533 *y*lens_ipow(dy, 5) + 1.67755e-06 *lens_ipow(y, 4)*lens_ipow(dx, 2) + -1.16413e-08 *lens_ipow(y, 5)*dy + 0.000140704 *lens_ipow(x, 2)*lens_ipow(dy, 2)*lens_ipow(lambda, 2) + -2.32883e-08 *lens_ipow(x, 2)*lens_ipow(y, 3)*dy + 0.126423 *y*lens_ipow(dx, 2)*lens_ipow(dy, 3)*lambda + -0.000147472 *lens_ipow(y, 2)*lens_ipow(lambda, 5) + -0.208208 *lens_ipow(lambda, 10) + -0.0257642 *y*dy*lens_ipow(lambda, 8) + -3.35488e-12 *lens_ipow(x, 2)*lens_ipow(y, 6)*lens_ipow(lambda, 2)+0.0f;
const float dx12 =  + 0.0104643 *dy + 0.000107669 *y + 90.8995 *dx*dy + 2.07157 *y*dx + 1.51862 *x*dy + -0.00064025 *x*dx + 0.0242141 *x*y + -0.000194699 *x*lens_ipow(y, 2)*dy + 6.71019e-07 *lens_ipow(y, 5)*dx + 0.126423 *lens_ipow(y, 2)*dx*lens_ipow(dy, 3)*lambda+0.0f;
const float dx13 =  + 84.7872  + 0.00372818 *dy + 0.0104643 *dx + 137.292 *lens_ipow(dy, 2) + 45.4497 *lens_ipow(dx, 2) + 5.15954 *y*dy + 0.0333983 *lens_ipow(y, 2) + 1.51862 *x*dx + 1.75895e-05 *x*y + 0.00964075 *lens_ipow(x, 2) + 5.08034 *lens_ipow(lambda, 4) + 0.00325202 *lens_ipow(y, 2)*lens_ipow(lambda, 2) + -0.000194699 *x*lens_ipow(y, 2)*dx + 0.513833 *lens_ipow(y, 2)*lens_ipow(dy, 4) + -1.94021e-09 *lens_ipow(y, 6) + 0.000281407 *lens_ipow(x, 2)*y*dy*lens_ipow(lambda, 2) + -5.82207e-09 *lens_ipow(x, 2)*lens_ipow(y, 4) + 0.189635 *lens_ipow(y, 2)*lens_ipow(dx, 2)*lens_ipow(dy, 2)*lambda + 6.51342e-05 *lens_ipow(x, 4)*lens_ipow(dy, 4) + -23.5022 *lens_ipow(lambda, 10) + -194.141 *lens_ipow(dy, 4)*lens_ipow(lambda, 6) + -231291 *lens_ipow(dy, 10) + -0.0128821 *lens_ipow(y, 2)*lens_ipow(lambda, 8)+0.0f;
const float dx14 =  + 0.119139 *y*lens_ipow(lambda, 2) + 4.46737e-05 *lens_ipow(y, 3) + 20.3213 *dy*lens_ipow(lambda, 3) + 0.00650403 *lens_ipow(y, 2)*dy*lambda + 0.000281407 *lens_ipow(x, 2)*y*lens_ipow(dy, 2)*lambda + 0.0632117 *lens_ipow(y, 2)*lens_ipow(dx, 2)*lens_ipow(dy, 3) + -0.000245787 *lens_ipow(y, 3)*lens_ipow(lambda, 4) + -235.022 *dy*lens_ipow(lambda, 9) + -232.969 *lens_ipow(dy, 5)*lens_ipow(lambda, 5) + -2.08208 *y*lens_ipow(lambda, 9) + -0.103057 *lens_ipow(y, 2)*dy*lens_ipow(lambda, 7) + -9.58537e-13 *lens_ipow(x, 2)*lens_ipow(y, 7)*lambda+0.0f;
const float dx20 =  + 0.00234223  + 9.7763e-06 *dy + -7.70532e-08 *y + 0.0619327 *lens_ipow(dy, 2) + 0.176116 *lens_ipow(dx, 2) + 0.00136958 *y*dy + 7.91819e-06 *lens_ipow(y, 2) + 1.26598e-06 *x*dy + 0.00412064 *x*dx + 2.37183e-05 *lens_ipow(x, 2) + -2.24639e-06 *y*dx*lambda + 1.0506e-05 *x*lens_ipow(dy, 2) + 1.62664e-05 *x*dx*dy + -0.00388383 *lens_ipow(lambda, 4) + -4.06258e-08 *lens_ipow(y, 3)*dx*dy + -1.60534e-07 *x*lens_ipow(y, 2)*dx*lens_ipow(lambda, 2) + 7.49763e-05 *x*y*lens_ipow(dx, 5)*dy + -2.60261e-09 *lens_ipow(x, 5)*dx*lens_ipow(lambda, 2) + 0.0194614 *lens_ipow(lambda, 10) + -0.000246845 *lens_ipow(x, 2)*lens_ipow(dx, 2)*lens_ipow(lambda, 6)+0.0f;
const float dx21 =  + -1.28695e-06  + -7.70532e-08 *x + 0.114401 *dx*dy + 3.07983e-07 *y*lambda + 0.00137019 *y*dx + 0.00136958 *x*dy + 1.58364e-05 *x*y + -2.24639e-06 *x*dx*lambda + 0.000388753 *y*lens_ipow(dx, 3) + -1.21877e-07 *x*lens_ipow(y, 2)*dx*dy + -0.0148125 *lens_ipow(dx, 6) + -1.60534e-07 *lens_ipow(x, 2)*y*dx*lens_ipow(lambda, 2) + 0.0664658 *y*lens_ipow(dx, 3)*lens_ipow(dy, 4) + 3.74882e-05 *lens_ipow(x, 2)*lens_ipow(dx, 5)*dy+0.0f;
const float dx22 =  + 1.61554  + -0.000527854 *dx + 6.01015 *lens_ipow(dy, 2) + 17.7322 *lens_ipow(dx, 2) + 0.114401 *y*dy + 0.000685094 *lens_ipow(y, 2) + 0.352233 *x*dx + 0.00206032 *lens_ipow(x, 2) + -0.371705 *lens_ipow(lambda, 3) + -2.24639e-06 *x*y*lambda + 8.13322e-06 *lens_ipow(x, 2)*dy + 1.12403 *lens_ipow(dy, 4) + 6.11123 *lens_ipow(dx, 2)*lens_ipow(dy, 2) + 11.9812 *lens_ipow(dx, 4) + 0.000583129 *lens_ipow(y, 2)*lens_ipow(dx, 2) + -4.06258e-08 *x*lens_ipow(y, 3)*dy + -0.3833 *lens_ipow(dy, 2)*lens_ipow(lambda, 4) + -0.0888748 *y*lens_ipow(dx, 5) + -8.02669e-08 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(lambda, 2) + -50.0668 *lens_ipow(dx, 4)*lens_ipow(lambda, 4) + 0.0996987 *lens_ipow(y, 2)*lens_ipow(dx, 2)*lens_ipow(dy, 4) + 0.000187441 *lens_ipow(x, 2)*y*lens_ipow(dx, 4)*dy + -4.33768e-10 *lens_ipow(x, 6)*lens_ipow(lambda, 2) + 2.05805 *lens_ipow(lambda, 10) + -0.000164563 *lens_ipow(x, 3)*dx*lens_ipow(lambda, 6)+0.0f;
const float dx23 =  + -0.000101673  + 9.7763e-06 *x + 12.0203 *dx*dy + 0.114401 *y*dx + 0.123865 *x*dy + 0.00136958 *x*y + 6.3299e-07 *lens_ipow(x, 2) + 1.0506e-05 *lens_ipow(x, 2)*dy + 8.13322e-06 *lens_ipow(x, 2)*dx + 4.49611 *dx*lens_ipow(dy, 3) + 4.07415 *lens_ipow(dx, 3)*dy + -4.06258e-08 *x*lens_ipow(y, 3)*dx + -0.766599 *dx*dy*lens_ipow(lambda, 4) + 0.132932 *lens_ipow(y, 2)*lens_ipow(dx, 3)*lens_ipow(dy, 3) + 3.74882e-05 *lens_ipow(x, 2)*y*lens_ipow(dx, 5)+0.0f;
const float dx24 =  + 1.53991e-07 *lens_ipow(y, 2) + -1.11511 *dx*lens_ipow(lambda, 2) + -2.24639e-06 *x*y*dx + -0.0155353 *x*lens_ipow(lambda, 3) + -1.5332 *dx*lens_ipow(dy, 2)*lens_ipow(lambda, 3) + -1.60534e-07 *lens_ipow(x, 2)*lens_ipow(y, 2)*dx*lambda + -40.0535 *lens_ipow(dx, 5)*lens_ipow(lambda, 3) + -8.67536e-10 *lens_ipow(x, 6)*dx*lambda + 20.5805 *dx*lens_ipow(lambda, 9) + 0.194614 *x*lens_ipow(lambda, 9) + -0.000493689 *lens_ipow(x, 3)*lens_ipow(dx, 2)*lens_ipow(lambda, 5)+0.0f;
const float dx30 =  + -8.1798e-07  + 9.88781e-06 *dy + 7.47659e-08 *y + 0.11402 *dx*dy + -6.23818e-07 *y*dy + 0.00136732 *y*dx + 0.00138292 *x*dy + 1.5824e-05 *x*y + -0.000182479 *dx*lens_ipow(dy, 2) + 4.38151e-05 *x*dx*lens_ipow(dy, 2) + -2.93552e-07 *x*lens_ipow(y, 2)*dy + 0.0100929 *x*lens_ipow(dy, 5) + -7.28436e-12 *x*lens_ipow(y, 5) + 1.85331e-08 *lens_ipow(x, 3)*y*lens_ipow(dx, 2) + 2.12372 *dx*lens_ipow(dy, 7) + -6.14608e-07 *lens_ipow(x, 4)*dx*lens_ipow(dy, 3)+0.0f;
const float dx31 =  + 0.00229889  + 7.47659e-08 *x + 0.177488 *lens_ipow(dy, 2) + 0.0617436 *lens_ipow(dx, 2) + 0.00414549 *y*dy + 2.38243e-05 *lens_ipow(y, 2) + -6.23818e-07 *x*dy + 0.00136732 *x*dx + 7.91198e-06 *lens_ipow(x, 2) + 0.00176033 *y*lens_ipow(dx, 2)*dy + -2.93552e-07 *lens_ipow(x, 2)*y*dy + -0.00666373 *lens_ipow(lambda, 5) + 1.04694e-07 *lens_ipow(y, 4)*lens_ipow(dx, 2) + -1.82109e-11 *lens_ipow(x, 2)*lens_ipow(y, 4) + 4.63327e-09 *lens_ipow(x, 4)*lens_ipow(dx, 2) + -1.54611e-09 *lens_ipow(y, 5)*dy*lambda + -2.31679 *lens_ipow(dy, 8) + 0.0277972 *lens_ipow(lambda, 10) + 0.0105427 *y*lens_ipow(dy, 3)*lens_ipow(lambda, 6)+0.0f;
const float dx32 =  + -5.74518e-05  + 12.0286 *dx*dy + 0.123487 *y*dx + 0.11402 *x*dy + 0.00136732 *x*y + -0.0123546 *lens_ipow(dy, 3) + -0.000182479 *x*lens_ipow(dy, 2) + 0.00176033 *lens_ipow(y, 2)*dx*dy + 2.19076e-05 *lens_ipow(x, 2)*lens_ipow(dy, 2) + -0.847652 *dx*dy*lens_ipow(lambda, 4) + 61.1138 *lens_ipow(dx, 5)*dy + 4.18776e-08 *lens_ipow(y, 5)*dx + 9.26654e-09 *lens_ipow(x, 4)*y*dx + 985.999 *lens_ipow(dx, 3)*lens_ipow(dy, 5) + 2.12372 *x*lens_ipow(dy, 7) + -1.22922e-07 *lens_ipow(x, 5)*lens_ipow(dy, 3)+0.0f;
const float dx33 =  + 1.61351  + 9.88781e-06 *x + 18.3253 *lens_ipow(dy, 2) + 6.01428 *lens_ipow(dx, 2) + 0.354977 *y*dy + 0.00207274 *lens_ipow(y, 2) + 0.11402 *x*dx + -6.23818e-07 *x*y + 0.000691459 *lens_ipow(x, 2) + -0.373649 *lens_ipow(lambda, 3) + -0.0370639 *dx*lens_ipow(dy, 2) + -0.000364958 *x*dx*dy + 0.000880163 *lens_ipow(y, 2)*lens_ipow(dx, 2) + 4.38151e-05 *lens_ipow(x, 2)*dx*dy + -1.46776e-07 *lens_ipow(x, 2)*lens_ipow(y, 2) + -0.423826 *lens_ipow(dx, 2)*lens_ipow(lambda, 4) + 10.1856 *lens_ipow(dx, 6) + 0.0252323 *lens_ipow(x, 2)*lens_ipow(dy, 4) + -2.57684e-10 *lens_ipow(y, 6)*lambda + -58.7129 *lens_ipow(dy, 4)*lens_ipow(lambda, 4) + 1232.5 *lens_ipow(dx, 4)*lens_ipow(dy, 4) + -18.5343 *y*lens_ipow(dy, 7) + 14.866 *x*dx*lens_ipow(dy, 6) + -3.68765e-07 *lens_ipow(x, 5)*dx*lens_ipow(dy, 2) + 2.07253 *lens_ipow(lambda, 10) + 0.015814 *lens_ipow(y, 2)*lens_ipow(dy, 2)*lens_ipow(lambda, 6)+0.0f;
const float dx34 =  + -1.12095 *dy*lens_ipow(lambda, 2) + -0.0333187 *y*lens_ipow(lambda, 4) + -1.6953 *lens_ipow(dx, 2)*dy*lens_ipow(lambda, 3) + -2.57684e-10 *lens_ipow(y, 6)*dy + -46.9703 *lens_ipow(dy, 5)*lens_ipow(lambda, 3) + 20.7253 *dy*lens_ipow(lambda, 9) + 0.277972 *y*lens_ipow(lambda, 9) + 0.031628 *lens_ipow(y, 2)*lens_ipow(dy, 3)*lens_ipow(lambda, 5)+0.0f;
const float dx40 =  + -3.09353e-07  + -0.00233739 *dx + 1.01446e-07 *y + -2.57919e-05 *x + -3.34568e-07 *y*dx + 9.72458e-05 *y*dx*dy + -3.32348e-08 *x*lens_ipow(y, 2) + 0.000207728 *dx*lens_ipow(lambda, 4) + 0.00020672 *x*lens_ipow(dy, 2)*lens_ipow(lambda, 2) + 4.83398e-07 *x*lens_ipow(y, 2)*lens_ipow(dy, 2) + 1.1622e-06 *lens_ipow(x, 2)*y*dx*dy + 1.41263e-06 *lens_ipow(x, 3)*lens_ipow(dx, 2) + -1.52267e-10 *lens_ipow(x, 5) + 0.250315 *y*dx*lens_ipow(dy, 7) + 0.151795 *y*lens_ipow(dx, 7)*dy + -1.07893e-11 *lens_ipow(x, 5)*lens_ipow(y, 2)*lens_ipow(dx, 2) + 0.025646 *x*lens_ipow(dx, 6)*lens_ipow(lambda, 3) + -7.05886e-12 *x*lens_ipow(y, 6)*lens_ipow(dy, 2)*lambda + -0.000198842 *lens_ipow(x, 3)*lens_ipow(dy, 6)*lambda+0.0f;
const float dx41 =  + 1.92069e-07  + -0.00230642 *dy + 6.99894e-06 *dx + -2.58289e-05 *y + 1.01446e-07 *x + -2.35388e-05 *dx*dy + -3.34568e-07 *x*dx + 0.000125198 *y*lens_ipow(dx, 2) + 9.72458e-05 *x*dx*dy + -3.32348e-08 *lens_ipow(x, 2)*y + 0.000228339 *dy*lens_ipow(lambda, 4) + 1.41476e-06 *lens_ipow(y, 3)*lens_ipow(dy, 2) + -1.47383e-10 *lens_ipow(y, 5) + 4.83398e-07 *lens_ipow(x, 2)*y*lens_ipow(dy, 2) + 3.874e-07 *lens_ipow(x, 3)*dx*dy + 0.250315 *x*dx*lens_ipow(dy, 7) + 0.151795 *x*lens_ipow(dx, 7)*dy + -3.59644e-12 *lens_ipow(x, 6)*y*lens_ipow(dx, 2) + 0.0262395 *y*lens_ipow(dy, 6)*lens_ipow(lambda, 3) + -2.11766e-11 *lens_ipow(x, 2)*lens_ipow(y, 5)*lens_ipow(dy, 2)*lambda+0.0f;
const float dx42 =  + -0.234879 *dx + 6.99894e-06 *y + -0.00233739 *x + -2.35388e-05 *y*dy + -3.34568e-07 *x*y + -1.76886 *dx*lens_ipow(dy, 2) + 0.000125198 *lens_ipow(y, 2)*dx + 9.72458e-05 *x*y*dy + -19.201 *dx*lens_ipow(dy, 4) + -39.2521 *lens_ipow(dx, 3)*lens_ipow(dy, 2) + -41.5799 *lens_ipow(dx, 5) + 0.000207728 *x*lens_ipow(lambda, 4) + 3.874e-07 *lens_ipow(x, 3)*y*dy + 7.06313e-07 *lens_ipow(x, 4)*dx + 0.250315 *x*y*lens_ipow(dy, 7) + 1.06257 *x*y*lens_ipow(dx, 6)*dy + -3.59644e-12 *lens_ipow(x, 6)*lens_ipow(y, 2)*dx + 0.076938 *lens_ipow(x, 2)*lens_ipow(dx, 5)*lens_ipow(lambda, 3)+0.0f;
const float dx43 =  + -0.230509 *dy + -0.00230642 *y + -2.35388e-05 *y*dx + -1.76886 *lens_ipow(dx, 2)*dy + 9.72458e-05 *x*y*dx + -43.0087 *lens_ipow(dy, 5) + -38.4019 *lens_ipow(dx, 2)*lens_ipow(dy, 3) + -19.626 *lens_ipow(dx, 4)*dy + 0.000228339 *y*lens_ipow(lambda, 4) + 7.07381e-07 *lens_ipow(y, 4)*dy + 0.00020672 *lens_ipow(x, 2)*dy*lens_ipow(lambda, 2) + 4.83398e-07 *lens_ipow(x, 2)*lens_ipow(y, 2)*dy + 3.874e-07 *lens_ipow(x, 3)*y*dx + 1.75221 *x*y*dx*lens_ipow(dy, 6) + 0.151795 *x*y*lens_ipow(dx, 7) + 0.0787186 *lens_ipow(y, 2)*lens_ipow(dy, 5)*lens_ipow(lambda, 3) + -7.05886e-12 *lens_ipow(x, 2)*lens_ipow(y, 6)*dy*lambda + -0.000298262 *lens_ipow(x, 4)*lens_ipow(dy, 5)*lambda+0.0f;
const float dx44 =  + 0.16531  + -0.416111 *lens_ipow(lambda, 2) + 0.000913357 *y*dy*lens_ipow(lambda, 3) + 0.000830914 *x*dx*lens_ipow(lambda, 3) + 0.00020672 *lens_ipow(x, 2)*lens_ipow(dy, 2)*lambda + 2.77485 *lens_ipow(lambda, 10) + 0.0393593 *lens_ipow(y, 2)*lens_ipow(dy, 6)*lens_ipow(lambda, 2) + 0.038469 *lens_ipow(x, 2)*lens_ipow(dx, 6)*lens_ipow(lambda, 2) + -3.52943e-12 *lens_ipow(x, 2)*lens_ipow(y, 6)*lens_ipow(dy, 2) + -4.97104e-05 *lens_ipow(x, 4)*lens_ipow(dy, 6)+0.0f;
