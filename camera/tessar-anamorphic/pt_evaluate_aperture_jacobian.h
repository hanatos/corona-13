const float dx00 =  + 0.746116  + 5.7837e-05 *dy + 1.17879 *lens_ipow(dy, 2) + 3.00255 *lens_ipow(dx, 2) + 0.025828 *y*dy + 0.000118713 *lens_ipow(y, 2) + 0.0730285 *x*dx + 0.000358858 *lens_ipow(x, 2) + 0.0419482 *lens_ipow(lambda, 3) + 1.74611e-05 *lens_ipow(x, 2)*lens_ipow(lambda, 2) + 0.000769934 *lens_ipow(x, 2)*lens_ipow(dy, 2)*lambda + 1.64312e-05 *lens_ipow(y, 2)*lens_ipow(lambda, 4) + 2.23226e-07 *lens_ipow(y, 4)*lens_ipow(dy, 2) + 0.295307 *x*dx*lens_ipow(dy, 4) + -0.012539 *x*y*lens_ipow(dx, 3)*dy + -5.52512e-08 *lens_ipow(x, 3)*lens_ipow(y, 2)*dx + -4.59973e-10 *lens_ipow(x, 4)*lens_ipow(y, 2) + -2.10946 *lens_ipow(dy, 6)*lambda + 0.0677619 *x*lens_ipow(dx, 3)*lens_ipow(lambda, 3) + -0.226439 *lens_ipow(lambda, 10) + 0.00770434 *y*dy*lens_ipow(lambda, 8)+0.0f;
const float dx01 =  + -3.59444e-06  + 1.79321 *dx*dy + 0.0189432 *y*dx + 0.025828 *x*dy + 0.000237427 *x*y + 0.00418244 *y*dx*lambda + 0.0417051 *dx*dy*lens_ipow(lambda, 4) + 6.98577e-05 *lens_ipow(y, 3)*dx*lens_ipow(dy, 2) + 3.28623e-05 *x*y*lens_ipow(lambda, 4) + 8.92906e-07 *x*lens_ipow(y, 3)*lens_ipow(dy, 2) + -0.00626952 *lens_ipow(x, 2)*lens_ipow(dx, 3)*dy + -2.76256e-08 *lens_ipow(x, 4)*y*dx + -1.83989e-10 *lens_ipow(x, 5)*y + -0.00633286 *y*dx*lens_ipow(lambda, 5) + -1.20888 *y*lens_ipow(dx, 7)*lens_ipow(lambda, 2) + 0.00770434 *x*dy*lens_ipow(lambda, 8)+0.0f;
const float dx02 =  + 88.7549  + 0.00558861 *dx + 61.8435 *lens_ipow(dy, 2) + 187.439 *lens_ipow(dx, 2) + 1.79321 *y*dy + 0.00947158 *lens_ipow(y, 2) + 6.00511 *x*dx + 0.0365142 *lens_ipow(x, 2) + 3.65264 *lens_ipow(lambda, 3) + 0.00209122 *lens_ipow(y, 2)*lambda + 0.0417051 *y*dy*lens_ipow(lambda, 4) + 1.74644e-05 *lens_ipow(y, 4)*lens_ipow(dy, 2) + 0.147654 *lens_ipow(x, 2)*lens_ipow(dy, 4) + -0.0188085 *lens_ipow(x, 2)*y*lens_ipow(dx, 2)*dy + -1.38128e-08 *lens_ipow(x, 4)*lens_ipow(y, 2) + -0.00316643 *lens_ipow(y, 2)*lens_ipow(lambda, 5) + 0.101643 *lens_ipow(x, 2)*lens_ipow(dx, 2)*lens_ipow(lambda, 3) + -36275.5 *lens_ipow(dx, 2)*lens_ipow(dy, 6) + -64216.2 *lens_ipow(dx, 6)*lens_ipow(dy, 2) + -19.1573 *lens_ipow(lambda, 10) + -10969.3 *lens_ipow(dx, 6)*lens_ipow(lambda, 4) + -4.2311 *lens_ipow(y, 2)*lens_ipow(dx, 6)*lens_ipow(lambda, 2)+0.0f;
const float dx03 =  + 0.00786728 *dy + 5.7837e-05 *x + 123.687 *dx*dy + 1.79321 *y*dx + 2.35758 *x*dy + 0.025828 *x*y + 0.000513289 *lens_ipow(x, 3)*dy*lambda + 0.0417051 *y*dx*lens_ipow(lambda, 4) + 3.49289e-05 *lens_ipow(y, 4)*dx*dy + 4.46453e-07 *x*lens_ipow(y, 4)*dy + 0.590615 *lens_ipow(x, 2)*dx*lens_ipow(dy, 3) + -0.00626952 *lens_ipow(x, 2)*y*lens_ipow(dx, 3) + -12.6568 *x*lens_ipow(dy, 5)*lambda + -72551 *lens_ipow(dx, 3)*lens_ipow(dy, 5) + -18347.5 *lens_ipow(dx, 7)*dy + 0.00770434 *x*y*lens_ipow(lambda, 8)+0.0f;
const float dx04 =  + 10.9579 *dx*lens_ipow(lambda, 2) + 0.00209122 *lens_ipow(y, 2)*dx + 0.125845 *x*lens_ipow(lambda, 2) + 1.16407e-05 *lens_ipow(x, 3)*lambda + 0.000256645 *lens_ipow(x, 3)*lens_ipow(dy, 2) + 0.16682 *y*dx*dy*lens_ipow(lambda, 3) + 6.57246e-05 *x*lens_ipow(y, 2)*lens_ipow(lambda, 3) + -0.0158322 *lens_ipow(y, 2)*dx*lens_ipow(lambda, 4) + -2.10946 *x*lens_ipow(dy, 6) + 0.101643 *lens_ipow(x, 2)*lens_ipow(dx, 3)*lens_ipow(lambda, 2) + -191.573 *dx*lens_ipow(lambda, 9) + -6268.16 *lens_ipow(dx, 7)*lens_ipow(lambda, 3) + -1.20888 *lens_ipow(y, 2)*lens_ipow(dx, 7)*lambda + -2.26439 *x*lens_ipow(lambda, 9) + 0.0616347 *x*y*dy*lens_ipow(lambda, 7)+0.0f;
const float dx10 =  + -9.31942e-06  + -7.22249e-07 *y + 1.79209 *dx*dy + 1.52149e-05 *y*dy + 0.0259214 *y*dx + 2.28374e-07 *lens_ipow(y, 2) + 0.0208651 *x*dy + 0.000242133 *x*y + 0.000296887 *dx*lens_ipow(lambda, 2) + 4.52989e-05 *y*dx*dy + 0.0224844 *y*dx*lens_ipow(dy, 2)*lambda + 2.24037e-06 *lens_ipow(x, 3)*y*lens_ipow(dy, 2) + 0.00315905 *lens_ipow(x, 3)*lens_ipow(dx, 2)*lens_ipow(dy, 3) + 3.37255e-06 *lens_ipow(x, 3)*y*lens_ipow(dx, 2)*lens_ipow(lambda, 2) + -4.0154e-13 *lens_ipow(x, 3)*lens_ipow(y, 5) + -5.97977e-13 *lens_ipow(x, 7)*y+0.0f;
const float dx11 =  + 0.747015  + -7.22249e-07 *x + 3.00018 *lens_ipow(dy, 2) + 1.18915 *lens_ipow(dx, 2) + 0.0734816 *y*dy + 0.00036385 *lens_ipow(y, 2) + 1.52149e-05 *x*dy + 0.0259214 *x*dx + 4.56748e-07 *x*y + 0.000121067 *lens_ipow(x, 2) + 4.52989e-05 *x*dx*dy + 0.068935 *lens_ipow(lambda, 4) + 0.0224844 *x*dx*lens_ipow(dy, 2)*lambda + 0.369606 *y*lens_ipow(dx, 2)*lens_ipow(dy, 3) + -0.00794992 *lens_ipow(y, 2)*lens_ipow(dy, 4) + 1.38871e-06 *lens_ipow(y, 4)*lens_ipow(dx, 2) + -5.44427e-08 *lens_ipow(y, 5)*dy + 5.60093e-07 *lens_ipow(x, 4)*lens_ipow(dy, 2) + 8.43138e-07 *lens_ipow(x, 4)*lens_ipow(dx, 2)*lens_ipow(lambda, 2) + -5.01925e-13 *lens_ipow(x, 4)*lens_ipow(y, 4) + -7.47471e-14 *lens_ipow(x, 8) + -15.2218 *lens_ipow(dx, 6)*lens_ipow(lambda, 3) + -0.335914 *lens_ipow(lambda, 10) + -5.00628e-15 *lens_ipow(y, 10)+0.0f;
const float dx12 =  + 0.00496451 *dy + 122.914 *dx*dy + 2.3783 *y*dx + 1.79209 *x*dy + 0.0259214 *x*y + 0.000296887 *x*lens_ipow(lambda, 2) + 4.52989e-05 *x*y*dy + 0.0224844 *x*y*lens_ipow(dy, 2)*lambda + 0.369606 *lens_ipow(y, 2)*dx*lens_ipow(dy, 3) + 5.55484e-07 *lens_ipow(y, 5)*dx + -273.532 *dx*lens_ipow(dy, 3)*lens_ipow(lambda, 4) + 0.00157952 *lens_ipow(x, 4)*dx*lens_ipow(dy, 3) + 1.68628e-06 *lens_ipow(x, 4)*y*dx*lens_ipow(lambda, 2) + -91.3306 *y*lens_ipow(dx, 5)*lens_ipow(lambda, 3) + -178247 *dx*lens_ipow(dy, 9) + -1.47558e+06 *lens_ipow(dx, 5)*lens_ipow(dy, 5)+0.0f;
const float dx13 =  + 88.7444  + 0.00496451 *dx + 184.716 *lens_ipow(dy, 2) + 61.457 *lens_ipow(dx, 2) + 6.00036 *y*dy + 0.0367408 *lens_ipow(y, 2) + 1.79209 *x*dx + 1.52149e-05 *x*y + 0.0104326 *lens_ipow(x, 2) + 3.8145 *lens_ipow(lambda, 3) + 4.52989e-05 *x*y*dx + 0.0449689 *x*y*dx*dy*lambda + 0.554409 *lens_ipow(y, 2)*lens_ipow(dx, 2)*lens_ipow(dy, 2) + -0.0105999 *lens_ipow(y, 3)*lens_ipow(dy, 3) + -9.07379e-09 *lens_ipow(y, 6) + 1.12019e-06 *lens_ipow(x, 4)*y*dy + -410.297 *lens_ipow(dx, 2)*lens_ipow(dy, 2)*lens_ipow(lambda, 4) + 0.00236929 *lens_ipow(x, 4)*lens_ipow(dx, 2)*lens_ipow(dy, 2) + -20.1841 *lens_ipow(lambda, 10) + -317485 *lens_ipow(dy, 10) + -802113 *lens_ipow(dx, 2)*lens_ipow(dy, 8) + -1.22965e+06 *lens_ipow(dx, 6)*lens_ipow(dy, 4)+0.0f;
const float dx14 =  + 11.4435 *dy*lens_ipow(lambda, 2) + 0.000593774 *x*dx*lambda + 0.27574 *y*lens_ipow(lambda, 3) + 0.0224844 *x*y*dx*lens_ipow(dy, 2) + -547.063 *lens_ipow(dx, 2)*lens_ipow(dy, 3)*lens_ipow(lambda, 3) + 1.68628e-06 *lens_ipow(x, 4)*y*lens_ipow(dx, 2)*lambda + -45.6653 *y*lens_ipow(dx, 6)*lens_ipow(lambda, 2) + -201.841 *dy*lens_ipow(lambda, 9) + -3.35914 *y*lens_ipow(lambda, 9)+0.0f;
const float dx20 =  + 0.00229368  + -4.89523e-06 *dx + -6.89611e-08 *y + 0.0699194 *lens_ipow(dy, 2) + 0.200911 *lens_ipow(dx, 2) + 0.00145803 *y*dy + 7.93143e-06 *lens_ipow(y, 2) + 0.00438372 *x*dx + 2.36462e-05 *lens_ipow(x, 2) + 9.82923e-08 *lens_ipow(y, 2)*dx + 0.00145134 *x*dx*lens_ipow(dy, 2) + -0.00641514 *lens_ipow(lambda, 5) + -0.0684174 *lens_ipow(dy, 6) + 7.22507e-09 *lens_ipow(y, 4)*lens_ipow(dy, 2) + -5.25081e-07 *x*lens_ipow(y, 2)*dx*lens_ipow(lambda, 2) + 9.57584e-08 *lens_ipow(x, 4)*lens_ipow(dy, 2) + -1.17912e-11 *lens_ipow(x, 4)*lens_ipow(y, 2) + -8.56841e-07 *lens_ipow(x, 3)*dx*lens_ipow(lambda, 3) + -0.010986 *lens_ipow(dy, 2)*lens_ipow(lambda, 6) + 0.0546309 *x*lens_ipow(dx, 7) + -0.0813517 *lens_ipow(dx, 2)*lens_ipow(lambda, 7) + 0.027659 *lens_ipow(lambda, 10)+0.0f;
const float dx21 =  + -9.35555e-07  + 2.04063e-07 *y + -6.89611e-08 *x + 0.130373 *dx*dy + 0.0014691 *y*dx + 0.00145803 *x*dy + 1.58629e-05 *x*y + 1.96585e-07 *x*y*dx + 2.89003e-08 *x*lens_ipow(y, 3)*lens_ipow(dy, 2) + -5.25081e-07 *lens_ipow(x, 2)*y*dx*lens_ipow(lambda, 2) + -4.71646e-12 *lens_ipow(x, 5)*y + 0.0723921 *y*lens_ipow(dx, 7) + -0.000247787 *lens_ipow(y, 2)*dx*dy*lens_ipow(lambda, 6)+0.0f;
const float dx22 =  + 1.61219  + 0.00059096 *dy + -4.89523e-06 *x + 7.1422 *lens_ipow(dy, 2) + 21.8475 *lens_ipow(dx, 2) + 0.130373 *y*dy + 0.000734551 *lens_ipow(y, 2) + 0.401823 *x*dx + 0.00219186 *lens_ipow(x, 2) + 9.82923e-08 *x*lens_ipow(y, 2) + -0.528422 *lens_ipow(lambda, 4) + 0.902644 *lens_ipow(dy, 4) + 3.1646 *lens_ipow(dx, 2)*lens_ipow(dy, 2) + 0.000725671 *lens_ipow(x, 2)*lens_ipow(dy, 2) + -0.986041 *lens_ipow(dy, 2)*lens_ipow(lambda, 3) + -4.27918 *lens_ipow(dx, 2)*lens_ipow(lambda, 3) + -2.62541e-07 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(lambda, 2) + -2.1421e-07 *lens_ipow(x, 4)*lens_ipow(lambda, 3) + 0.253372 *lens_ipow(y, 2)*lens_ipow(dx, 6) + 0.191208 *lens_ipow(x, 2)*lens_ipow(dx, 6) + -0.162703 *x*dx*lens_ipow(lambda, 7) + 2.54923 *lens_ipow(lambda, 10) + -8.25955e-05 *lens_ipow(y, 3)*dy*lens_ipow(lambda, 6)+0.0f;
const float dx23 =  + -6.80065e-05  + 0.00059096 *dx + 14.2844 *dx*dy + 0.130373 *y*dx + 0.139839 *x*dy + 0.00145803 *x*y + 3.61058 *dx*lens_ipow(dy, 3) + 2.10973 *lens_ipow(dx, 3)*dy + 0.00145134 *lens_ipow(x, 2)*dx*dy + -1.97208 *dx*dy*lens_ipow(lambda, 3) + -0.410505 *x*lens_ipow(dy, 5) + 1.44501e-08 *x*lens_ipow(y, 4)*dy + 3.83034e-08 *lens_ipow(x, 5)*dy + -0.0219721 *x*dy*lens_ipow(lambda, 6) + -8.25955e-05 *lens_ipow(y, 3)*dx*lens_ipow(lambda, 6)+0.0f;
const float dx24 =  + -2.11369 *dx*lens_ipow(lambda, 3) + -2.95812 *dx*lens_ipow(dy, 2)*lens_ipow(lambda, 2) + -4.27918 *lens_ipow(dx, 3)*lens_ipow(lambda, 2) + -0.0320757 *x*lens_ipow(lambda, 4) + -5.25081e-07 *lens_ipow(x, 2)*lens_ipow(y, 2)*dx*lambda + -6.4263e-07 *lens_ipow(x, 4)*dx*lens_ipow(lambda, 2) + -0.0659162 *x*lens_ipow(dy, 2)*lens_ipow(lambda, 5) + -0.569462 *x*lens_ipow(dx, 2)*lens_ipow(lambda, 6) + 25.4923 *dx*lens_ipow(lambda, 9) + -0.000495573 *lens_ipow(y, 3)*dx*dy*lens_ipow(lambda, 5) + 0.27659 *x*lens_ipow(lambda, 9)+0.0f;
const float dx30 =  + -1.69386e-07  + 1.49325e-05 *dy + 8.8508e-08 *y + -2.38833e-05 *lens_ipow(dy, 2) + 0.129823 *dx*dy + 0.0014561 *y*dx + 0.00146051 *x*dy + 1.22003e-06 *x*dx + 1.58471e-05 *x*y + 3.20074e-06 *y*dx*dy + 3.37169e-08 *lens_ipow(x, 3)*y*lens_ipow(dx, 2) + -3.50519e-08 *lens_ipow(x, 4)*dx*dy + 2.60165e-06 *x*y*lens_ipow(dy, 2)*lens_ipow(lambda, 3) + 0.121168 *x*lens_ipow(dy, 7) + 1.59847e-10 *x*lens_ipow(y, 5)*lens_ipow(dy, 2) + -9.77827e-17 *lens_ipow(x, 3)*lens_ipow(y, 7)+0.0f;
const float dx31 =  + 0.0023379  + -6.28315e-06 *dy + -1.94745e-07 *y + 8.8508e-08 *x + 0.200098 *lens_ipow(dy, 2) + 0.0695817 *lens_ipow(dx, 2) + 0.0043666 *y*dy + 2.37044e-05 *lens_ipow(y, 2) + 0.0014561 *x*dx + 7.92354e-06 *lens_ipow(x, 2) + 3.20074e-06 *x*dx*dy + -0.00386528 *lens_ipow(lambda, 4) + 0.0016003 *y*lens_ipow(dx, 2)*dy + -0.000206298 *lens_ipow(y, 2)*lens_ipow(dy, 4) + 9.45083e-08 *lens_ipow(y, 4)*lens_ipow(dx, 2) + 8.42922e-09 *lens_ipow(x, 4)*lens_ipow(dx, 2) + 1.30083e-06 *lens_ipow(x, 2)*lens_ipow(dy, 2)*lens_ipow(lambda, 3) + -0.00783016 *lens_ipow(dx, 2)*lens_ipow(lambda, 6) + 3.99616e-10 *lens_ipow(x, 2)*lens_ipow(y, 4)*lens_ipow(dy, 2) + 0.0194209 *lens_ipow(lambda, 10) + -1.7112e-16 *lens_ipow(x, 4)*lens_ipow(y, 6)+0.0f;
const float dx32 =  + -0.000117436  + 14.1624 *dx*dy + 0.139163 *y*dx + 0.129823 *x*dy + 0.0014561 *x*y + 6.10017e-07 *lens_ipow(x, 2) + 3.20074e-06 *x*y*dy + 0.0016003 *lens_ipow(y, 2)*dx*dy + -1.08749 *dx*dy*lens_ipow(lambda, 3) + 73.7985 *lens_ipow(dx, 5)*dy + 3.78033e-08 *lens_ipow(y, 5)*dx + 1.68584e-08 *lens_ipow(x, 4)*y*dx + -7.01038e-09 *lens_ipow(x, 5)*dy + 975.421 *lens_ipow(dx, 3)*lens_ipow(dy, 5) + -0.0156603 *y*dx*lens_ipow(lambda, 6)+0.0f;
const float dx33 =  + 1.62761  + -6.28315e-06 *y + 1.49325e-05 *x + 21.2793 *lens_ipow(dy, 2) + 7.0812 *lens_ipow(dx, 2) + 0.400196 *y*dy + 0.0021833 *lens_ipow(y, 2) + -4.77667e-05 *x*dy + 0.129823 *x*dx + 0.000730256 *lens_ipow(x, 2) + -0.384908 *lens_ipow(lambda, 3) + 3.20074e-06 *x*y*dx + 0.000800151 *lens_ipow(y, 2)*lens_ipow(dx, 2) + -0.543746 *lens_ipow(dx, 2)*lens_ipow(lambda, 3) + 12.2998 *lens_ipow(dx, 6) + -0.000275064 *lens_ipow(y, 3)*lens_ipow(dy, 3) + -7.01038e-09 *lens_ipow(x, 5)*dx + 2.60165e-06 *lens_ipow(x, 2)*y*dy*lens_ipow(lambda, 3) + -26.0996 *lens_ipow(dy, 4)*lens_ipow(lambda, 4) + 1219.28 *lens_ipow(dx, 4)*lens_ipow(dy, 4) + 0.424087 *lens_ipow(x, 2)*lens_ipow(dy, 6) + 1.59847e-10 *lens_ipow(x, 2)*lens_ipow(y, 5)*dy + 2.16074 *lens_ipow(lambda, 10)+0.0f;
const float dx34 =  + -1.15472 *dy*lens_ipow(lambda, 2) + -0.0154611 *y*lens_ipow(lambda, 3) + -1.63124 *lens_ipow(dx, 2)*dy*lens_ipow(lambda, 2) + 3.90248e-06 *lens_ipow(x, 2)*y*lens_ipow(dy, 2)*lens_ipow(lambda, 2) + -20.8797 *lens_ipow(dy, 5)*lens_ipow(lambda, 3) + -0.046981 *y*lens_ipow(dx, 2)*lens_ipow(lambda, 5) + 21.6074 *dy*lens_ipow(lambda, 9) + 0.194209 *y*lens_ipow(lambda, 9)+0.0f;
const float dx40 =  + -6.19535e-07  + 1.13423e-05 *dy + -0.00247576 *dx + 1.95652e-07 *y + -2.62298e-05 *x + -1.44272e-05 *lens_ipow(dy, 2) + -3.91438e-07 *x*dx + 9.49629e-09 *x*y + 8.25607e-05 *y*dx*dy + -2.73545e-08 *x*lens_ipow(y, 2) + 0.000166838 *x*lens_ipow(dy, 2)*lambda + 0.000211309 *dx*lens_ipow(lambda, 4) + 9.54536e-08 *x*lens_ipow(y, 2)*lens_ipow(dy, 2) + 1.21632e-06 *lens_ipow(x, 2)*y*dx*dy + 1.3657e-06 *lens_ipow(x, 3)*lens_ipow(dx, 2) + -1.32617e-10 *lens_ipow(x, 5) + 0.216067 *y*lens_ipow(dx, 7)*dy + -0.00503713 *lens_ipow(y, 2)*dx*lens_ipow(dy, 6) + -1.14893e-11 *lens_ipow(x, 5)*lens_ipow(y, 2)*lens_ipow(dx, 2) + 0.0308057 *x*lens_ipow(dx, 6)*lens_ipow(lambda, 3) + -0.000198447 *lens_ipow(x, 3)*lens_ipow(dy, 6)*lambda+0.0f;
const float dx41 =  + 7.15764e-08  + -0.00243182 *dy + 1.6837e-05 *dx + -2.59832e-05 *y + 1.95652e-07 *x + 4.74815e-09 *lens_ipow(x, 2) + 0.000127614 *y*lens_ipow(dx, 2) + 8.25607e-05 *x*dx*dy + -2.73545e-08 *lens_ipow(x, 2)*y + 1.3386e-06 *lens_ipow(y, 3)*lens_ipow(dy, 2) + -1.37063e-10 *lens_ipow(y, 5) + 9.54536e-08 *lens_ipow(x, 2)*y*lens_ipow(dy, 2) + 4.05441e-07 *lens_ipow(x, 3)*dx*dy + 0.216067 *x*lens_ipow(dx, 7)*dy + -0.0100743 *x*y*dx*lens_ipow(dy, 6) + -3.82975e-12 *lens_ipow(x, 6)*y*lens_ipow(dx, 2)+0.0f;
const float dx42 =  + 0.00104942 *dy + -0.257885 *dx + 1.6837e-05 *y + -0.00247576 *x + -1.95719e-07 *lens_ipow(x, 2) + -1.82762 *dx*lens_ipow(dy, 2) + 0.000127614 *lens_ipow(y, 2)*dx + 8.25607e-05 *x*y*dy + -21.6453 *dx*lens_ipow(dy, 4) + -50.2974 *lens_ipow(dx, 3)*lens_ipow(dy, 2) + -48.517 *lens_ipow(dx, 5) + 0.000211309 *x*lens_ipow(lambda, 4) + 4.05441e-07 *lens_ipow(x, 3)*y*dy + 6.82852e-07 *lens_ipow(x, 4)*dx + 1.51247 *x*y*lens_ipow(dx, 6)*dy + -0.00503713 *x*lens_ipow(y, 2)*lens_ipow(dy, 6) + -3.82975e-12 *lens_ipow(x, 6)*lens_ipow(y, 2)*dx + 0.0924172 *lens_ipow(x, 2)*lens_ipow(dx, 5)*lens_ipow(lambda, 3)+0.0f;
const float dx43 =  + -0.255751 *dy + 0.00104942 *dx + -0.00243182 *y + 1.13423e-05 *x + -2.88543e-05 *x*dy + -1.82762 *lens_ipow(dx, 2)*dy + 8.25607e-05 *x*y*dx + 0.000166838 *lens_ipow(x, 2)*dy*lambda + -44.8014 *lens_ipow(dy, 5) + -43.2907 *lens_ipow(dx, 2)*lens_ipow(dy, 3) + -25.1487 *lens_ipow(dx, 4)*dy + 6.693e-07 *lens_ipow(y, 4)*dy + 9.54536e-08 *lens_ipow(x, 2)*lens_ipow(y, 2)*dy + 4.05441e-07 *lens_ipow(x, 3)*y*dx + 0.216067 *x*y*lens_ipow(dx, 7) + -0.0302228 *x*lens_ipow(y, 2)*dx*lens_ipow(dy, 5) + -0.000297671 *lens_ipow(x, 4)*lens_ipow(dy, 5)*lambda+0.0f;
const float dx44 =  + 0.164963  + -0.415018 *lens_ipow(lambda, 2) + 8.3419e-05 *lens_ipow(x, 2)*lens_ipow(dy, 2) + 0.000845236 *x*dx*lens_ipow(lambda, 3) + 2.75928 *lens_ipow(lambda, 10) + 0.0462086 *lens_ipow(x, 2)*lens_ipow(dx, 6)*lens_ipow(lambda, 2) + -4.96118e-05 *lens_ipow(x, 4)*lens_ipow(dy, 6)+0.0f;