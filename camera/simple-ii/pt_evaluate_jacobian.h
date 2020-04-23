const float dx00 =  + 0.494773  + 2.1368e-06 *y + 0.00761577 *y*dy + 5.13835e-05 *lens_ipow(y, 2) + 0.0288834 *x*dx + 0.000208539 *lens_ipow(x, 2) + -0.0338336 *lens_ipow(dy, 4) + 0.268621 *lens_ipow(dx, 2)*lens_ipow(dy, 2) + 0.00118461 *lens_ipow(x, 2)*lens_ipow(dx, 2) + 0.0964549 *lens_ipow(lambda, 5) + 0.000238725 *lens_ipow(x, 2)*lens_ipow(dy, 2)*lens_ipow(lambda, 2) + 2.75395e-07 *lens_ipow(y, 4)*lens_ipow(dx, 2)*lens_ipow(dy, 2) + 7.91825e-11 *lens_ipow(x, 5)*lens_ipow(y, 2)*dx + 4.29096e-10 *lens_ipow(x, 7)*dx + -0.0927918 *lens_ipow(dy, 2)*lens_ipow(lambda, 8) + -0.00679733 *y*lens_ipow(dy, 9) + 0.0102953 *y*lens_ipow(dx, 4)*dy*lens_ipow(lambda, 4) + 3.69317e-08 *lens_ipow(y, 4)*lens_ipow(lambda, 6) + -0.00130922 *lens_ipow(x, 2)*lens_ipow(lambda, 8) + -3.88829e-05 *lens_ipow(x, 2)*y*lens_ipow(dx, 6)*dy + 2.75535e-10 *lens_ipow(x, 7)*dx*lens_ipow(dy, 2) + 1.60515e-11 *lens_ipow(x, 8)*lens_ipow(lambda, 2)+0.0f;
const float dx01 =  + 4.18341e-05  + 2.1368e-06 *x + 0.0086323 *y*dx + 0.00761577 *x*dy + 0.000102767 *x*y + 0.204007 *lens_ipow(dx, 3)*dy + -0.00663019 *y*dx*lens_ipow(dy, 2) + 1.10158e-06 *x*lens_ipow(y, 3)*lens_ipow(dx, 2)*lens_ipow(dy, 2) + 2.63942e-11 *lens_ipow(x, 6)*y*dx + 1.13367e-07 *lens_ipow(y, 5)*lens_ipow(dx, 5) + -0.00679733 *x*lens_ipow(dy, 9) + 0.0102953 *x*lens_ipow(dx, 4)*dy*lens_ipow(lambda, 4) + 1.47727e-07 *x*lens_ipow(y, 3)*lens_ipow(lambda, 6) + -1.2961e-05 *lens_ipow(x, 3)*lens_ipow(dx, 6)*dy+0.0f;
const float dx02 =  + 30.421  + -12.3511 *lens_ipow(dy, 2) + -41.2386 *lens_ipow(dx, 2) + 0.00431615 *lens_ipow(y, 2) + 0.0144417 *lens_ipow(x, 2) + 3.27288 *lens_ipow(lambda, 4) + 52.1837 *lens_ipow(dx, 2)*lens_ipow(dy, 2) + 30.3204 *lens_ipow(dx, 4) + 0.61202 *y*lens_ipow(dx, 2)*dy + -0.00331509 *lens_ipow(y, 2)*lens_ipow(dy, 2) + 0.537243 *x*dx*lens_ipow(dy, 2) + 0.000789738 *lens_ipow(x, 3)*dx + 5.74207 *lens_ipow(dy, 6) + -25.9346 *lens_ipow(dx, 4)*lens_ipow(dy, 2)*lambda + -4.74837 *lens_ipow(dx, 2)*lens_ipow(dy, 6) + -6.74876 *lens_ipow(dx, 8) + 5.5079e-07 *x*lens_ipow(y, 4)*dx*lens_ipow(dy, 2) + 1.31971e-11 *lens_ipow(x, 6)*lens_ipow(y, 2) + 5.3637e-11 *lens_ipow(x, 8) + -10.377 *lens_ipow(lambda, 10) + 9.44729e-08 *lens_ipow(y, 6)*lens_ipow(dx, 4) + 0.0411812 *x*y*lens_ipow(dx, 3)*dy*lens_ipow(lambda, 4) + -7.77658e-05 *lens_ipow(x, 3)*y*lens_ipow(dx, 5)*dy + 3.44419e-11 *lens_ipow(x, 8)*lens_ipow(dy, 2)+0.0f;
const float dx03 =  + 0.0105258 *lens_ipow(dy, 2) + -24.7022 *dx*dy + 0.00761577 *x*y + 34.7891 *lens_ipow(dx, 3)*dy + 0.204007 *y*lens_ipow(dx, 3) + -0.00663019 *lens_ipow(y, 2)*dx*dy + -0.135334 *x*lens_ipow(dy, 3) + 0.537243 *x*lens_ipow(dx, 2)*dy + 34.4524 *dx*lens_ipow(dy, 5) + 0.00015915 *lens_ipow(x, 3)*dy*lens_ipow(lambda, 2) + -10.3738 *lens_ipow(dx, 5)*dy*lambda + -9.49674 *lens_ipow(dx, 3)*lens_ipow(dy, 5) + 5.5079e-07 *x*lens_ipow(y, 4)*lens_ipow(dx, 2)*dy + -0.185584 *x*dy*lens_ipow(lambda, 8) + -0.061176 *x*y*lens_ipow(dy, 8) + 0.0102953 *x*y*lens_ipow(dx, 4)*lens_ipow(lambda, 4) + -1.2961e-05 *lens_ipow(x, 3)*y*lens_ipow(dx, 6) + 6.88838e-11 *lens_ipow(x, 8)*dx*dy+0.0f;
const float dx04 =  + 13.0915 *dx*lens_ipow(lambda, 3) + 0.482275 *x*lens_ipow(lambda, 4) + 0.00015915 *lens_ipow(x, 3)*lens_ipow(dy, 2)*lambda + -5.18691 *lens_ipow(dx, 5)*lens_ipow(dy, 2) + -103.77 *dx*lens_ipow(lambda, 9) + -0.742334 *x*lens_ipow(dy, 2)*lens_ipow(lambda, 7) + 0.0411812 *x*y*lens_ipow(dx, 4)*dy*lens_ipow(lambda, 3) + 2.2159e-07 *x*lens_ipow(y, 4)*lens_ipow(lambda, 5) + -0.00349124 *lens_ipow(x, 3)*lens_ipow(lambda, 7) + 3.56699e-12 *lens_ipow(x, 9)*lambda+0.0f;
const float dx10 =  + 0.00775022 *y*dx + 0.00780782 *x*dy + 7.99405e-05 *x*y + -0.00974809 *y*dx*lens_ipow(dy, 2) + 0.00136618 *x*lens_ipow(dy, 5) + -4.16407e-05 *lens_ipow(x, 2)*y*dx*lens_ipow(lambda, 2) + -3.92535e-05 *lens_ipow(x, 3)*lens_ipow(dx, 2)*dy + -7.98515e-13 *x*lens_ipow(y, 7) + -2.24236e-09 *lens_ipow(x, 3)*lens_ipow(y, 3)*lens_ipow(dx, 2) + -4.92814e-12 *lens_ipow(x, 7)*y + 0.0245688 *y*dx*lens_ipow(lambda, 8) + -0.035008 *y*lens_ipow(dx, 5)*lens_ipow(dy, 4) + -0.00375033 *y*lens_ipow(dx, 9) + 0.000523625 *lens_ipow(y, 2)*dx*lens_ipow(dy, 7) + -8.11677e-09 *lens_ipow(x, 2)*lens_ipow(y, 4)*dx*lens_ipow(dy, 3) + -4.33541e-06 *lens_ipow(x, 4)*dx*dy*lens_ipow(lambda, 4)+0.0f;
const float dx11 =  + 0.495423  + 0.027811 *y*dy + 0.000216138 *lens_ipow(y, 2) + 0.00775022 *x*dx + 3.99703e-05 *lens_ipow(x, 2) + -0.0243531 *lens_ipow(dx, 4) + 0.00084875 *lens_ipow(y, 2)*lens_ipow(dy, 2) + -0.00974809 *x*dx*lens_ipow(dy, 2) + 0.0793546 *lens_ipow(lambda, 5) + -8.69587e-07 *lens_ipow(y, 2)*lens_ipow(dy, 3) + 0.000407622 *lens_ipow(y, 2)*lens_ipow(dx, 2)*lens_ipow(lambda, 2) + -1.38802e-05 *lens_ipow(x, 3)*dx*lens_ipow(lambda, 2) + -6.55101e-12 *lens_ipow(y, 8) + -2.7948e-12 *lens_ipow(x, 2)*lens_ipow(y, 6) + -1.68177e-09 *lens_ipow(x, 4)*lens_ipow(y, 2)*lens_ipow(dx, 2) + -6.16018e-13 *lens_ipow(x, 8) + 0.367653 *lens_ipow(dx, 8)*lens_ipow(dy, 2) + 0.0245688 *x*dx*lens_ipow(lambda, 8) + -0.035008 *x*lens_ipow(dx, 5)*lens_ipow(dy, 4) + -0.00375033 *x*lens_ipow(dx, 9) + 0.00104725 *x*y*dx*lens_ipow(dy, 7) + -1.08224e-08 *lens_ipow(x, 3)*lens_ipow(y, 3)*dx*lens_ipow(dy, 3)+0.0f;
const float dx12 =  + -23.5608 *dx*dy + 0.00775022 *x*y + -0.0974125 *y*lens_ipow(dx, 3) + -0.00974809 *x*y*lens_ipow(dy, 2) + 33.1835 *dx*lens_ipow(dy, 3)*lambda + -7.52106 *dx*dy*lens_ipow(lambda, 4) + 10.5599 *dx*lens_ipow(dy, 5) + 27.6121 *lens_ipow(dx, 5)*dy + 0.000271748 *lens_ipow(y, 3)*dx*lens_ipow(lambda, 2) + -1.38802e-05 *lens_ipow(x, 3)*y*lens_ipow(lambda, 2) + -1.96268e-05 *lens_ipow(x, 4)*dx*dy + 65.4936 *lens_ipow(dx, 5)*lens_ipow(dy, 3) + -1.12118e-09 *lens_ipow(x, 4)*lens_ipow(y, 3)*dx + -69.1008 *dx*lens_ipow(dy, 5)*lens_ipow(lambda, 4) + 2.94122 *y*lens_ipow(dx, 7)*lens_ipow(dy, 2) + 0.0245688 *x*y*lens_ipow(lambda, 8) + -0.17504 *x*y*lens_ipow(dx, 4)*lens_ipow(dy, 4) + -0.0337529 *x*y*lens_ipow(dx, 8) + 0.000523625 *x*lens_ipow(y, 2)*lens_ipow(dy, 7) + -2.70559e-09 *lens_ipow(x, 3)*lens_ipow(y, 4)*lens_ipow(dy, 3) + -8.67083e-07 *lens_ipow(x, 5)*dy*lens_ipow(lambda, 4)+0.0f;
const float dx13 =  + 30.4512  + -40.3168 *lens_ipow(dy, 2) + -11.7804 *lens_ipow(dx, 2) + 0.0139055 *lens_ipow(y, 2) + 0.00390391 *lens_ipow(x, 2) + 2.65911 *lens_ipow(lambda, 4) + 26.2899 *lens_ipow(dy, 4) + 0.000565833 *lens_ipow(y, 3)*dy + -0.0194962 *x*y*dx*dy + 49.7753 *lens_ipow(dx, 2)*lens_ipow(dy, 2)*lambda + -8.69587e-07 *lens_ipow(y, 3)*lens_ipow(dy, 2) + -3.76053 *lens_ipow(dx, 2)*lens_ipow(lambda, 4) + 26.3998 *lens_ipow(dx, 2)*lens_ipow(dy, 4) + 4.60202 *lens_ipow(dx, 6) + 0.00341544 *lens_ipow(x, 2)*lens_ipow(dy, 4) + -9.81338e-06 *lens_ipow(x, 4)*lens_ipow(dx, 2) + -2.52316 *lens_ipow(dy, 8) + 32.7468 *lens_ipow(dx, 6)*lens_ipow(dy, 2) + -5.10244 *lens_ipow(lambda, 10) + -172.752 *lens_ipow(dx, 2)*lens_ipow(dy, 4)*lens_ipow(lambda, 4) + 0.735305 *y*lens_ipow(dx, 8)*dy + -0.140032 *x*y*lens_ipow(dx, 5)*lens_ipow(dy, 3) + 0.00366537 *x*lens_ipow(y, 2)*dx*lens_ipow(dy, 6) + -8.11677e-09 *lens_ipow(x, 3)*lens_ipow(y, 4)*dx*lens_ipow(dy, 2) + -8.67083e-07 *lens_ipow(x, 5)*dx*lens_ipow(lambda, 4)+0.0f;
const float dx14 =  + 10.6364 *dy*lens_ipow(lambda, 3) + 16.5918 *lens_ipow(dx, 2)*lens_ipow(dy, 3) + 0.396773 *y*lens_ipow(lambda, 4) + -15.0421 *lens_ipow(dx, 2)*dy*lens_ipow(lambda, 3) + 0.000271748 *lens_ipow(y, 3)*lens_ipow(dx, 2)*lambda + -2.77605e-05 *lens_ipow(x, 3)*y*dx*lambda + -51.0244 *dy*lens_ipow(lambda, 9) + -138.202 *lens_ipow(dx, 2)*lens_ipow(dy, 5)*lens_ipow(lambda, 3) + 0.196551 *x*y*dx*lens_ipow(lambda, 7) + -3.46833e-06 *lens_ipow(x, 5)*dx*dy*lens_ipow(lambda, 3)+0.0f;
const float dx20 =  + -0.0492502  + -0.0215981 *lens_ipow(dy, 2) + -0.0411757 *lens_ipow(dx, 2) + -0.00103873 *y*dy + -4.47595e-06 *lens_ipow(y, 2) + -0.00203225 *x*dx + -1.68619e-06 *lens_ipow(x, 2) + 0.0061882 *lens_ipow(lambda, 4) + 0.000480974 *y*lens_ipow(dx, 2)*dy + -0.000537056 *x*dx*lens_ipow(dy, 2) + 2.66538e-07 *lens_ipow(y, 3)*lens_ipow(dy, 3) + 2.47423e-05 *x*lens_ipow(dx, 5) + 1.96129e-08 *x*lens_ipow(y, 3)*dx*dy + -4.8192e-08 *lens_ipow(x, 4)*lens_ipow(dy, 2) + -8.9827e-14 *lens_ipow(x, 2)*lens_ipow(y, 6) + -1.43176e-06 *lens_ipow(x, 3)*dx*lens_ipow(dy, 4) + -1.47181e-10 *lens_ipow(x, 4)*lens_ipow(y, 2)*lens_ipow(dx, 2) + -0.0286649 *lens_ipow(lambda, 10) + 0.00390668 *lens_ipow(dy, 10) + -4.09758e-08 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(lambda, 6) + -1.28498e-08 *lens_ipow(x, 4)*y*lens_ipow(dx, 4)*dy+0.0f;
const float dx21 =  + 1.19235e-06  + -0.0468609 *dx*dy + -0.00113546 *y*dx + -0.00103873 *x*dy + -8.95191e-06 *x*y + 7.99839e-07 *y*lens_ipow(dx, 2) + -0.000580469 *y*lens_ipow(dx, 3) + 0.000480974 *x*lens_ipow(dx, 2)*dy + 7.99613e-07 *x*lens_ipow(y, 2)*lens_ipow(dy, 3) + 2.94194e-08 *lens_ipow(x, 2)*lens_ipow(y, 2)*dx*dy + 3.96842e-10 *lens_ipow(y, 6)*dx*dy + -1.79654e-13 *lens_ipow(x, 3)*lens_ipow(y, 5) + -5.88723e-11 *lens_ipow(x, 5)*y*lens_ipow(dx, 2) + 5.11651e-05 *lens_ipow(y, 2)*dx*dy*lens_ipow(lambda, 5) + -0.004286 *dx*lens_ipow(dy, 9) + 2.51777e-11 *lens_ipow(y, 7)*dx*lens_ipow(lambda, 2) + -2.73172e-08 *lens_ipow(x, 3)*y*lens_ipow(lambda, 6) + -2.56995e-09 *lens_ipow(x, 5)*lens_ipow(dx, 4)*dy+0.0f;
const float dx22 =  + -1.02651  + -0.902276 *lens_ipow(dy, 2) + -1.65977 *lens_ipow(dx, 2) + -0.0468609 *y*dy + -0.000567729 *lens_ipow(y, 2) + -0.0823514 *x*dx + -0.00101613 *lens_ipow(x, 2) + 0.113201 *lens_ipow(lambda, 3) + 7.99839e-07 *lens_ipow(y, 2)*dx + -0.000870704 *lens_ipow(y, 2)*lens_ipow(dx, 2) + 0.000961949 *x*y*dx*dy + -0.000268528 *lens_ipow(x, 2)*lens_ipow(dy, 2) + 0.315461 *lens_ipow(dx, 4)*lens_ipow(dy, 2) + 6.18557e-05 *lens_ipow(x, 2)*lens_ipow(dx, 4) + 9.80647e-09 *lens_ipow(x, 2)*lens_ipow(y, 3)*dy + 5.66917e-11 *lens_ipow(y, 7)*dy + -3.5794e-07 *lens_ipow(x, 4)*lens_ipow(dy, 4) + -5.88723e-11 *lens_ipow(x, 5)*lens_ipow(y, 2)*dx + 1.7055e-05 *lens_ipow(y, 3)*dy*lens_ipow(lambda, 5) + -0.573908 *lens_ipow(lambda, 10) + 0.0828084 *lens_ipow(dy, 2)*lens_ipow(lambda, 8) + -0.0628251 *lens_ipow(dx, 2)*lens_ipow(lambda, 8) + -0.004286 *y*lens_ipow(dy, 9) + 3.14721e-12 *lens_ipow(y, 8)*lens_ipow(lambda, 2) + -1.02798e-08 *lens_ipow(x, 5)*y*lens_ipow(dx, 3)*dy+0.0f;
const float dx23 =  + 2.87302e-05  + -1.80455 *dx*dy + -0.0468609 *y*dx + -0.0431961 *x*dy + -0.00103873 *x*y + 0.000480974 *x*y*lens_ipow(dx, 2) + -0.000537056 *lens_ipow(x, 2)*dx*dy + 0.126185 *lens_ipow(dx, 5)*dy + 7.99613e-07 *x*lens_ipow(y, 3)*lens_ipow(dy, 2) + 9.80647e-09 *lens_ipow(x, 2)*lens_ipow(y, 3)*dx + -1.92768e-08 *lens_ipow(x, 5)*dy + 5.66917e-11 *lens_ipow(y, 7)*dx + -1.43176e-06 *lens_ipow(x, 4)*dx*lens_ipow(dy, 3) + 1.7055e-05 *lens_ipow(y, 3)*dx*lens_ipow(lambda, 5) + 0.165617 *dx*dy*lens_ipow(lambda, 8) + -0.038574 *y*dx*lens_ipow(dy, 8) + 0.0390668 *x*lens_ipow(dy, 9) + -2.56995e-09 *lens_ipow(x, 5)*y*lens_ipow(dx, 4)+0.0f;
const float dx24 =  + 0.339603 *dx*lens_ipow(lambda, 2) + 0.0247528 *x*lens_ipow(lambda, 3) + 8.52752e-05 *lens_ipow(y, 3)*dx*dy*lens_ipow(lambda, 4) + -5.73908 *dx*lens_ipow(lambda, 9) + 0.662468 *dx*lens_ipow(dy, 2)*lens_ipow(lambda, 7) + -0.167534 *lens_ipow(dx, 3)*lens_ipow(lambda, 7) + 6.29442e-12 *lens_ipow(y, 8)*dx*lambda + -0.286649 *x*lens_ipow(lambda, 9) + -8.19516e-08 *lens_ipow(x, 3)*lens_ipow(y, 2)*lens_ipow(lambda, 5)+0.0f;
const float dx30 =  + -0.000276425 *y*dx + -4.58158e-05 *x*dy + 9.13239e-06 *x*y + -0.0150238 *dx*lens_ipow(dy, 3) + -0.000637033 *x*lens_ipow(dy, 3) + 0.000325284 *x*lens_ipow(dx, 2)*dy + 2.71891e-07 *x*lens_ipow(y, 2)*dy + -0.00327368 *lens_ipow(dx, 5)*dy + -8.80251e-07 *lens_ipow(x, 3)*lens_ipow(dy, 3) + -3.45194e-10 *lens_ipow(x, 5)*dy + -0.000188792 *x*lens_ipow(dy, 7) + 4.9347e-12 *x*lens_ipow(y, 5)*lens_ipow(dx, 2) + -2.73561e-08 *x*lens_ipow(y, 3)*lens_ipow(lambda, 6) + -2.64994e-05 *lens_ipow(x, 2)*dx*lens_ipow(dy, 7) + -1.19652e-13 *lens_ipow(x, 3)*lens_ipow(y, 5)*lens_ipow(dy, 2) + -2.66441e-10 *lens_ipow(x, 5)*y*lens_ipow(dx, 4) + -2.35921e-11 *lens_ipow(x, 7)*lens_ipow(dx, 2)*dy+0.0f;
const float dx31 =  + -0.0492014  + -0.0421999 *lens_ipow(dy, 2) + -3.10808e-06 *dx*dy + -0.00967712 *lens_ipow(dx, 2) + -0.00211281 *y*dy + -3.15439e-06 *lens_ipow(y, 2) + -0.000276425 *x*dx + 4.5662e-06 *lens_ipow(x, 2) + -3.02066e-08 *lens_ipow(y, 2)*dy + 0.00626454 *lens_ipow(lambda, 4) + 2.71891e-07 *lens_ipow(x, 2)*y*dy + 0.00695589 *lens_ipow(dx, 2)*lens_ipow(dy, 4) + -1.91092e-06 *lens_ipow(y, 2)*lens_ipow(dy, 4) + 0.00323231 *lens_ipow(dx, 8) + -8.44922e-10 *lens_ipow(y, 5)*lens_ipow(dx, 2)*dy + 1.23368e-11 *lens_ipow(x, 2)*lens_ipow(y, 4)*lens_ipow(dx, 2) + -0.0293551 *lens_ipow(lambda, 10) + -0.00020412 *y*lens_ipow(dx, 8)*dy + -4.10342e-08 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(lambda, 6) + -1.49564e-13 *lens_ipow(x, 4)*lens_ipow(y, 4)*lens_ipow(dy, 2) + -4.44069e-11 *lens_ipow(x, 6)*lens_ipow(dx, 4)+0.0f;
const float dx32 =  + -0.407311 *dx*dy + -3.10808e-06 *y*dy + -0.0193542 *y*dx + -0.000276425 *x*y + -0.0150238 *x*lens_ipow(dy, 3) + 0.000325284 *lens_ipow(x, 2)*dx*dy + -0.6219 *lens_ipow(dx, 3)*lens_ipow(dy, 3) + 0.0139118 *y*dx*lens_ipow(dy, 4) + -0.0163684 *x*lens_ipow(dx, 4)*dy + 0.0258585 *y*lens_ipow(dx, 7) + -2.81641e-10 *lens_ipow(y, 6)*dx*dy + 4.9347e-12 *lens_ipow(x, 2)*lens_ipow(y, 5)*dx + -0.176341 *dx*lens_ipow(dy, 3)*lens_ipow(lambda, 6) + -0.000816481 *lens_ipow(y, 2)*lens_ipow(dx, 7)*dy + -8.83314e-06 *lens_ipow(x, 3)*lens_ipow(dy, 7) + -1.77627e-10 *lens_ipow(x, 6)*y*lens_ipow(dx, 3) + -5.89802e-12 *lens_ipow(x, 8)*dx*dy+0.0f;
const float dx33 =  + -1.02313  + -1.71059 *lens_ipow(dy, 2) + -0.203656 *lens_ipow(dx, 2) + -0.0843999 *y*dy + -3.10808e-06 *y*dx + -0.0010564 *lens_ipow(y, 2) + -2.29079e-05 *lens_ipow(x, 2) + 0.115364 *lens_ipow(lambda, 3) + -1.00689e-08 *lens_ipow(y, 3) + -0.0450714 *x*dx*lens_ipow(dy, 2) + -0.00095555 *lens_ipow(x, 2)*lens_ipow(dy, 2) + 0.000162642 *lens_ipow(x, 2)*lens_ipow(dx, 2) + 1.35945e-07 *lens_ipow(x, 2)*lens_ipow(y, 2) + -0.466425 *lens_ipow(dx, 4)*lens_ipow(dy, 2) + 0.0278236 *y*lens_ipow(dx, 2)*lens_ipow(dy, 3) + -2.54789e-06 *lens_ipow(y, 3)*lens_ipow(dy, 3) + -0.00327368 *x*lens_ipow(dx, 5) + -6.60189e-07 *lens_ipow(x, 4)*lens_ipow(dy, 2) + -5.75323e-11 *lens_ipow(x, 6) + 0.124802 *lens_ipow(dy, 8) + -1.4082e-10 *lens_ipow(y, 6)*lens_ipow(dx, 2) + -0.000660773 *lens_ipow(x, 2)*lens_ipow(dy, 6) + -0.594001 *lens_ipow(lambda, 10) + -0.264511 *lens_ipow(dx, 2)*lens_ipow(dy, 2)*lens_ipow(lambda, 6) + -0.00010206 *lens_ipow(y, 2)*lens_ipow(dx, 8) + -6.1832e-05 *lens_ipow(x, 3)*dx*lens_ipow(dy, 6) + -5.98258e-14 *lens_ipow(x, 4)*lens_ipow(y, 5)*dy + -2.94901e-12 *lens_ipow(x, 8)*lens_ipow(dx, 2)+0.0f;
const float dx34 =  + 0.346093 *dy*lens_ipow(lambda, 2) + 0.0250582 *y*lens_ipow(lambda, 3) + -5.94001 *dy*lens_ipow(lambda, 9) + -0.529022 *lens_ipow(dx, 2)*lens_ipow(dy, 3)*lens_ipow(lambda, 5) + -0.293551 *y*lens_ipow(lambda, 9) + -8.20684e-08 *lens_ipow(x, 2)*lens_ipow(y, 3)*lens_ipow(lambda, 5)+0.0f;
const float dx40 =  + -4.54686e-07  + -9.82938e-07 *dy + -0.00153956 *dx + -5.71243e-05 *x + 3.87742e-05 *y*dx*dy + 2.65231e-05 *x*lens_ipow(dy, 2) + -4.59541e-06 *lens_ipow(x, 2)*dx*lambda + -7.31443e-10 *lens_ipow(x, 3)*lens_ipow(y, 2) + -1.1046e-07 *x*lens_ipow(y, 2)*lens_ipow(lambda, 5) + -0.00306392 *lens_ipow(dx, 5)*lens_ipow(dy, 4) + 0.000331899 *x*lens_ipow(dx, 8) + 2.73107e-08 *x*lens_ipow(y, 3)*lens_ipow(dy, 5) + 3.26415e-16 *x*lens_ipow(y, 8) + 1.36549e-10 *lens_ipow(x, 3)*lens_ipow(y, 3)*lens_ipow(dx, 2)*dy + 1.74216e-11 *lens_ipow(x, 5)*lens_ipow(dy, 4) + -5.69493e-12 *lens_ipow(x, 7)*lens_ipow(dx, 2)+0.0f;
const float dx41 =  + -0.00152821 *dy + 1.15717e-06 *dx + -5.23832e-05 *y + 3.87742e-05 *x*dx*dy + 4.87492e-05 *y*lens_ipow(dy, 2)*lambda + 5.27922e-05 *y*lens_ipow(dx, 2)*lambda + -3.65721e-10 *lens_ipow(x, 4)*y + -5.24724e-07 *lens_ipow(y, 3)*lens_ipow(dx, 4) + -1.1046e-07 *lens_ipow(x, 2)*y*lens_ipow(lambda, 5) + -0.00319302 *lens_ipow(dy, 9) + -0.00631391 *lens_ipow(dx, 2)*lens_ipow(dy, 7) + -0.00179212 *lens_ipow(dx, 8)*dy + -5.04147e-15 *lens_ipow(y, 9) + 4.09661e-08 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(dy, 5) + 1.30566e-15 *lens_ipow(x, 2)*lens_ipow(y, 7) + 1.02412e-10 *lens_ipow(x, 4)*lens_ipow(y, 2)*lens_ipow(dx, 2)*dy+0.0f;
const float dx42 =  + -0.0736607 *dx + 1.15717e-06 *y + -0.00153956 *x + -0.202951 *dx*lens_ipow(dy, 2) + 3.87742e-05 *x*y*dy + 5.27922e-05 *lens_ipow(y, 2)*dx*lambda + -1.5318e-06 *lens_ipow(x, 3)*lambda + -0.152275 *dx*lens_ipow(dy, 4) + -0.177931 *lens_ipow(dx, 3)*lens_ipow(dy, 2) + -0.554743 *lens_ipow(dx, 5) + -0.0267286 *dx*lens_ipow(lambda, 5) + -5.24724e-07 *lens_ipow(y, 4)*lens_ipow(dx, 3) + -0.0126278 *y*dx*lens_ipow(dy, 7) + -0.0143369 *y*lens_ipow(dx, 7)*dy + -0.0153196 *x*lens_ipow(dx, 4)*lens_ipow(dy, 4) + 0.0013276 *lens_ipow(x, 2)*lens_ipow(dx, 7) + 6.82747e-11 *lens_ipow(x, 4)*lens_ipow(y, 3)*dx*dy + -1.42373e-12 *lens_ipow(x, 8)*dx+0.0f;
const float dx43 =  + -9.98441e-06  + -0.0736102 *dy + -0.00152821 *y + -9.82938e-07 *x + -0.202951 *lens_ipow(dx, 2)*dy + 3.87742e-05 *x*y*dx + 2.65231e-05 *lens_ipow(x, 2)*dy + 4.87492e-05 *lens_ipow(y, 2)*dy*lambda + -0.580729 *lens_ipow(dy, 5) + -0.30455 *lens_ipow(dx, 2)*lens_ipow(dy, 3) + -0.0889656 *lens_ipow(dx, 4)*dy + -0.0191465 *dy*lens_ipow(lambda, 5) + -0.0287372 *y*lens_ipow(dy, 8) + -0.0441974 *y*lens_ipow(dx, 2)*lens_ipow(dy, 6) + -0.00179212 *y*lens_ipow(dx, 8) + -0.0122557 *x*lens_ipow(dx, 5)*lens_ipow(dy, 3) + 6.82768e-08 *lens_ipow(x, 2)*lens_ipow(y, 3)*lens_ipow(dy, 4) + 3.41374e-11 *lens_ipow(x, 4)*lens_ipow(y, 3)*lens_ipow(dx, 2) + 1.16144e-11 *lens_ipow(x, 6)*lens_ipow(dy, 3)+0.0f;
const float dx44 =  + 0.0952603  + -0.239561 *lens_ipow(lambda, 2) + 2.43746e-05 *lens_ipow(y, 2)*lens_ipow(dy, 2) + 2.63961e-05 *lens_ipow(y, 2)*lens_ipow(dx, 2) + -1.5318e-06 *lens_ipow(x, 3)*dx + -0.0478662 *lens_ipow(dy, 2)*lens_ipow(lambda, 4) + -0.0668214 *lens_ipow(dx, 2)*lens_ipow(lambda, 4) + -2.76151e-07 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(lambda, 4) + 1.75289 *lens_ipow(lambda, 10)+0.0f;
