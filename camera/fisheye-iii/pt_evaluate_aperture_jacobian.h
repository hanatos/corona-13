const float dx00 =  + 0.116556  + -1.12507e-05 *x + -0.0167978 *lens_ipow(dy, 2) + 0.889901 *lens_ipow(dx, 2) + 0.0112782 *y*dy + 0.0343099 *x*dx + -5.31713e-06 *lens_ipow(y, 2)*dx + -0.234674 *lens_ipow(lambda, 4) + 2.9518e-07 *lens_ipow(y, 4) + 0.00159433 *x*y*dx*dy + 6.15312e-05 *x*lens_ipow(y, 2)*dx + 2.08476e-06 *lens_ipow(x, 2)*lens_ipow(y, 2) + 1.48949e-06 *lens_ipow(x, 4) + -0.0827476 *y*lens_ipow(dy, 5) + 0.0060907 *lens_ipow(x, 2)*lens_ipow(dx, 2)*lens_ipow(dy, 2) + 2.62347e-07 *lens_ipow(x, 2)*lens_ipow(y, 3)*dy + -0.00495765 *y*dy*lens_ipow(lambda, 5) + -0.000472481 *x*lens_ipow(y, 2)*dx*lens_ipow(dy, 2)*lambda + -0.000465728 *lens_ipow(x, 3)*lens_ipow(dx, 3)*lambda + 7.6351e-11 *lens_ipow(y, 7)*dy + 7.07283e-10 *lens_ipow(x, 6)*y*dy + 1.13258e-09 *lens_ipow(x, 7)*dx + 1.14029 *lens_ipow(lambda, 10) + -1.84095e-08 *lens_ipow(x, 4)*lens_ipow(y, 2)*lens_ipow(lambda, 4)+0.0f;
const float dx01 =  + -1.91794e-05 *y + 0.838948 *dx*dy + 0.0122765 *y*dx + 0.0112782 *x*dy + -1.06343e-05 *x*y*dx + -0.0320092 *y*lens_ipow(dx, 3) + 1.18072e-06 *x*lens_ipow(y, 3) + 0.000797167 *lens_ipow(x, 2)*dx*dy + 6.15312e-05 *lens_ipow(x, 2)*y*dx + 1.38984e-06 *lens_ipow(x, 3)*y + 4.37404e-06 *lens_ipow(y, 4)*dx*dy + 9.41882e-08 *lens_ipow(y, 5)*dx + -0.0827476 *x*lens_ipow(dy, 5) + 2.62347e-07 *lens_ipow(x, 3)*lens_ipow(y, 2)*dy + -0.00495765 *x*dy*lens_ipow(lambda, 5) + -0.000472481 *lens_ipow(x, 2)*y*dx*lens_ipow(dy, 2)*lambda + 5.34457e-10 *x*lens_ipow(y, 6)*dy + 1.0104e-10 *lens_ipow(x, 7)*dy + -7.36378e-09 *lens_ipow(x, 5)*y*lens_ipow(lambda, 4)+0.0f;
const float dx02 =  + 37.6368  + -0.153206 *lens_ipow(dy, 2) + 3.51684 *lens_ipow(dx, 2) + 0.838948 *y*dy + 0.00613825 *lens_ipow(y, 2) + 1.7798 *x*dx + 0.0171549 *lens_ipow(x, 2) + -5.31713e-06 *x*lens_ipow(y, 2) + -7.4088 *lens_ipow(lambda, 4) + -0.0480137 *lens_ipow(y, 2)*lens_ipow(dx, 2) + 0.000797167 *lens_ipow(x, 2)*y*dy + 3.07656e-05 *lens_ipow(x, 2)*lens_ipow(y, 2) + 1931.5 *lens_ipow(dx, 4)*lens_ipow(dy, 2) + 8.74808e-07 *lens_ipow(y, 5)*dy + 1.5698e-08 *lens_ipow(y, 6) + 0.00406046 *lens_ipow(x, 3)*dx*lens_ipow(dy, 2) + -0.00023624 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(dy, 2)*lambda + -0.000349296 *lens_ipow(x, 4)*lens_ipow(dx, 2)*lambda + 538.864 *lens_ipow(dy, 8) + 7131.71 *lens_ipow(dx, 2)*lens_ipow(dy, 6) + 3843.51 *lens_ipow(dx, 8) + 1.41572e-10 *lens_ipow(x, 8) + 35.8328 *lens_ipow(lambda, 10)+0.0f;
const float dx03 =  + -0.306411 *dx*dy + 0.838948 *y*dx + -0.0335956 *x*dy + 0.0112782 *x*y + 0.000797167 *lens_ipow(x, 2)*y*dx + 772.598 *lens_ipow(dx, 5)*dy + 8.74808e-07 *lens_ipow(y, 5)*dx + -0.413738 *x*y*lens_ipow(dy, 4) + 0.00406046 *lens_ipow(x, 3)*lens_ipow(dx, 2)*dy + 8.74489e-08 *lens_ipow(x, 3)*lens_ipow(y, 3) + -0.00495765 *x*y*lens_ipow(lambda, 5) + -0.000472481 *lens_ipow(x, 2)*lens_ipow(y, 2)*dx*dy*lambda + 4310.91 *dx*lens_ipow(dy, 7) + 14263.4 *lens_ipow(dx, 3)*lens_ipow(dy, 5) + 7.6351e-11 *x*lens_ipow(y, 7) + 1.0104e-10 *lens_ipow(x, 7)*y+0.0f;
const float dx04 =  + -29.6352 *dx*lens_ipow(lambda, 3) + -0.938696 *x*lens_ipow(lambda, 3) + -0.0247882 *x*y*dy*lens_ipow(lambda, 4) + -0.00023624 *lens_ipow(x, 2)*lens_ipow(y, 2)*dx*lens_ipow(dy, 2) + -0.000116432 *lens_ipow(x, 4)*lens_ipow(dx, 3) + 358.328 *dx*lens_ipow(lambda, 9) + 11.4029 *x*lens_ipow(lambda, 9) + -1.47276e-08 *lens_ipow(x, 5)*lens_ipow(y, 2)*lens_ipow(lambda, 3)+0.0f;
const float dx10 =  + 0.918219 *dx*dy + 0.0108813 *y*dx + 0.0130531 *x*dy + 0.000200154 *x*y + 4.95863e-05 *x*lens_ipow(y, 2)*dy + 4.71671e-05 *lens_ipow(x, 2)*y*dx + 0.00387435 *lens_ipow(y, 2)*dx*lens_ipow(dy, 3) + 0.00554742 *lens_ipow(y, 2)*lens_ipow(dx, 3)*dy + -0.000303021 *x*y*lens_ipow(lambda, 4) + 0.00427207 *x*y*lens_ipow(dx, 4) + 2.81116e-09 *x*lens_ipow(y, 5) + -0.000157263 *lens_ipow(x, 3)*lens_ipow(dx, 2)*dy + 3.27332e-09 *lens_ipow(x, 5)*y + -0.474885 *dx*dy*lens_ipow(lambda, 5) + -0.0152121 *y*dx*lens_ipow(lambda, 5) + -0.0965704 *x*lens_ipow(dy, 3)*lens_ipow(lambda, 3) + 1.25035e-07 *lens_ipow(y, 5)*dx*lens_ipow(lambda, 2) + 2.14763e-05 *lens_ipow(x, 3)*y*lens_ipow(dy, 4) + 2.9301e-10 *lens_ipow(x, 7)*dy+0.0f;
const float dx11 =  + 0.114645  + 0.885817 *lens_ipow(dy, 2) + 0.0125457 *lens_ipow(dx, 2) + 0.0328423 *y*dy + 0.000272429 *lens_ipow(y, 2) + 0.0108813 *x*dx + 0.000100077 *lens_ipow(x, 2) + -0.152055 *lens_ipow(lambda, 3) + -0.000844463 *lens_ipow(y, 2)*lens_ipow(dx, 2) + 4.95863e-05 *lens_ipow(x, 2)*y*dy + 1.57224e-05 *lens_ipow(x, 3)*dx + -0.167334 *y*lens_ipow(dx, 2)*dy*lens_ipow(lambda, 2) + 3.19007e-07 *lens_ipow(y, 5)*dy + 0.0077487 *x*y*dx*lens_ipow(dy, 3) + 0.0110948 *x*y*lens_ipow(dx, 3)*dy + -0.00015151 *lens_ipow(x, 2)*lens_ipow(lambda, 4) + 0.00213603 *lens_ipow(x, 2)*lens_ipow(dx, 4) + 7.0279e-09 *lens_ipow(x, 2)*lens_ipow(y, 4) + 5.45554e-10 *lens_ipow(x, 6) + -0.0152121 *x*dx*lens_ipow(lambda, 5) + 5.70141e-05 *lens_ipow(y, 4)*lens_ipow(dy, 4) + 6.25177e-07 *x*lens_ipow(y, 4)*dx*lens_ipow(lambda, 2) + 5.36908e-06 *lens_ipow(x, 4)*lens_ipow(dy, 4) + 0.810168 *lens_ipow(lambda, 10) + 3.23138e-14 *lens_ipow(y, 10)+0.0f;
const float dx12 =  + 3.44368 *dx*dy + 0.0250915 *y*dx + 0.918219 *x*dy + 0.0108813 *x*y + -0.000562975 *lens_ipow(y, 3)*dx + 1.57224e-05 *lens_ipow(x, 3)*y + -0.167334 *lens_ipow(y, 2)*dx*dy*lens_ipow(lambda, 2) + 0.00387435 *x*lens_ipow(y, 2)*lens_ipow(dy, 3) + 0.0166423 *x*lens_ipow(y, 2)*lens_ipow(dx, 2)*dy + 0.00854413 *lens_ipow(x, 2)*y*lens_ipow(dx, 3) + -7.86316e-05 *lens_ipow(x, 4)*dx*dy + -0.474885 *x*dy*lens_ipow(lambda, 5) + -0.0152121 *x*y*lens_ipow(lambda, 5) + 32052.7 *lens_ipow(dx, 3)*lens_ipow(dy, 5) + 6293.93 *lens_ipow(dx, 7)*dy + 1.25035e-07 *x*lens_ipow(y, 5)*lens_ipow(lambda, 2)+0.0f;
const float dx13 =  + 37.6385  + 2.31613 *lens_ipow(dy, 2) + 1.72184 *lens_ipow(dx, 2) + 1.77163 *y*dy + 0.0164211 *lens_ipow(y, 2) + 0.918219 *x*dx + 0.00652653 *lens_ipow(x, 2) + -7.44186 *lens_ipow(lambda, 4) + 2.47932e-05 *lens_ipow(x, 2)*lens_ipow(y, 2) + -0.0836668 *lens_ipow(y, 2)*lens_ipow(dx, 2)*lens_ipow(lambda, 2) + 5.31678e-08 *lens_ipow(y, 6) + 0.011623 *x*lens_ipow(y, 2)*dx*lens_ipow(dy, 2) + 0.00554742 *x*lens_ipow(y, 2)*lens_ipow(dx, 3) + -3.93158e-05 *lens_ipow(x, 4)*lens_ipow(dx, 2) + -0.474885 *x*dx*lens_ipow(lambda, 5) + -0.144856 *lens_ipow(x, 2)*lens_ipow(dy, 2)*lens_ipow(lambda, 3) + 5313.89 *lens_ipow(dy, 8) + 40065.9 *lens_ipow(dx, 4)*lens_ipow(dy, 4) + 786.742 *lens_ipow(dx, 8) + 4.56113e-05 *lens_ipow(y, 5)*lens_ipow(dy, 3) + 2.14763e-05 *lens_ipow(x, 4)*y*lens_ipow(dy, 3) + 3.66262e-11 *lens_ipow(x, 8) + 35.8743 *lens_ipow(lambda, 10)+0.0f;
const float dx14 =  + -0.456165 *y*lens_ipow(lambda, 2) + -29.7675 *dy*lens_ipow(lambda, 3) + -0.167334 *lens_ipow(y, 2)*lens_ipow(dx, 2)*dy*lambda + -0.000606041 *lens_ipow(x, 2)*y*lens_ipow(lambda, 3) + -2.37442 *x*dx*dy*lens_ipow(lambda, 4) + -0.0760605 *x*y*dx*lens_ipow(lambda, 4) + -0.144856 *lens_ipow(x, 2)*lens_ipow(dy, 3)*lens_ipow(lambda, 2) + 2.50071e-07 *x*lens_ipow(y, 5)*dx*lambda + 358.743 *dy*lens_ipow(lambda, 9) + 8.10168 *y*lens_ipow(lambda, 9)+0.0f;
const float dx20 =  + -0.0285989  + 0.000398916 *lens_ipow(dy, 2) + -0.0109949 *lens_ipow(dx, 2) + -0.000213514 *y*dy + -9.9664e-06 *lens_ipow(y, 2) + -0.00139911 *x*dx + -2.98883e-05 *lens_ipow(x, 2) + 0.00112127 *y*lens_ipow(dx, 2)*dy + 2.38682e-05 *lens_ipow(x, 2)*lens_ipow(dy, 2) + -0.000508211 *lens_ipow(lambda, 5) + -0.0870952 *lens_ipow(dx, 2)*lens_ipow(lambda, 4) + 5.53917e-05 *lens_ipow(y, 2)*lens_ipow(dx, 4) + -3.15492e-06 *lens_ipow(y, 3)*lens_ipow(dy, 3) + -6.08308e-12 *lens_ipow(y, 6) + 2.01544e-07 *x*lens_ipow(y, 3)*dx*dy + 0.000141905 *lens_ipow(x, 2)*lens_ipow(dx, 2)*lens_ipow(dy, 2) + 2.70122e-06 *lens_ipow(x, 3)*y*dx*lens_ipow(dy, 3) + 7.76681e-11 *lens_ipow(x, 6)*lens_ipow(dx, 2) + -4.20676e-13 *lens_ipow(x, 6)*lens_ipow(y, 2) + -0.195445 *lens_ipow(dy, 4)*lens_ipow(lambda, 5) + -2.69019e-10 *lens_ipow(x, 2)*lens_ipow(y, 4)*lens_ipow(lambda, 3) + -2.67153e-08 *lens_ipow(x, 3)*lens_ipow(y, 2)*dx*lens_ipow(lambda, 3) + -0.00247517 *y*dy*lens_ipow(lambda, 8) + 2.40364e-05 *lens_ipow(y, 3)*lens_ipow(dx, 2)*dy*lens_ipow(lambda, 4) + -3.11635e-05 *lens_ipow(x, 3)*dx*lens_ipow(lambda, 6) + -4.71401e-10 *lens_ipow(x, 6)*lens_ipow(lambda, 4)+0.0f;
const float dx21 =  + -0.0169563 *dx*dy + -0.0010155 *y*dx + -0.000213514 *x*dy + -1.99328e-05 *x*y + 0.00112127 *x*lens_ipow(dx, 2)*dy + 0.00888972 *y*dx*lens_ipow(dy, 4) + 0.000110783 *x*y*lens_ipow(dx, 4) + -9.46476e-06 *x*lens_ipow(y, 2)*lens_ipow(dy, 3) + -3.64985e-11 *x*lens_ipow(y, 5) + 3.02316e-07 *lens_ipow(x, 2)*lens_ipow(y, 2)*dx*dy + -0.0566703 *dx*dy*lens_ipow(lambda, 5) + 6.75305e-07 *lens_ipow(x, 4)*dx*lens_ipow(dy, 3) + -1.20193e-13 *lens_ipow(x, 7)*y + -4.00872e-05 *lens_ipow(y, 3)*lens_ipow(dx, 3)*lens_ipow(lambda, 3) + -3.58692e-10 *lens_ipow(x, 3)*lens_ipow(y, 3)*lens_ipow(lambda, 3) + -1.33577e-08 *lens_ipow(x, 4)*y*dx*lens_ipow(lambda, 3) + -0.00247517 *x*dy*lens_ipow(lambda, 8) + 7.21093e-05 *x*lens_ipow(y, 2)*lens_ipow(dx, 2)*dy*lens_ipow(lambda, 4)+0.0f;
const float dx22 =  + -0.667314  + -1.44854 *lens_ipow(dy, 2) + -3.30876 *lens_ipow(dx, 2) + -0.0169563 *y*dy + -0.00050775 *lens_ipow(y, 2) + -0.0219899 *x*dx + -0.000699554 *lens_ipow(x, 2) + -0.236851 *lens_ipow(lambda, 3) + -3.91141 *lens_ipow(dx, 2)*lens_ipow(lambda, 2) + 0.00224255 *x*y*dx*dy + 0.00444486 *lens_ipow(y, 2)*lens_ipow(dy, 4) + -0.17419 *x*dx*lens_ipow(lambda, 4) + 0.000221567 *x*lens_ipow(y, 2)*lens_ipow(dx, 3) + 1.00772e-07 *lens_ipow(x, 2)*lens_ipow(y, 3)*dy + 9.46032e-05 *lens_ipow(x, 3)*dx*lens_ipow(dy, 2) + -0.0566703 *y*dy*lens_ipow(lambda, 5) + 6.75305e-07 *lens_ipow(x, 4)*y*lens_ipow(dy, 3) + 2.21909e-11 *lens_ipow(x, 7)*dx + -16.5736 *lens_ipow(dy, 4)*lens_ipow(lambda, 5) + -3.00654e-05 *lens_ipow(y, 4)*lens_ipow(dx, 2)*lens_ipow(lambda, 3) + -6.67883e-09 *lens_ipow(x, 4)*lens_ipow(y, 2)*lens_ipow(lambda, 3) + 1.11761 *lens_ipow(lambda, 10) + -267.502 *lens_ipow(dx, 4)*lens_ipow(dy, 2)*lens_ipow(lambda, 4) + 4.80729e-05 *x*lens_ipow(y, 3)*dx*dy*lens_ipow(lambda, 4) + -7.79088e-06 *lens_ipow(x, 4)*lens_ipow(lambda, 6)+0.0f;
const float dx23 =  + -2.89709 *dx*dy + -0.0169563 *y*dx + 0.000797833 *x*dy + -0.000213514 *x*y + 0.00112127 *x*y*lens_ipow(dx, 2) + 1.59121e-05 *lens_ipow(x, 3)*dy + 0.0177794 *lens_ipow(y, 2)*dx*lens_ipow(dy, 3) + -9.46476e-06 *x*lens_ipow(y, 3)*lens_ipow(dy, 2) + 1.00772e-07 *lens_ipow(x, 2)*lens_ipow(y, 3)*dx + 9.46032e-05 *lens_ipow(x, 3)*lens_ipow(dx, 2)*dy + -0.0566703 *y*dx*lens_ipow(lambda, 5) + 2.02591e-06 *lens_ipow(x, 4)*y*dx*lens_ipow(dy, 2) + -66.2946 *dx*lens_ipow(dy, 3)*lens_ipow(lambda, 5) + -0.781779 *x*lens_ipow(dy, 3)*lens_ipow(lambda, 5) + -107.001 *lens_ipow(dx, 5)*dy*lens_ipow(lambda, 4) + -0.00247517 *x*y*lens_ipow(lambda, 8) + 2.40364e-05 *x*lens_ipow(y, 3)*lens_ipow(dx, 2)*lens_ipow(lambda, 4)+0.0f;
const float dx24 =  + -0.710552 *dx*lens_ipow(lambda, 2) + -2.60761 *lens_ipow(dx, 3)*lambda + -0.00254105 *x*lens_ipow(lambda, 4) + -0.348381 *x*lens_ipow(dx, 2)*lens_ipow(lambda, 3) + -0.283352 *y*dx*dy*lens_ipow(lambda, 4) + -82.8682 *dx*lens_ipow(dy, 4)*lens_ipow(lambda, 4) + -3.00654e-05 *lens_ipow(y, 4)*lens_ipow(dx, 3)*lens_ipow(lambda, 2) + -0.977223 *x*lens_ipow(dy, 4)*lens_ipow(lambda, 4) + -2.69019e-10 *lens_ipow(x, 3)*lens_ipow(y, 4)*lens_ipow(lambda, 2) + -2.00365e-08 *lens_ipow(x, 4)*lens_ipow(y, 2)*dx*lens_ipow(lambda, 2) + 11.1761 *dx*lens_ipow(lambda, 9) + -214.002 *lens_ipow(dx, 5)*lens_ipow(dy, 2)*lens_ipow(lambda, 3) + -0.0198014 *x*y*dy*lens_ipow(lambda, 7) + 9.61458e-05 *x*lens_ipow(y, 3)*lens_ipow(dx, 2)*dy*lens_ipow(lambda, 3) + -4.67453e-05 *lens_ipow(x, 4)*dx*lens_ipow(lambda, 5) + -2.69372e-10 *lens_ipow(x, 7)*lens_ipow(lambda, 3)+0.0f;
const float dx30 =  + -0.0183763 *dx*dy + -0.000251102 *y*dx + -0.00101541 *x*dy + -1.97462e-05 *x*y + -1.07008e-05 *y*lens_ipow(dx, 2) + 0.000927268 *y*dx*lens_ipow(dy, 2) + 0.000164545 *x*y*lens_ipow(dx, 4) + -1.26629e-06 *x*lens_ipow(y, 2)*lens_ipow(dy, 3) + -3.85225e-06 *x*lens_ipow(y, 2)*lens_ipow(dx, 2)*dy + -6.41862e-11 *x*lens_ipow(y, 5) + -0.000206813 *lens_ipow(x, 2)*lens_ipow(dx, 3)*dy + -2.3606e-09 *lens_ipow(x, 2)*lens_ipow(y, 3)*dx + -5.07909e-11 *lens_ipow(x, 5)*y + -0.0313921 *dx*dy*lens_ipow(lambda, 5) + -1.37552e-05 *lens_ipow(x, 3)*lens_ipow(dy, 3)*lambda + -1.24907e-08 *x*lens_ipow(y, 4)*dy*lens_ipow(lambda, 2) + -8.01994e-05 *x*y*lens_ipow(dy, 2)*lens_ipow(lambda, 5) + -1.37717e-06 *lens_ipow(y, 3)*dx*lens_ipow(lambda, 6)+0.0f;
const float dx31 =  + -0.0285768  + -3.91201e-05 *dy + -1.96027e-05 *dx + -0.0108042 *lens_ipow(dy, 2) + -0.000295841 *lens_ipow(dx, 2) + -0.0013562 *y*dy + -2.91995e-05 *lens_ipow(y, 2) + -0.000251102 *x*dx + -9.87308e-06 *lens_ipow(x, 2) + -1.07008e-05 *x*lens_ipow(dx, 2) + -0.000665779 *lens_ipow(lambda, 4) + 0.000927268 *x*dx*lens_ipow(dy, 2) + -0.0902573 *lens_ipow(dy, 2)*lens_ipow(lambda, 4) + -1.55445e-05 *lens_ipow(y, 3)*lens_ipow(dx, 2)*dy + 8.22724e-05 *lens_ipow(x, 2)*lens_ipow(dx, 4) + -1.26629e-06 *lens_ipow(x, 2)*y*lens_ipow(dy, 3) + -3.85225e-06 *lens_ipow(x, 2)*y*lens_ipow(dx, 2)*dy + -1.60466e-10 *lens_ipow(x, 2)*lens_ipow(y, 4) + -2.3606e-09 *lens_ipow(x, 3)*lens_ipow(y, 2)*dx + -8.46515e-12 *lens_ipow(x, 6) + -2.49814e-08 *lens_ipow(x, 2)*lens_ipow(y, 3)*dy*lens_ipow(lambda, 2) + -4.00997e-05 *lens_ipow(x, 2)*lens_ipow(dy, 2)*lens_ipow(lambda, 5) + -0.00923224 *y*dy*lens_ipow(lambda, 8) + -1.2379e-07 *lens_ipow(y, 5)*lens_ipow(dy, 5) + -3.20643e-10 *lens_ipow(y, 6)*lens_ipow(lambda, 4) + -4.1315e-06 *x*lens_ipow(y, 2)*dx*lens_ipow(lambda, 6)+0.0f;
const float dx32 =  + -1.96027e-05 *y + -2.94907 *dx*dy + -0.000591682 *y*dx + -0.0183763 *x*dy + -0.000251102 *x*y + -2.14016e-05 *x*y*dx + 0.000927268 *x*y*lens_ipow(dy, 2) + -7.77227e-06 *lens_ipow(y, 4)*dx*dy + 0.00032909 *lens_ipow(x, 2)*y*lens_ipow(dx, 3) + -3.85225e-06 *lens_ipow(x, 2)*lens_ipow(y, 2)*dx*dy + -0.000206813 *lens_ipow(x, 3)*lens_ipow(dx, 2)*dy + -7.86867e-10 *lens_ipow(x, 3)*lens_ipow(y, 3) + -0.0313921 *x*dy*lens_ipow(lambda, 5) + -30.8725 *lens_ipow(dx, 3)*dy*lens_ipow(lambda, 5) + -1.37717e-06 *x*lens_ipow(y, 3)*lens_ipow(lambda, 6)+0.0f;
const float dx33 =  + -0.667848  + -0.003186 *dy + -3.91201e-05 *y + -3.25945 *lens_ipow(dy, 2) + -1.47453 *lens_ipow(dx, 2) + -0.0216084 *y*dy + -0.0006781 *lens_ipow(y, 2) + -0.0183763 *x*dx + -0.000507707 *lens_ipow(x, 2) + -0.240532 *lens_ipow(lambda, 3) + -4.16633 *lens_ipow(dy, 2)*lens_ipow(lambda, 2) + 0.00185454 *x*y*dx*dy + -0.180515 *y*dy*lens_ipow(lambda, 4) + -3.88614e-06 *lens_ipow(y, 4)*lens_ipow(dx, 2) + -1.89943e-06 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(dy, 2) + -1.92613e-06 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(dx, 2) + -6.89377e-05 *lens_ipow(x, 3)*lens_ipow(dx, 3) + -0.0313921 *x*dx*lens_ipow(lambda, 5) + -1.03164e-05 *lens_ipow(x, 4)*lens_ipow(dy, 2)*lambda + -6.24534e-09 *lens_ipow(x, 2)*lens_ipow(y, 4)*lens_ipow(lambda, 2) + -7.71813 *lens_ipow(dx, 4)*lens_ipow(lambda, 5) + -8.01994e-05 *lens_ipow(x, 2)*y*dy*lens_ipow(lambda, 5) + 1.50561 *lens_ipow(lambda, 10) + -0.00461612 *lens_ipow(y, 2)*lens_ipow(lambda, 8) + -1.03158e-07 *lens_ipow(y, 6)*lens_ipow(dy, 4)+0.0f;
const float dx34 =  + -0.721595 *dy*lens_ipow(lambda, 2) + -2.77755 *lens_ipow(dy, 3)*lambda + -0.00266311 *y*lens_ipow(lambda, 3) + -0.361029 *y*lens_ipow(dy, 2)*lens_ipow(lambda, 3) + -0.156961 *x*dx*dy*lens_ipow(lambda, 4) + -3.43879e-06 *lens_ipow(x, 4)*lens_ipow(dy, 3) + -1.24907e-08 *lens_ipow(x, 2)*lens_ipow(y, 4)*dy*lambda + -38.5907 *lens_ipow(dx, 4)*dy*lens_ipow(lambda, 4) + -0.000200498 *lens_ipow(x, 2)*y*lens_ipow(dy, 2)*lens_ipow(lambda, 4) + 15.0561 *dy*lens_ipow(lambda, 9) + -0.0369289 *lens_ipow(y, 2)*dy*lens_ipow(lambda, 7) + -1.83224e-10 *lens_ipow(y, 7)*lens_ipow(lambda, 3) + -8.263e-06 *x*lens_ipow(y, 3)*dx*lens_ipow(lambda, 5)+0.0f;
const float dx40 =  + -0.000396702 *dx + -6.70034e-06 *x + -0.00215847 *dx*lens_ipow(dy, 2) + -3.39926e-05 *y*dx*dy + -1.93753e-06 *lens_ipow(y, 2)*dx + -5.40813e-05 *x*lens_ipow(dy, 2) + 4.31424e-05 *x*lens_ipow(dx, 2) + -3.75778e-06 *x*y*dy + -1.30183e-07 *x*lens_ipow(y, 2) + 2.52368e-08 *x*lens_ipow(y, 2)*lens_ipow(dy, 2) + -1.77254e-08 *lens_ipow(x, 4)*dx + -4.43911e-10 *lens_ipow(x, 5) + 0.00726618 *y*dx*lens_ipow(dy, 5) + 0.0146502 *y*lens_ipow(dx, 3)*lens_ipow(dy, 3) + 0.00888838 *y*lens_ipow(dx, 5)*dy + 8.19731e-08 *lens_ipow(y, 4)*dx*lens_ipow(dy, 2) + 2.4331e-07 *lens_ipow(x, 3)*y*lens_ipow(dx, 2)*dy + -0.000251479 *y*dx*dy*lens_ipow(lambda, 5)+0.0f;
const float dx41 =  + 5.47844e-07  + -0.000376043 *dy + -7.06492e-07 *dx + -6.90047e-06 *y + -1.05147e-05 *y*lens_ipow(dx, 2) + -3.56498e-06 *lens_ipow(y, 2)*dy + -3.39926e-05 *x*dx*dy + -3.87506e-06 *x*y*dx + -1.87889e-06 *lens_ipow(x, 2)*dy + -1.30183e-07 *lens_ipow(x, 2)*y + -1.38401e-07 *lens_ipow(y, 3)*lens_ipow(dy, 2) + -4.47473e-10 *lens_ipow(y, 5) + 2.52368e-08 *lens_ipow(x, 2)*y*lens_ipow(dy, 2) + -0.00385523 *y*lens_ipow(dx, 6) + 0.00726618 *x*dx*lens_ipow(dy, 5) + 0.0146502 *x*lens_ipow(dx, 3)*lens_ipow(dy, 3) + 0.00888838 *x*lens_ipow(dx, 5)*dy + 3.27892e-07 *x*lens_ipow(y, 3)*dx*lens_ipow(dy, 2) + 6.08274e-08 *lens_ipow(x, 4)*lens_ipow(dx, 2)*dy + -0.000251479 *x*dx*dy*lens_ipow(lambda, 5) + -0.000857763 *lens_ipow(y, 2)*lens_ipow(dy, 7) + -9.14938e-10 *lens_ipow(y, 6)*lens_ipow(dx, 2)*dy+0.0f;
const float dx42 =  + -0.0554564 *dx + -7.06492e-07 *y + -0.000396702 *x + -1.05147e-05 *lens_ipow(y, 2)*dx + -0.00215847 *x*lens_ipow(dy, 2) + -3.39926e-05 *x*y*dy + -1.93753e-06 *x*lens_ipow(y, 2) + 4.31424e-05 *lens_ipow(x, 2)*dx + -9.58965 *dx*lens_ipow(dy, 4) + -20.3855 *lens_ipow(dx, 3)*lens_ipow(dy, 2) + -7.65898 *lens_ipow(dx, 5) + -3.54509e-09 *lens_ipow(x, 5) + -0.0115657 *lens_ipow(y, 2)*lens_ipow(dx, 5) + 0.00726618 *x*y*lens_ipow(dy, 5) + 0.0439505 *x*y*lens_ipow(dx, 2)*lens_ipow(dy, 3) + 0.0444419 *x*y*lens_ipow(dx, 4)*dy + 8.19731e-08 *x*lens_ipow(y, 4)*lens_ipow(dy, 2) + 1.21655e-07 *lens_ipow(x, 4)*y*dx*dy + -0.000251479 *x*y*dy*lens_ipow(lambda, 5) + -2.61411e-10 *lens_ipow(y, 7)*dx*dy+0.0f;
const float dx43 =  + 2.19793e-05  + -0.0482325 *dy + -0.000376043 *y + -1.18833e-06 *lens_ipow(y, 3) + -0.00431695 *x*dx*dy + -3.39926e-05 *x*y*dx + -5.40813e-05 *lens_ipow(x, 2)*dy + -1.87889e-06 *lens_ipow(x, 2)*y + -9.80974 *lens_ipow(dy, 5) + -19.1793 *lens_ipow(dx, 2)*lens_ipow(dy, 3) + -10.1928 *lens_ipow(dx, 4)*dy + -6.92004e-08 *lens_ipow(y, 4)*dy + 2.52368e-08 *lens_ipow(x, 2)*lens_ipow(y, 2)*dy + 0.0363309 *x*y*dx*lens_ipow(dy, 4) + 0.0439505 *x*y*lens_ipow(dx, 3)*lens_ipow(dy, 2) + 0.00888838 *x*y*lens_ipow(dx, 5) + 1.63946e-07 *x*lens_ipow(y, 4)*dx*dy + 6.08274e-08 *lens_ipow(x, 4)*y*lens_ipow(dx, 2) + -0.000251479 *x*y*dx*lens_ipow(lambda, 5) + -0.00200145 *lens_ipow(y, 3)*lens_ipow(dy, 6) + -1.30705e-10 *lens_ipow(y, 7)*lens_ipow(dx, 2)+0.0f;
const float dx44 =  + 0.27353  + -0.686048 *lens_ipow(lambda, 2) + -0.0012574 *x*y*dx*dy*lens_ipow(lambda, 4) + 4.56155 *lens_ipow(lambda, 10)+0.0f;