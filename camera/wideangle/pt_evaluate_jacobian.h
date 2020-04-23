const float dx00 =  + -1.82627  + 0.859651 *lens_ipow(dy, 2) + 2.4542 *lens_ipow(dx, 2) + 0.0583735 *y*dy + 0.000158755 *y*dx + 0.0031524 *lens_ipow(y, 2) + 0.22372 *x*dx + 0.00960063 *lens_ipow(x, 2) + -0.287397 *lens_ipow(lambda, 4) + 1.42363 *lens_ipow(dx, 4) + -0.315245 *y*lens_ipow(dx, 2)*dy + 0.00797403 *lens_ipow(y, 2)*lens_ipow(dx, 2) + -0.020144 *x*y*dx*dy + 0.0120823 *lens_ipow(x, 2)*lens_ipow(dy, 2) + -0.899481 *y*lens_ipow(dy, 5) + -28.9797 *lens_ipow(dx, 2)*lens_ipow(dy, 2)*lens_ipow(lambda, 3) + 1752.5 *lens_ipow(dx, 2)*lens_ipow(dy, 6) + 0.196608 *x*lens_ipow(y, 2)*dx*lens_ipow(dy, 4) + -1.7311e-09 *lens_ipow(x, 2)*lens_ipow(y, 6) + -1.12339e-08 *lens_ipow(x, 5)*lens_ipow(y, 2)*dx + -2.71833e-09 *lens_ipow(x, 6)*lens_ipow(y, 2) + -6.56453e-10 *lens_ipow(x, 8) + 1.24533 *lens_ipow(lambda, 10) + 3.51606e-05 *lens_ipow(x, 3)*lens_ipow(y, 3)*dx*dy*lens_ipow(lambda, 2) + 4.16768e-06 *lens_ipow(x, 6)*lens_ipow(dx, 2)*lens_ipow(lambda, 2)+0.0f;
const float dx01 =  + 0.000854024 *dx + 1.58758 *dx*dy + 0.104797 *y*dx + 4.32597e-06 *lens_ipow(y, 2) + 0.0583735 *x*dy + 0.000158755 *x*dx + 0.00630479 *x*y + 0.0018585 *y*lens_ipow(dx, 2) + -0.315245 *x*lens_ipow(dx, 2)*dy + 0.0159481 *x*y*lens_ipow(dx, 2) + -0.010072 *lens_ipow(x, 2)*dx*dy + 4.41184 *y*lens_ipow(dx, 5) + 0.00540816 *lens_ipow(y, 3)*lens_ipow(dx, 3) + 9.68594e-05 *lens_ipow(y, 4)*dx*dy + -0.899481 *x*lens_ipow(dy, 5) + 0.196608 *lens_ipow(x, 2)*y*dx*lens_ipow(dy, 4) + -3.4622e-09 *lens_ipow(x, 3)*lens_ipow(y, 5) + -3.74462e-09 *lens_ipow(x, 6)*y*dx + -7.76664e-10 *lens_ipow(x, 7)*y + 2.63705e-05 *lens_ipow(x, 4)*lens_ipow(y, 2)*dx*dy*lens_ipow(lambda, 2)+0.0f;
const float dx02 =  + 15.983  + 0.000854024 *y + 28.0752 *lens_ipow(dy, 2) + 82.8344 *lens_ipow(dx, 2) + 1.58758 *y*dy + 0.0523987 *lens_ipow(y, 2) + 4.9084 *x*dx + 0.000158755 *x*y + 0.11186 *lens_ipow(x, 2) + 0.0018585 *lens_ipow(y, 2)*dx + 5.6945 *x*lens_ipow(dx, 3) + -0.63049 *x*y*dx*dy + 0.0159481 *x*lens_ipow(y, 2)*dx + -0.010072 *lens_ipow(x, 2)*y*dy + 297.43 *lens_ipow(dx, 4)*lens_ipow(lambda, 2) + 11.0296 *lens_ipow(y, 2)*lens_ipow(dx, 4) + 0.00405612 *lens_ipow(y, 4)*lens_ipow(dx, 2) + 1.93719e-05 *lens_ipow(y, 5)*dy + -57.9595 *x*dx*lens_ipow(dy, 2)*lens_ipow(lambda, 3) + -4.85736 *lens_ipow(lambda, 8) + 3505.01 *x*dx*lens_ipow(dy, 6) + 0.0983038 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(dy, 4) + -1.87231e-09 *lens_ipow(x, 6)*lens_ipow(y, 2) + -5.44756e+06 *lens_ipow(dx, 8)*lens_ipow(dy, 2) + 8.79015e-06 *lens_ipow(x, 4)*lens_ipow(y, 3)*dy*lens_ipow(lambda, 2) + 1.19077e-06 *lens_ipow(x, 7)*dx*lens_ipow(lambda, 2)+0.0f;
const float dx03 =  + 56.1503 *dx*dy + 1.58758 *y*dx + 1.7193 *x*dy + 0.0583735 *x*y + -0.315245 *x*y*lens_ipow(dx, 2) + -0.010072 *lens_ipow(x, 2)*y*dx + 0.00805487 *lens_ipow(x, 3)*dy + 1.93719e-05 *lens_ipow(y, 5)*dx + -4.49741 *x*y*lens_ipow(dy, 4) + -57.9595 *x*lens_ipow(dx, 2)*dy*lens_ipow(lambda, 3) + 10515 *x*lens_ipow(dx, 2)*lens_ipow(dy, 5) + 0.393215 *lens_ipow(x, 2)*lens_ipow(y, 2)*dx*lens_ipow(dy, 3) + -1.21057e+06 *lens_ipow(dx, 9)*dy + 8.79015e-06 *lens_ipow(x, 4)*lens_ipow(y, 3)*dx*lens_ipow(lambda, 2)+0.0f;
const float dx04 =  + -1.14959 *x*lens_ipow(lambda, 3) + 118.972 *lens_ipow(dx, 5)*lambda + -86.9392 *x*lens_ipow(dx, 2)*lens_ipow(dy, 2)*lens_ipow(lambda, 2) + -38.8589 *dx*lens_ipow(lambda, 7) + 12.4533 *x*lens_ipow(lambda, 9) + 1.75803e-05 *lens_ipow(x, 4)*lens_ipow(y, 3)*dx*dy*lambda + 1.19077e-06 *lens_ipow(x, 7)*lens_ipow(dx, 2)*lambda+0.0f;
const float dx10 =  + 4.4789e-05  + 1.6031 *dx*dy + 0.0599259 *y*dx + 0.102875 *x*dy + 0.00651622 *x*y + -0.146538 *y*dx*lens_ipow(dy, 2) + -0.00723343 *lens_ipow(y, 2)*dx*dy + 0.0115327 *x*y*lens_ipow(dy, 2) + 10.7047 *x*lens_ipow(dy, 5) + -8.00964e-08 *x*lens_ipow(y, 5) + -0.0141743 *lens_ipow(x, 2)*y*dx*lens_ipow(dy, 2) + -1.12059e-07 *lens_ipow(x, 3)*lens_ipow(y, 3) + -8.69128e-08 *lens_ipow(x, 5)*y + -0.0291065 *lens_ipow(y, 3)*lens_ipow(dx, 5) + 0.131818 *x*lens_ipow(y, 2)*lens_ipow(dx, 4)*dy + -1.9663e-09 *lens_ipow(x, 3)*lens_ipow(y, 4)*dy + 1.80589e-06 *lens_ipow(y, 6)*dx*dy*lens_ipow(lambda, 2) + 55.2566 *x*y*lens_ipow(dy, 8)+0.0f;
const float dx11 =  + -1.82256  + 2.54778 *lens_ipow(dy, 2) + 0.940597 *lens_ipow(dx, 2) + 0.226771 *y*dy + 0.00968289 *lens_ipow(y, 2) + 0.0599259 *x*dx + 0.00325811 *lens_ipow(x, 2) + -0.186441 *lens_ipow(lambda, 3) + -1.33388 *lens_ipow(dy, 4) + 0.00781577 *lens_ipow(y, 2)*lens_ipow(dx, 2) + -0.146538 *x*dx*lens_ipow(dy, 2) + -0.0144669 *x*y*dx*dy + 0.00576635 *lens_ipow(x, 2)*lens_ipow(dy, 2) + 0.0600501 *y*lens_ipow(dy, 4) + 0.000147127 *lens_ipow(y, 4)*lens_ipow(dy, 2) + -2.00241e-07 *lens_ipow(x, 2)*lens_ipow(y, 4) + -0.00472476 *lens_ipow(x, 3)*dx*lens_ipow(dy, 2) + -8.40445e-08 *lens_ipow(x, 4)*lens_ipow(y, 2) + -1.44855e-08 *lens_ipow(x, 6) + -27.0863 *lens_ipow(dx, 2)*lens_ipow(dy, 2)*lens_ipow(lambda, 4) + -7.63851e-10 *lens_ipow(y, 8) + -0.0873195 *x*lens_ipow(y, 2)*lens_ipow(dx, 5) + 0.131818 *lens_ipow(x, 2)*y*lens_ipow(dx, 4)*dy + -1.9663e-09 *lens_ipow(x, 4)*lens_ipow(y, 3)*dy + 0.880054 *lens_ipow(lambda, 10) + 2.3037 *lens_ipow(y, 3)*lens_ipow(dx, 2)*lens_ipow(dy, 5) + 1.08354e-05 *x*lens_ipow(y, 5)*dx*dy*lens_ipow(lambda, 2) + 27.6283 *lens_ipow(x, 2)*lens_ipow(dy, 8)+0.0f;
const float dx12 =  + 0.0451305 *dx + 56.2831 *dx*dy + 1.88119 *y*dx + 1.6031 *x*dy + 0.0599259 *x*y + 0.00521051 *lens_ipow(y, 3)*dx + -0.146538 *x*y*lens_ipow(dy, 2) + -0.00723343 *x*lens_ipow(y, 2)*dy + -0.00472476 *lens_ipow(x, 3)*y*lens_ipow(dy, 2) + -54.1725 *y*dx*lens_ipow(dy, 2)*lens_ipow(lambda, 4) + -0.145532 *x*lens_ipow(y, 3)*lens_ipow(dx, 4) + 0.263635 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(dx, 3)*dy + -1.75491e+06 *dx*lens_ipow(dy, 9) + 1.15185 *lens_ipow(y, 4)*dx*lens_ipow(dy, 5) + 1.80589e-06 *x*lens_ipow(y, 6)*dy*lens_ipow(lambda, 2)+0.0f;
const float dx13 =  + 16.0201  + -0.124323 *dy + 78.6313 *lens_ipow(dy, 2) + 28.1415 *lens_ipow(dx, 2) + 5.09555 *y*dy + 0.113385 *lens_ipow(y, 2) + 1.6031 *x*dx + 0.0514373 *lens_ipow(x, 2) + -5.33551 *y*lens_ipow(dy, 3) + -0.293075 *x*y*dx*dy + -0.00723343 *x*lens_ipow(y, 2)*dx + 0.0115327 *lens_ipow(x, 2)*y*dy + 0.1201 *lens_ipow(y, 2)*lens_ipow(dy, 3) + 5.88507e-05 *lens_ipow(y, 5)*dy + 26.7617 *lens_ipow(x, 2)*lens_ipow(dy, 4) + -0.00944951 *lens_ipow(x, 3)*y*dx*dy + -3.23319 *lens_ipow(lambda, 7) + -54.1725 *y*lens_ipow(dx, 2)*dy*lens_ipow(lambda, 4) + 0.0659088 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(dx, 4) + -4.91574e-10 *lens_ipow(x, 4)*lens_ipow(y, 4) + -7.8971e+06 *lens_ipow(dx, 2)*lens_ipow(dy, 8) + 2.87962 *lens_ipow(y, 4)*lens_ipow(dx, 2)*lens_ipow(dy, 4) + 1.80589e-06 *x*lens_ipow(y, 6)*dx*lens_ipow(lambda, 2) + 221.026 *lens_ipow(x, 2)*y*lens_ipow(dy, 7)+0.0f;
const float dx14 =  + -0.559322 *y*lens_ipow(lambda, 2) + -22.6323 *dy*lens_ipow(lambda, 6) + -108.345 *y*lens_ipow(dx, 2)*lens_ipow(dy, 2)*lens_ipow(lambda, 3) + 8.80054 *y*lens_ipow(lambda, 9) + 3.61179e-06 *x*lens_ipow(y, 6)*dx*dy*lambda+0.0f;
const float dx20 =  + -0.0234805  + 1.05778e-07 *y + 0.00768814 *lens_ipow(dx, 2) + -0.000225225 *y*dy + -9.86015e-06 *lens_ipow(y, 2) + -0.000649462 *x*dx + 3.54495e-05 *lens_ipow(x, 2) + 0.00505966 *lens_ipow(lambda, 3) + 0.00338762 *lens_ipow(dy, 2)*lambda + -0.00312758 *y*lens_ipow(dx, 2)*dy + -0.000204932 *x*y*dx*dy + -6.38312e-06 *lens_ipow(x, 2)*y*dy + -0.000529579 *lens_ipow(x, 2)*lens_ipow(dx, 4) + 3.46763e-06 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(dx, 2) + 0.00921361 *lens_ipow(x, 2)*lens_ipow(dy, 6) + -1.17033e-08 *lens_ipow(x, 6)*lens_ipow(dy, 2) + -1.89157e-11 *lens_ipow(x, 6)*lens_ipow(y, 2) + -9.29412e-12 *lens_ipow(x, 8) + -2.10319e-07 *lens_ipow(x, 4)*lens_ipow(lambda, 5) + 4.27959e-08 *lens_ipow(x, 6)*lens_ipow(dx, 2)*lambda + -0.0274803 *lens_ipow(lambda, 10) + -0.300797 *lens_ipow(dx, 4)*lens_ipow(lambda, 6) + -1.78428e-05 *lens_ipow(y, 5)*lens_ipow(dx, 4)*dy + 9.31965e-08 *lens_ipow(x, 5)*y*dx*dy*lens_ipow(lambda, 2)+0.0f;
const float dx21 =  + 2.55585e-07 *y + 1.05778e-07 *x + 0.0128846 *dx*dy + 2.73165e-05 *lens_ipow(dx, 2) + 0.000351702 *y*dx + -0.000225225 *x*dy + -1.97203e-05 *x*y + -0.0410086 *lens_ipow(dx, 3)*dy + -0.000159894 *y*dx*lens_ipow(lambda, 2) + -0.00312758 *x*lens_ipow(dx, 2)*dy + -0.000102466 *lens_ipow(x, 2)*dx*dy + -2.12771e-06 *lens_ipow(x, 3)*dy + 0.0362047 *y*lens_ipow(dx, 5) + 2.31175e-06 *lens_ipow(x, 3)*y*lens_ipow(dx, 2) + 1.11412e-06 *lens_ipow(y, 5)*lens_ipow(dx, 3) + -5.40447e-12 *lens_ipow(x, 7)*y + -8.92138e-05 *x*lens_ipow(y, 4)*lens_ipow(dx, 4)*dy + 1.55327e-08 *lens_ipow(x, 6)*dx*dy*lens_ipow(lambda, 2)+0.0f;
const float dx22 =  + -0.339988  + 0.41412 *lens_ipow(dy, 2) + 1.40538 *lens_ipow(dx, 2) + 0.0128846 *y*dy + 5.46329e-05 *y*dx + 0.000175851 *lens_ipow(y, 2) + 0.0153763 *x*dx + -0.000324731 *lens_ipow(x, 2) + -0.0219382 *lens_ipow(lambda, 3) + -0.123026 *y*lens_ipow(dx, 2)*dy + -7.99468e-05 *lens_ipow(y, 2)*lens_ipow(lambda, 2) + -0.00625515 *x*y*dx*dy + -0.000102466 *lens_ipow(x, 2)*y*dy + 0.0905117 *lens_ipow(y, 2)*lens_ipow(dx, 4) + -0.000706105 *lens_ipow(x, 3)*lens_ipow(dx, 3) + 2.31175e-06 *lens_ipow(x, 3)*lens_ipow(y, 2)*dx + -2328.51 *lens_ipow(dx, 6)*lens_ipow(dy, 2) + -965.936 *lens_ipow(dx, 8) + 5.57058e-07 *lens_ipow(y, 6)*lens_ipow(dx, 2) + 1.09275 *lens_ipow(dy, 2)*lens_ipow(lambda, 7) + 1.22274e-08 *lens_ipow(x, 7)*dx*lambda + -1.20319 *x*lens_ipow(dx, 3)*lens_ipow(lambda, 6) + -7.13711e-05 *x*lens_ipow(y, 5)*lens_ipow(dx, 3)*dy + 1.55327e-08 *lens_ipow(x, 6)*y*dy*lens_ipow(lambda, 2)+0.0f;
const float dx23 =  + 0.82824 *dx*dy + 0.0128846 *y*dx + -0.000225225 *x*y + 0.00677524 *x*dy*lambda + -0.0410086 *y*lens_ipow(dx, 3) + -0.00312758 *x*y*lens_ipow(dx, 2) + -0.000102466 *lens_ipow(x, 2)*y*dx + -2.12771e-06 *lens_ipow(x, 3)*y + -665.289 *lens_ipow(dx, 7)*dy + 0.0184272 *lens_ipow(x, 3)*lens_ipow(dy, 5) + -3.34379e-09 *lens_ipow(x, 7)*dy + 2.1855 *dx*dy*lens_ipow(lambda, 7) + -1.78428e-05 *x*lens_ipow(y, 5)*lens_ipow(dx, 4) + 1.55327e-08 *lens_ipow(x, 6)*y*dx*lens_ipow(lambda, 2)+0.0f;
const float dx24 =  + -0.0658147 *dx*lens_ipow(lambda, 2) + 0.015179 *x*lens_ipow(lambda, 2) + 0.00338762 *x*lens_ipow(dy, 2) + -0.000159894 *lens_ipow(y, 2)*dx*lambda + 7.64926 *dx*lens_ipow(dy, 2)*lens_ipow(lambda, 6) + -2.10319e-07 *lens_ipow(x, 5)*lens_ipow(lambda, 4) + 6.11371e-09 *lens_ipow(x, 7)*lens_ipow(dx, 2) + -0.274803 *x*lens_ipow(lambda, 9) + -1.80478 *x*lens_ipow(dx, 4)*lens_ipow(lambda, 5) + 3.10655e-08 *lens_ipow(x, 6)*y*dx*dy*lambda+0.0f;
const float dx30 =  + 0.00698629 *dx*dy + 2.60451e-06 *y*dy + -0.000159778 *y*dx + -0.000771513 *x*dy + 6.54481e-05 *x*y + 3.03038e-05 *y*dx*dy + 2.785e-06 *y*dx*lens_ipow(dy, 2) + -0.00126211 *x*lens_ipow(dx, 2)*dy + -0.00039787 *lens_ipow(y, 2)*lens_ipow(dx, 3)*dy + -0.000370777 *lens_ipow(x, 2)*y*dx*lens_ipow(dy, 2) + 0.000102194 *lens_ipow(x, 3)*lens_ipow(dy, 3) + -2.80415e-09 *lens_ipow(x, 5)*y + 9.16997e-09 *lens_ipow(y, 6)*dx*dy + 3.23717e-05 *x*lens_ipow(y, 3)*lens_ipow(dy, 4) + -5.7409e-11 *lens_ipow(x, 3)*lens_ipow(y, 5) + -0.0338558 *x*y*lens_ipow(dx, 2)*lens_ipow(dy, 2)*lens_ipow(lambda, 4) + 3.69147e-12 *x*lens_ipow(y, 8)*dy+0.0f;
const float dx31 =  + -0.0234193  + 0.00672364 *lens_ipow(dy, 2) + -0.00635551 *lens_ipow(dx, 2) + -0.000760741 *y*dy + 3.13018e-05 *lens_ipow(y, 2) + 2.60451e-06 *x*dy + -0.000159778 *x*dx + 3.2724e-05 *lens_ipow(x, 2) + 0.00504652 *lens_ipow(lambda, 3) + 3.03038e-05 *x*dx*dy + 0.00438687 *y*lens_ipow(dx, 2)*dy + 2.785e-06 *x*dx*lens_ipow(dy, 2) + 4.48499e-05 *lens_ipow(y, 2)*lens_ipow(dy, 4) + 0.00181484 *lens_ipow(y, 2)*lens_ipow(dx, 4) + -0.00079574 *x*y*lens_ipow(dx, 3)*dy + -0.000123592 *lens_ipow(x, 3)*dx*lens_ipow(dy, 2) + -4.67359e-10 *lens_ipow(x, 6) + -8.09074e-12 *lens_ipow(y, 8) + 5.50198e-08 *x*lens_ipow(y, 5)*dx*dy + 4.85576e-05 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(dy, 4) + -7.17613e-11 *lens_ipow(x, 4)*lens_ipow(y, 4) + 2.5712e-08 *lens_ipow(y, 6)*lens_ipow(dy, 2)*lambda + -0.0279742 *lens_ipow(lambda, 10) + -0.187182 *lens_ipow(dy, 4)*lens_ipow(lambda, 6) + -0.0169279 *lens_ipow(x, 2)*lens_ipow(dx, 2)*lens_ipow(dy, 2)*lens_ipow(lambda, 4) + 1.47659e-11 *lens_ipow(x, 2)*lens_ipow(y, 7)*dy+0.0f;
const float dx32 =  + 7.91718e-06  + -0.000271803 *dy + -0.00481802 *lens_ipow(dy, 2) + 0.95314 *dx*dy + -0.012711 *y*dx + 0.00698629 *x*dy + -0.000159778 *x*y + 3.03038e-05 *x*y*dy + 0.00438687 *lens_ipow(y, 2)*dx*dy + 2.785e-06 *x*y*lens_ipow(dy, 2) + -0.00126211 *lens_ipow(x, 2)*dx*dy + 0.00241979 *lens_ipow(y, 3)*lens_ipow(dx, 3) + -0.00119361 *x*lens_ipow(y, 2)*lens_ipow(dx, 2)*dy + -0.000123592 *lens_ipow(x, 3)*y*lens_ipow(dy, 2) + 9.16997e-09 *x*lens_ipow(y, 6)*dy + -76037.3 *lens_ipow(dx, 3)*lens_ipow(dy, 7) + -0.0338558 *lens_ipow(x, 2)*y*dx*lens_ipow(dy, 2)*lens_ipow(lambda, 4)+0.0f;
const float dx33 =  + -0.3394  + -0.000917504 *dy + -0.000271803 *dx + 1.3373 *lens_ipow(dy, 2) + -0.00963604 *dx*dy + 0.47657 *lens_ipow(dx, 2) + 0.0134473 *y*dy + -0.000380371 *lens_ipow(y, 2) + 0.00698629 *x*dx + 2.60451e-06 *x*y + -0.000385757 *lens_ipow(x, 2) + -0.019918 *lens_ipow(lambda, 3) + 3.03038e-05 *x*y*dx + 0.00219344 *lens_ipow(y, 2)*lens_ipow(dx, 2) + 5.56999e-06 *x*y*dx*dy + -0.000631056 *lens_ipow(x, 2)*lens_ipow(dx, 2) + 5.97999e-05 *lens_ipow(y, 3)*lens_ipow(dy, 3) + -0.00039787 *x*lens_ipow(y, 2)*lens_ipow(dx, 3) + -0.000247185 *lens_ipow(x, 3)*y*dx*dy + 7.66458e-05 *lens_ipow(x, 4)*lens_ipow(dy, 2) + 9.16997e-09 *x*lens_ipow(y, 6)*dx + 6.47435e-05 *lens_ipow(x, 2)*lens_ipow(y, 3)*lens_ipow(dy, 3) + 7.34628e-09 *lens_ipow(y, 7)*dy*lambda + -133065 *lens_ipow(dx, 4)*lens_ipow(dy, 6) + -0.748728 *y*lens_ipow(dy, 3)*lens_ipow(lambda, 6) + -0.0338558 *lens_ipow(x, 2)*y*lens_ipow(dx, 2)*dy*lens_ipow(lambda, 4) + 1.84573e-12 *lens_ipow(x, 2)*lens_ipow(y, 8)+0.0f;
const float dx34 =  + -0.0597539 *dy*lens_ipow(lambda, 2) + 0.0151396 *y*lens_ipow(lambda, 2) + 3.67314e-09 *lens_ipow(y, 7)*lens_ipow(dy, 2) + -0.279742 *y*lens_ipow(lambda, 9) + -1.12309 *y*lens_ipow(dy, 4)*lens_ipow(lambda, 5) + -0.0677116 *lens_ipow(x, 2)*y*lens_ipow(dx, 2)*lens_ipow(dy, 2)*lens_ipow(lambda, 3)+0.0f;
const float dx40 =  + -1.40242e-07  + -0.00111004 *dx + -4.75352e-05 *x + -0.0524655 *dx*lens_ipow(dy, 2) + -0.0506835 *lens_ipow(dx, 3) + -0.00136887 *y*dx*dy + -2.33081e-06 *lens_ipow(y, 2)*dx + 0.0609999 *y*dx*lens_ipow(dy, 3) + 0.0857349 *y*lens_ipow(dx, 3)*dy + 3.67271e-08 *x*lens_ipow(y, 3)*dy + -0.000899009 *lens_ipow(x, 2)*lens_ipow(dx, 3) + 3.50856e-09 *lens_ipow(x, 5) + 1.33942e-05 *lens_ipow(y, 4)*dx*lens_ipow(dy, 2) + -0.137456 *x*lens_ipow(dy, 6) + -0.00241226 *x*lens_ipow(y, 2)*lens_ipow(dx, 2)*lens_ipow(dy, 2) + 0.000863948 *lens_ipow(x, 3)*lens_ipow(dx, 2)*lens_ipow(dy, 2) + 5.69851e-05 *lens_ipow(x, 3)*y*lens_ipow(dx, 2)*dy + -9.98436e-08 *lens_ipow(x, 5)*lens_ipow(dy, 2) + -5.27261e-07 *lens_ipow(x, 5)*lens_ipow(dx, 2)+0.0f;
const float dx41 =  + -0.00206827 *dy + -4.23935e-05 *y + 0.00108764 *y*lens_ipow(dy, 2) + -0.00136887 *x*dx*dy + -4.66162e-06 *x*y*dx + -2.45605 *lens_ipow(dx, 2)*lens_ipow(dy, 3) + -0.0197757 *y*lens_ipow(dx, 4) + 5.86348e-07 *lens_ipow(y, 4)*dy + 0.0609999 *x*dx*lens_ipow(dy, 3) + 0.0857349 *x*lens_ipow(dx, 3)*dy + 5.50907e-08 *lens_ipow(x, 2)*lens_ipow(y, 2)*dy + 0.36325 *y*lens_ipow(dy, 6) + 2.44343e-07 *lens_ipow(y, 5)*lens_ipow(dx, 2) + 5.35769e-05 *x*lens_ipow(y, 3)*dx*lens_ipow(dy, 2) + -0.00241226 *lens_ipow(x, 2)*y*lens_ipow(dx, 2)*lens_ipow(dy, 2) + 1.42463e-05 *lens_ipow(x, 4)*lens_ipow(dx, 2)*dy + -0.032145 *lens_ipow(y, 2)*lens_ipow(dx, 6)*dy+0.0f;
const float dx42 =  + -0.146641 *dx + -0.00111004 *x + -0.00279243 *lens_ipow(dy, 3) + -0.0524655 *x*lens_ipow(dy, 2) + -0.152051 *x*lens_ipow(dx, 2) + -0.00136887 *x*y*dy + -2.33081e-06 *x*lens_ipow(y, 2) + -253.741 *dx*lens_ipow(dy, 4) + -391.838 *lens_ipow(dx, 3)*lens_ipow(dy, 2) + -180.158 *lens_ipow(dx, 5) + -4.9121 *y*dx*lens_ipow(dy, 3) + -0.0395514 *lens_ipow(y, 2)*lens_ipow(dx, 3) + 0.0609999 *x*y*lens_ipow(dy, 3) + 0.257205 *x*y*lens_ipow(dx, 2)*dy + -0.000899009 *lens_ipow(x, 3)*lens_ipow(dx, 2) + 8.14476e-08 *lens_ipow(y, 6)*dx + 1.33942e-05 *x*lens_ipow(y, 4)*lens_ipow(dy, 2) + -0.00241226 *lens_ipow(x, 2)*lens_ipow(y, 2)*dx*lens_ipow(dy, 2) + 0.000431974 *lens_ipow(x, 4)*dx*lens_ipow(dy, 2) + 2.84926e-05 *lens_ipow(x, 4)*y*dx*dy + -1.75754e-07 *lens_ipow(x, 6)*dx + -0.0642901 *lens_ipow(y, 3)*lens_ipow(dx, 5)*dy+0.0f;
const float dx43 =  + 3.57585e-05  + -0.137336 *dy + -0.00206827 *y + -0.00837728 *dx*lens_ipow(dy, 2) + 0.00108764 *lens_ipow(y, 2)*dy + -0.104931 *x*dx*dy + -0.00136887 *x*y*dx + -166.477 *lens_ipow(dy, 5) + -507.482 *lens_ipow(dx, 2)*lens_ipow(dy, 3) + -195.919 *lens_ipow(dx, 4)*dy + -7.36815 *y*lens_ipow(dx, 2)*lens_ipow(dy, 2) + 1.1727e-07 *lens_ipow(y, 5) + 0.183 *x*y*dx*lens_ipow(dy, 2) + 0.0857349 *x*y*lens_ipow(dx, 3) + 1.83636e-08 *lens_ipow(x, 2)*lens_ipow(y, 3) + 1.08975 *lens_ipow(y, 2)*lens_ipow(dy, 5) + 2.67884e-05 *x*lens_ipow(y, 4)*dx*dy + -0.412369 *lens_ipow(x, 2)*lens_ipow(dy, 5) + -0.00241226 *lens_ipow(x, 2)*lens_ipow(y, 2)*lens_ipow(dx, 2)*dy + 0.000431974 *lens_ipow(x, 4)*lens_ipow(dx, 2)*dy + 1.42463e-05 *lens_ipow(x, 4)*y*lens_ipow(dx, 2) + -3.32812e-08 *lens_ipow(x, 6)*dy + -0.010715 *lens_ipow(y, 3)*lens_ipow(dx, 6)+0.0f;
const float dx44 =  + 0.248005  + -0.610383 *lens_ipow(lambda, 2) + 3.97792 *lens_ipow(lambda, 10)+0.0f;