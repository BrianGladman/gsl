#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_gbmv () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int trans = 111;
   int M = 3;
   int N = 2;
   int KL = 1;
   int KU = 1;
   int lda = 3;
   float alpha = 1;
   float beta = 0;
   float A[] = { -0.001, 0.601, -0.21, 0.233, 0.879, 0.627, 0.907, 0.944, 0.021 };
   float X[] = { 0.876, -0.599 };
   int incX = 1;
   float Y[] = { 0.387, 0.016, -0.611 };
   int incY = -1;
   float y_expected[] = { -0.565456, 0.268011, -0.155891 };
   cblas_sgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sgbmv(case 540)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int M = 3;
   int N = 2;
   int KL = 1;
   int KU = 1;
   int lda = 3;
   float alpha = 1;
   float beta = 0;
   float A[] = { -0.001, 0.601, -0.21, 0.233, 0.879, 0.627, 0.907, 0.944, 0.021 };
   float X[] = { 0.876, -0.599 };
   int incX = 1;
   float Y[] = { 0.387, 0.016, -0.611 };
   int incY = -1;
   float y_expected[] = { -0.375573, -0.710481, 0.386909 };
   cblas_sgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sgbmv(case 541)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int M = 3;
   int N = 2;
   int KL = 1;
   int KU = 1;
   int lda = 3;
   float alpha = -0.3;
   float beta = 1;
   float A[] = { -0.159, 0.655, 0.338, -0.465, -0.128, 0.126, 0.055, 0.122, 0.77 };
   float X[] = { 0.992, 0.006, -0.247 };
   int incX = 1;
   float Y[] = { 0.68, -0.423 };
   int incY = -1;
   float y_expected[] = { 0.4943426, -0.284715 };
   cblas_sgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sgbmv(case 542)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int M = 3;
   int N = 2;
   int KL = 1;
   int KU = 1;
   int lda = 3;
   float alpha = -0.3;
   float beta = 1;
   float A[] = { -0.159, 0.655, 0.338, -0.465, -0.128, 0.126, 0.055, 0.122, 0.77 };
   float X[] = { 0.992, 0.006, -0.247 };
   int incX = 1;
   float Y[] = { 0.68, -0.423 };
   int incY = -1;
   float y_expected[] = { 0.827951, -0.6185364 };
   cblas_sgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sgbmv(case 543)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int M = 3;
   int N = 2;
   int KL = 1;
   int KU = 1;
   int lda = 3;
   double alpha = -0.3;
   double beta = -1;
   double A[] = { -0.502, -0.773, -0.709, 0.09, 0.261, 0.84, -0.112, -0.319, 0.254 };
   double X[] = { 0.237, -0.034 };
   int incX = 1;
   double Y[] = { -0.56, 0.729, -0.498 };
   int incY = -1;
   double y_expected[] = { 0.5567462, -0.7183746, 0.4837164 };
   cblas_dgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dgbmv(case 544)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int M = 3;
   int N = 2;
   int KL = 1;
   int KU = 1;
   int lda = 3;
   double alpha = -0.3;
   double beta = -1;
   double A[] = { -0.502, -0.773, -0.709, 0.09, 0.261, 0.84, -0.112, -0.319, 0.254 };
   double X[] = { 0.237, -0.034 };
   int incX = 1;
   double Y[] = { -0.56, 0.729, -0.498 };
   int incY = -1;
   double y_expected[] = { 0.568568, -0.6759279, 0.5538783 };
   cblas_dgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dgbmv(case 545)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int M = 3;
   int N = 2;
   int KL = 1;
   int KU = 1;
   int lda = 3;
   double alpha = 0.1;
   double beta = 1;
   double A[] = { 0.238, -0.379, 0.894, 0.282, 0.24, 0.321, -0.952, 0.095, 0.761 };
   double X[] = { 0.76, 0.95, 0.323 };
   int incX = 1;
   double Y[] = { 0.593, -0.763 };
   int incY = -1;
   double y_expected[] = { 0.5900645, -0.832008 };
   cblas_dgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dgbmv(case 546)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int M = 3;
   int N = 2;
   int KL = 1;
   int KU = 1;
   int lda = 3;
   double alpha = 0.1;
   double beta = 1;
   double A[] = { 0.238, -0.379, 0.894, 0.282, 0.24, 0.321, -0.952, 0.095, 0.761 };
   double X[] = { 0.76, 0.95, 0.323 };
   int incX = 1;
   double Y[] = { 0.593, -0.763 };
   int incY = -1;
   double y_expected[] = { 0.6476003, -0.706874 };
   cblas_dgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dgbmv(case 547)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int M = 3;
   int N = 2;
   int KL = 1;
   int KU = 1;
   int lda = 3;
   float alpha[2] = {-1, 0};
   float beta[2] = {0, 0.1};
   float A[] = { 0.754, 0.195, 0.31, -0.826, -0.327, 0.41, 0.961, 0.412, -0.047, 0.585, -0.604, 0.482, -0.892, -0.586, -0.302, 0.965, -0.597, -0.322 };
   float X[] = { 0.242, 0.704, -0.862, -0.441 };
   int incX = 1;
   float Y[] = { 0.466, -0.02, -0.12, 0.282, 0.529, -0.722 };
   int incY = -1;
   float y_expected[] = { -0.683889, 0.745248, -0.523379, 1.241323, 0.761172, -1.29865 };
   cblas_cgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "cgbmv(case 548) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "cgbmv(case 548) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int M = 3;
   int N = 2;
   int KL = 1;
   int KU = 1;
   int lda = 3;
   float alpha[2] = {-1, 0};
   float beta[2] = {0, 0.1};
   float A[] = { 0.754, 0.195, 0.31, -0.826, -0.327, 0.41, 0.961, 0.412, -0.047, 0.585, -0.604, 0.482, -0.892, -0.586, -0.302, 0.965, -0.597, -0.322 };
   float X[] = { 0.242, 0.704, -0.862, -0.441 };
   int incX = 1;
   float Y[] = { 0.466, -0.02, -0.12, 0.282, 0.529, -0.722 };
   int incY = -1;
   float y_expected[] = { -0.73121, 0.19572, 0.041075, 0.602531, 0.062366, 0.813497 };
   cblas_cgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "cgbmv(case 549) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "cgbmv(case 549) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int M = 3;
   int N = 2;
   int KL = 1;
   int KU = 1;
   int lda = 3;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {-0.3, 0.1};
   float A[] = { -0.145, 0.845, -0.106, -0.018, -0.781, 0.342, -0.424, -0.412, 0.012, -0.583, 0.669, -0.529, 0.389, -0.528, -0.406, 0.059, -0.991, -0.945 };
   float X[] = { -0.117, 0.247, 0.293, 0.374, -0.006, 0.483 };
   int incX = 1;
   float Y[] = { 0.984, -0.923, 0.255, -0.269 };
   int incY = -1;
   float y_expected[] = { -0.2279176, 0.5125922, -0.1818721, 0.1722047 };
   cblas_cgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "cgbmv(case 550) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "cgbmv(case 550) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int M = 3;
   int N = 2;
   int KL = 1;
   int KU = 1;
   int lda = 3;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {-0.3, 0.1};
   float A[] = { -0.145, 0.845, -0.106, -0.018, -0.781, 0.342, -0.424, -0.412, 0.012, -0.583, 0.669, -0.529, 0.389, -0.528, -0.406, 0.059, -0.991, -0.945 };
   float X[] = { -0.117, 0.247, 0.293, 0.374, -0.006, 0.483 };
   int incX = 1;
   float Y[] = { 0.984, -0.923, 0.255, -0.269 };
   int incY = -1;
   float y_expected[] = { -0.4005715, 0.4067085, 0.0739643, 0.1369999 };
   cblas_cgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "cgbmv(case 551) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "cgbmv(case 551) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int M = 3;
   int N = 2;
   int KL = 1;
   int KU = 1;
   int lda = 3;
   float alpha[2] = {0, 1};
   float beta[2] = {-0.3, 0.1};
   float A[] = { 0.84, -0.806, -0.711, 0.818, 0.18, 0.237, -0.582, 0.192, -0.744, 0.188, 0.478, 0.464, -0.337, 0.001, 0.953, 0.663, 0.127, 0.534 };
   float X[] = { 0.009, -0.248, 0.686, 0.493, -0.368, 0.1 };
   int incX = 1;
   float Y[] = { 0.619, -0.928, -0.578, 0.48 };
   int incY = -1;
   float y_expected[] = { -0.10539, -0.571067, 0.149619, -0.485343 };
   cblas_cgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "cgbmv(case 552) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "cgbmv(case 552) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int M = 3;
   int N = 2;
   int KL = 1;
   int KU = 1;
   int lda = 3;
   float alpha[2] = {0, 1};
   float beta[2] = {-0.3, 0.1};
   float A[] = { 0.84, -0.806, -0.711, 0.818, 0.18, 0.237, -0.582, 0.192, -0.744, 0.188, 0.478, 0.464, -0.337, 0.001, 0.953, 0.663, 0.127, 0.534 };
   float X[] = { 0.009, -0.248, 0.686, 0.493, -0.368, 0.1 };
   int incX = 1;
   float Y[] = { 0.619, -0.928, -0.578, 0.48 };
   int incY = -1;
   float y_expected[] = { 0.0417, -0.259758, 0.030276, -0.170742 };
   cblas_cgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "cgbmv(case 553) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "cgbmv(case 553) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int M = 3;
   int N = 2;
   int KL = 1;
   int KU = 1;
   int lda = 3;
   double alpha[2] = {-1, 0};
   double beta[2] = {0, 0};
   double A[] = { 0.683, 0.955, -0.178, -0.954, 0.78, 0.123, -0.111, -0.84, -0.349, 0.81, -0.717, 0.496, 0.405, 0.627, 0.317, 0.679, -0.401, 0.85 };
   double X[] = { 0.841, -0.247, 0.583, -0.481 };
   int incX = 1;
   double Y[] = { 0.204, 0.938, -0.413, 0.453, 0.784, 0.028 };
   int incY = -1;
   double y_expected[] = { -0.51141, -0.24338, -0.681617, -1.067371, 0.863479, 1.149587 };
   cblas_zgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zgbmv(case 554) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zgbmv(case 554) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int M = 3;
   int N = 2;
   int KL = 1;
   int KU = 1;
   int lda = 3;
   double alpha[2] = {-1, 0};
   double beta[2] = {0, 0};
   double A[] = { 0.683, 0.955, -0.178, -0.954, 0.78, 0.123, -0.111, -0.84, -0.349, 0.81, -0.717, 0.496, 0.405, 0.627, 0.317, 0.679, -0.401, 0.85 };
   double X[] = { 0.841, -0.247, 0.583, -0.481 };
   int incX = 1;
   double Y[] = { 0.204, 0.938, -0.413, 0.453, 0.784, 0.028 };
   int incY = -1;
   double y_expected[] = { 0.179435, -0.634045, -0.872504, -0.550882, 0.854089, 1.194677 };
   cblas_zgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zgbmv(case 555) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zgbmv(case 555) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int M = 3;
   int N = 2;
   int KL = 1;
   int KU = 1;
   int lda = 3;
   double alpha[2] = {0, 0};
   double beta[2] = {-0.3, 0.1};
   double A[] = { -0.699, -0.057, -0.327, -0.929, -0.766, -0.953, 0.026, 0.644, 0.055, 0.434, -0.835, 0.461, -0.203, -0.286, -0.949, -0.236, 0.5, -0.868 };
   double X[] = { 0.629, 0.221, 0.386, 0.197, 0.442, 0.406 };
   int incX = 1;
   double Y[] = { -0.66, 0.553, 0.727, 0.02 };
   int incY = -1;
   double y_expected[] = { 0.1427, -0.2319, -0.2201, 0.0667 };
   cblas_zgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zgbmv(case 556) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zgbmv(case 556) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int M = 3;
   int N = 2;
   int KL = 1;
   int KU = 1;
   int lda = 3;
   double alpha[2] = {0, 0};
   double beta[2] = {-0.3, 0.1};
   double A[] = { -0.699, -0.057, -0.327, -0.929, -0.766, -0.953, 0.026, 0.644, 0.055, 0.434, -0.835, 0.461, -0.203, -0.286, -0.949, -0.236, 0.5, -0.868 };
   double X[] = { 0.629, 0.221, 0.386, 0.197, 0.442, 0.406 };
   int incX = 1;
   double Y[] = { -0.66, 0.553, 0.727, 0.02 };
   int incY = -1;
   double y_expected[] = { 0.1427, -0.2319, -0.2201, 0.0667 };
   cblas_zgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zgbmv(case 557) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zgbmv(case 557) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int M = 3;
   int N = 2;
   int KL = 1;
   int KU = 1;
   int lda = 3;
   double alpha[2] = {-1, 0};
   double beta[2] = {0, 1};
   double A[] = { -0.006, -0.192, 0.621, -0.974, 0.546, -0.558, 0.591, -0.938, 0.601, 0.167, -0.799, -0.456, 0.267, 0.745, -0.308, -0.236, -0.952, -0.712 };
   double X[] = { 0.299, 0.568, 0.824, 0.012, -0.137, 0.526 };
   int incX = 1;
   double Y[] = { 0.341, -0.142, -0.856, -0.1 };
   int incY = -1;
   double y_expected[] = { 0.094265, 0.021782, 0.227127, -0.861474 };
   cblas_zgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zgbmv(case 558) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zgbmv(case 558) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int M = 3;
   int N = 2;
   int KL = 1;
   int KU = 1;
   int lda = 3;
   double alpha[2] = {-1, 0};
   double beta[2] = {0, 1};
   double A[] = { -0.006, -0.192, 0.621, -0.974, 0.546, -0.558, 0.591, -0.938, 0.601, 0.167, -0.799, -0.456, 0.267, 0.745, -0.308, -0.236, -0.952, -0.712 };
   double X[] = { 0.299, 0.568, 0.824, 0.012, -0.137, 0.526 };
   int incX = 1;
   double Y[] = { 0.341, -0.142, -0.856, -0.1 };
   int incY = -1;
   double y_expected[] = { 0.13124, 0.337992, 0.024345, -1.966298 };
   cblas_zgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zgbmv(case 559) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zgbmv(case 559) imag");
     };
   };
  };


}
