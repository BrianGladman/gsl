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
   int N = 4;
   int KL = 1;
   int KU = 1;
   int lda = 4;
   float alpha = -1;
   float beta = -1;
   float A[] = { 0.423, -0.143, -0.182, -0.076, -0.855, 0.599, 0.389, -0.473, 0.493, -0.902, -0.889, -0.256, 0.112, 0.128, -0.277, -0.777 };
   float X[] = { 0.488, 0.029, -0.633, 0.84 };
   int incX = -1;
   float Y[] = { 0.874, 0.322, -0.477 };
   int incY = -1;
   float y_expected[] = { -0.101941, 0.764086, 0.481914 };
   cblas_sgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sgbmv(case 794)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int M = 3;
   int N = 4;
   int KL = 1;
   int KU = 1;
   int lda = 4;
   float alpha = -1;
   float beta = -1;
   float A[] = { 0.423, -0.143, -0.182, -0.076, -0.855, 0.599, 0.389, -0.473, 0.493, -0.902, -0.889, -0.256, 0.112, 0.128, -0.277, -0.777 };
   float X[] = { 0.488, 0.029, -0.633, 0.84 };
   int incX = -1;
   float Y[] = { 0.874, 0.322, -0.477 };
   int incY = -1;
   float y_expected[] = { -0.656261, 0.19575, 0.055905 };
   cblas_sgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sgbmv(case 795)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int M = 3;
   int N = 4;
   int KL = 1;
   int KU = 1;
   int lda = 4;
   float alpha = 0;
   float beta = 0.1;
   float A[] = { -0.066, -0.153, -0.619, 0.174, 0.777, 0.543, 0.614, -0.446, -0.138, -0.767, 0.725, 0.222, 0.165, -0.063, -0.047, 0.267 };
   float X[] = { -0.096, -0.007, -0.657 };
   int incX = -1;
   float Y[] = { -0.88, 0.102, -0.278, 0.403 };
   int incY = -1;
   float y_expected[] = { -0.088, 0.0102, -0.0278, 0.0403 };
   cblas_sgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sgbmv(case 796)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int M = 3;
   int N = 4;
   int KL = 1;
   int KU = 1;
   int lda = 4;
   float alpha = 0;
   float beta = 0.1;
   float A[] = { -0.066, -0.153, -0.619, 0.174, 0.777, 0.543, 0.614, -0.446, -0.138, -0.767, 0.725, 0.222, 0.165, -0.063, -0.047, 0.267 };
   float X[] = { -0.096, -0.007, -0.657 };
   int incX = -1;
   float Y[] = { -0.88, 0.102, -0.278, 0.403 };
   int incY = -1;
   float y_expected[] = { -0.088, 0.0102, -0.0278, 0.0403 };
   cblas_sgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sgbmv(case 797)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int M = 3;
   int N = 4;
   int KL = 1;
   int KU = 1;
   int lda = 4;
   double alpha = 0.1;
   double beta = 0;
   double A[] = { -0.688, 0.29, 0.442, -0.001, 0.313, -0.073, 0.991, -0.654, -0.12, 0.416, 0.571, 0.932, -0.179, -0.724, 0.492, -0.965 };
   double X[] = { 0.187, -0.338, -0.976, -0.052 };
   int incX = -1;
   double Y[] = { -0.101, 0.8, 0.026 };
   int incY = -1;
   double y_expected[] = { 0.0083289, -0.0279986, -0.0446472 };
   cblas_dgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dgbmv(case 798)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int M = 3;
   int N = 4;
   int KL = 1;
   int KU = 1;
   int lda = 4;
   double alpha = 0.1;
   double beta = 0;
   double A[] = { -0.688, 0.29, 0.442, -0.001, 0.313, -0.073, 0.991, -0.654, -0.12, 0.416, 0.571, 0.932, -0.179, -0.724, 0.492, -0.965 };
   double X[] = { 0.187, -0.338, -0.976, -0.052 };
   int incX = -1;
   double Y[] = { -0.101, 0.8, 0.026 };
   int incY = -1;
   double y_expected[] = { -0.1141297, 0.0088824, -0.0320568 };
   cblas_dgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dgbmv(case 799)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int M = 3;
   int N = 4;
   int KL = 1;
   int KU = 1;
   int lda = 4;
   double alpha = -0.3;
   double beta = -0.3;
   double A[] = { 0.746, 0.262, -0.449, -0.954, -0.093, 0.108, -0.496, 0.927, 0.177, 0.729, -0.92, -0.469, 0.87, -0.877, -0.308, -0.806 };
   double X[] = { 0.662, -0.887, 0.261 };
   int incX = -1;
   double Y[] = { 0.771, 0.637, -0.177, -0.018 };
   int incY = -1;
   double y_expected[] = { -0.048588, -0.467865, 0.0818433, -0.0398619 };
   cblas_dgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dgbmv(case 800)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int M = 3;
   int N = 4;
   int KL = 1;
   int KU = 1;
   int lda = 4;
   double alpha = -0.3;
   double beta = -0.3;
   double A[] = { 0.746, 0.262, -0.449, -0.954, -0.093, 0.108, -0.496, 0.927, 0.177, 0.729, -0.92, -0.469, 0.87, -0.877, -0.308, -0.806 };
   double X[] = { 0.662, -0.887, 0.261 };
   int incX = -1;
   double Y[] = { 0.771, 0.637, -0.177, -0.018 };
   int incY = -1;
   double y_expected[] = { -0.404082, -0.2887797, 0.1876263, -0.1345935 };
   cblas_dgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dgbmv(case 801)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int M = 3;
   int N = 4;
   int KL = 1;
   int KU = 1;
   int lda = 4;
   float alpha[2] = {0, 1};
   float beta[2] = {0, 0.1};
   float A[] = { -0.107, 0.926, -0.246, -0.555, -0.301, 0.276, 0.471, -0.084, -0.754, 0.082, -0.952, -0.394, 0.659, 0.054, 0.795, 0.923, 0.232, -0.788, 0.478, 0.775, -0.118, 0.691, -0.933, 0.809, 0.164, -0.263, -0.923, -0.88, 0.819, -0.521, -0.045, 0.034 };
   float X[] = { 0.407, 0.895, 0.301, 0.769, -0.269, -0.465, 0.455, -0.628 };
   int incX = -1;
   float Y[] = { -0.116, -0.744, -0.936, -0.064, -0.232, -0.665 };
   int incY = -1;
   float y_expected[] = { -0.806176, -1.559, -1.57611, -0.155463, 0.098816, -0.274361 };
   cblas_cgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "cgbmv(case 802) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "cgbmv(case 802) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int M = 3;
   int N = 4;
   int KL = 1;
   int KU = 1;
   int lda = 4;
   float alpha[2] = {0, 1};
   float beta[2] = {0, 0.1};
   float A[] = { -0.107, 0.926, -0.246, -0.555, -0.301, 0.276, 0.471, -0.084, -0.754, 0.082, -0.952, -0.394, 0.659, 0.054, 0.795, 0.923, 0.232, -0.788, 0.478, 0.775, -0.118, 0.691, -0.933, 0.809, 0.164, -0.263, -0.923, -0.88, 0.819, -0.521, -0.045, 0.034 };
   float X[] = { 0.407, 0.895, 0.301, 0.769, -0.269, -0.465, 0.455, -0.628 };
   int incX = -1;
   float Y[] = { -0.116, -0.744, -0.936, -0.064, -0.232, -0.665 };
   int incY = -1;
   float y_expected[] = { -0.245235, -0.313725, -0.798094, 0.691455, -0.164015, -0.242714 };
   cblas_cgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "cgbmv(case 803) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "cgbmv(case 803) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int M = 3;
   int N = 4;
   int KL = 1;
   int KU = 1;
   int lda = 4;
   float alpha[2] = {-1, 0};
   float beta[2] = {0, 0.1};
   float A[] = { -0.258, 0.838, -0.106, -0.066, 0.395, 0.982, -0.546, 0.565, 0.14, -0.18, 0.165, -0.186, 0.499, -0.038, -0.305, -0.653, -0.811, -0.466, -0.674, -0.013, -0.552, -0.807, -0.536, 0.864, -0.027, -0.606, 0.459, 0.564, -0.968, 0.717, -0.312, -0.485 };
   float X[] = { -0.399, 0.459, 0.398, 0.358, -0.161, -0.359 };
   int incX = -1;
   float Y[] = { 0.572, 0.293, -0.813, -0.096, -0.611, -0.717, 0.736, 0.259 };
   int incY = -1;
   float y_expected[] = { -0.619961, -0.011425, -0.477499, 0.059361, -0.886984, 0.44008, -0.139432, 0.04644 };
   cblas_cgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "cgbmv(case 804) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "cgbmv(case 804) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int M = 3;
   int N = 4;
   int KL = 1;
   int KU = 1;
   int lda = 4;
   float alpha[2] = {-1, 0};
   float beta[2] = {0, 0.1};
   float A[] = { -0.258, 0.838, -0.106, -0.066, 0.395, 0.982, -0.546, 0.565, 0.14, -0.18, 0.165, -0.186, 0.499, -0.038, -0.305, -0.653, -0.811, -0.466, -0.674, -0.013, -0.552, -0.807, -0.536, 0.864, -0.027, -0.606, 0.459, 0.564, -0.968, 0.717, -0.312, -0.485 };
   float X[] = { -0.399, 0.459, 0.398, 0.358, -0.161, -0.359 };
   int incX = -1;
   float Y[] = { 0.572, 0.293, -0.813, -0.096, -0.611, -0.717, 0.736, 0.259 };
   int incY = -1;
   float y_expected[] = { -0.318227, -0.172201, -0.109343, 0.698685, 0.208261, -0.269065, 0.175074, -0.507326 };
   cblas_cgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "cgbmv(case 805) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "cgbmv(case 805) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int M = 3;
   int N = 4;
   int KL = 1;
   int KU = 1;
   int lda = 4;
   float alpha[2] = {-1, 0};
   float beta[2] = {0, 0.1};
   float A[] = { -0.804, 0.232, -0.448, -0.558, -0.078, -0.056, -0.345, -0.379, 0.369, -0.662, -0.169, -0.391, -0.215, 0.467, 0.374, 0.889, -0.698, 0.734, 0.377, -0.955, 0.498, 0.151, -0.725, -0.728, -0.655, -0.581, 0.389, 0.949, -0.553, -0.434, 0.237, 0.641 };
   float X[] = { -0.262, -0.823, -0.357, -0.994, -0.347, -0.375 };
   int incX = -1;
   float Y[] = { -0.683, -0.87, -0.708, 0.071, 0.575, -0.575, 0.845, 0.032 };
   int incY = -1;
   float y_expected[] = { 0.341749, 0.301992, -0.306848, 0.109252, -0.018347, -0.747479, -0.894201, 0.713246 };
   cblas_cgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "cgbmv(case 806) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "cgbmv(case 806) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int M = 3;
   int N = 4;
   int KL = 1;
   int KU = 1;
   int lda = 4;
   float alpha[2] = {-1, 0};
   float beta[2] = {0, 0.1};
   float A[] = { -0.804, 0.232, -0.448, -0.558, -0.078, -0.056, -0.345, -0.379, 0.369, -0.662, -0.169, -0.391, -0.215, 0.467, 0.374, 0.889, -0.698, 0.734, 0.377, -0.955, 0.498, 0.151, -0.725, -0.728, -0.655, -0.581, 0.389, 0.949, -0.553, -0.434, 0.237, 0.641 };
   float X[] = { -0.262, -0.823, -0.357, -0.994, -0.347, -0.375 };
   int incX = -1;
   float Y[] = { -0.683, -0.87, -0.708, 0.071, 0.575, -0.575, 0.845, 0.032 };
   int incY = -1;
   float y_expected[] = { -0.562773, -0.455143, -0.213881, -0.466169, -0.183683, 0.097891, -0.451416, 0.052586 };
   cblas_cgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "cgbmv(case 807) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "cgbmv(case 807) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int M = 3;
   int N = 4;
   int KL = 1;
   int KU = 1;
   int lda = 4;
   double alpha[2] = {0, 0.1};
   double beta[2] = {1, 0};
   double A[] = { -0.919, -0.002, 0.105, -0.338, -0.358, -0.715, -0.157, 0.307, 0.334, 0.121, 0.366, 0.029, -0.006, -0.662, -0.314, 0.061, -0.322, -0.865, -0.586, 0.556, 0.507, 0.581, 0.855, -0.09, 0.836, -0.788, -0.209, -0.694, -0.695, 0.11, -0.234, 0.17 };
   double X[] = { 0.356, -0.76, -0.96, 0.437, -0.849, 0.397, -0.382, -0.826 };
   int incX = -1;
   double Y[] = { 0.288, -0.832, 0.889, 0.576, -0.809, 0.4 };
   int incY = -1;
   double y_expected[] = { 0.3241775, -0.6761577, 0.8458527, 0.5705165, -0.8597295, 0.4268499 };
   cblas_zgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zgbmv(case 808) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zgbmv(case 808) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int M = 3;
   int N = 4;
   int KL = 1;
   int KU = 1;
   int lda = 4;
   double alpha[2] = {0, 0.1};
   double beta[2] = {1, 0};
   double A[] = { -0.919, -0.002, 0.105, -0.338, -0.358, -0.715, -0.157, 0.307, 0.334, 0.121, 0.366, 0.029, -0.006, -0.662, -0.314, 0.061, -0.322, -0.865, -0.586, 0.556, 0.507, 0.581, 0.855, -0.09, 0.836, -0.788, -0.209, -0.694, -0.695, 0.11, -0.234, 0.17 };
   double X[] = { 0.356, -0.76, -0.96, 0.437, -0.849, 0.397, -0.382, -0.826 };
   int incX = -1;
   double Y[] = { 0.288, -0.832, 0.889, 0.576, -0.809, 0.4 };
   int incY = -1;
   double y_expected[] = { 0.4026074, -0.8033768, 0.7510795, 0.5671044, -0.8162255, 0.3349099 };
   cblas_zgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zgbmv(case 809) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zgbmv(case 809) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int M = 3;
   int N = 4;
   int KL = 1;
   int KU = 1;
   int lda = 4;
   double alpha[2] = {1, 0};
   double beta[2] = {1, 0};
   double A[] = { 0.511, -0.707, -0.906, 0.345, -0.524, -0.933, 0.154, -0.529, -0.651, -0.851, 0.104, 0.532, -0.297, 0.477, 0.511, 0.469, -0.888, -0.789, 0.656, 0.288, -0.749, 0.961, 0.571, 0.539, 0.465, 0.647, 0.653, -0.994, -0.515, 0.297, 0.35, -0.707 };
   double X[] = { -0.991, 0.658, -0.909, -0.99, -0.517, -0.071 };
   int incX = -1;
   double Y[] = { 0.451, 0.351, -0.113, -0.62, 0.983, 0.511, 0.142, -0.186 };
   int incY = -1;
   double y_expected[] = { 0.560921, -1.094193, -0.210397, -0.613323, 3.018979, 0.641612, 0.384166, 1.11801 };
   cblas_zgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zgbmv(case 810) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zgbmv(case 810) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int M = 3;
   int N = 4;
   int KL = 1;
   int KU = 1;
   int lda = 4;
   double alpha[2] = {1, 0};
   double beta[2] = {1, 0};
   double A[] = { 0.511, -0.707, -0.906, 0.345, -0.524, -0.933, 0.154, -0.529, -0.651, -0.851, 0.104, 0.532, -0.297, 0.477, 0.511, 0.469, -0.888, -0.789, 0.656, 0.288, -0.749, 0.961, 0.571, 0.539, 0.465, 0.647, 0.653, -0.994, -0.515, 0.297, 0.35, -0.707 };
   double X[] = { -0.991, 0.658, -0.909, -0.99, -0.517, -0.071 };
   int incX = -1;
   double Y[] = { 0.451, 0.351, -0.113, -0.62, 0.983, 0.511, 0.142, -0.186 };
   int incY = -1;
   double y_expected[] = { -0.435541, 0.015793, -0.926518, 1.122561, 1.671751, -0.257493, 0.187543, 1.066818 };
   cblas_zgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zgbmv(case 811) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zgbmv(case 811) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int M = 3;
   int N = 4;
   int KL = 1;
   int KU = 1;
   int lda = 4;
   double alpha[2] = {0, 0.1};
   double beta[2] = {-0.3, 0.1};
   double A[] = { 0.534, 0.67, -0.621, 0.143, -0.794, 0.073, 0.414, -0.9, 0.155, -0.368, 0.122, -0.583, 0.03, 0.646, -0.768, -0.892, -0.741, -0.397, 0.626, 0.004, -0.515, 0.355, 0.196, -0.989, -0.982, 0.985, 0.445, 0.63, -0.849, -0.528, 0.146, -0.319 };
   double X[] = { -0.199, -0.259, 0.386, -0.131, -0.867, 0.888 };
   int incX = -1;
   double Y[] = { 0.106, 0.874, 0.962, 0.636, -0.759, 0.415, -0.053, 0.315 };
   int incY = -1;
   double y_expected[] = { -0.139603, -0.250546, -0.3107376, -0.1144656, 0.2181809, -0.0877031, 0.0149724, -0.0224571 };
   cblas_zgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zgbmv(case 812) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zgbmv(case 812) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int M = 3;
   int N = 4;
   int KL = 1;
   int KU = 1;
   int lda = 4;
   double alpha[2] = {0, 0.1};
   double beta[2] = {-0.3, 0.1};
   double A[] = { 0.534, 0.67, -0.621, 0.143, -0.794, 0.073, 0.414, -0.9, 0.155, -0.368, 0.122, -0.583, 0.03, 0.646, -0.768, -0.892, -0.741, -0.397, 0.626, 0.004, -0.515, 0.355, 0.196, -0.989, -0.982, 0.985, 0.445, 0.63, -0.849, -0.528, 0.146, -0.319 };
   double X[] = { -0.199, -0.259, 0.386, -0.131, -0.867, 0.888 };
   int incX = -1;
   double Y[] = { 0.106, 0.874, 0.962, 0.636, -0.759, 0.415, -0.053, 0.315 };
   int incY = -1;
   double y_expected[] = { -0.1642353, -0.2575697, -0.3610975, -0.1305629, 0.1713576, -0.2514988, 0.0195631, -0.0648656 };
   cblas_zgbmv(order, trans, M, N, KU, KL, alpha, A,                                  lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zgbmv(case 813) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zgbmv(case 813) imag");
     };
   };
  };


}
