#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

#include "tests.h"

void
test_syrk (void) {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int uplo = 121;
   int trans = 111;
   int N = 2;
   int K = 1;
   float alpha = -1;
   float beta = 0.1;
   float A[] = { 0.412, -0.229 };
   int lda = 1;
   float C[] = { 0.628, -0.664, -0.268, 0.096 };
   int ldc = 2;
   float C_expected[] = { -0.106944, 0.027948, -0.268, -0.042841 };
   cblas_ssyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssyrk(case 1566)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 111;
   int N = 2;
   int K = 1;
   float alpha = -1;
   float beta = 0.1;
   float A[] = { 0.101, -0.653 };
   int lda = 2;
   float C[] = { 0.432, 0.107, -0.952, -0.532 };
   int ldc = 2;
   float C_expected[] = { 0.032999, 0.107, -0.029247, -0.479609 };
   cblas_ssyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssyrk(case 1567)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 112;
   int N = 2;
   int K = 1;
   float alpha = 1;
   float beta = 0.1;
   float A[] = { 0.79, 0.595 };
   int lda = 2;
   float C[] = { 0.257, 0.183, -0.021, -0.053 };
   int ldc = 2;
   float C_expected[] = { 0.6498, 0.48835, -0.021, 0.348725 };
   cblas_ssyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssyrk(case 1568)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 112;
   int N = 2;
   int K = 1;
   float alpha = 1;
   float beta = 0.1;
   float A[] = { -0.181, -0.654 };
   int lda = 1;
   float C[] = { -0.4, 0.615, 0.147, -0.163 };
   int ldc = 2;
   float C_expected[] = { -0.007239, 0.615, 0.133074, 0.411416 };
   cblas_ssyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssyrk(case 1569)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 111;
   int N = 2;
   int K = 1;
   float alpha = 0;
   float beta = -1;
   float A[] = { -0.191, 0.584 };
   int lda = 1;
   float C[] = { -0.719, -0.681, -0.003, 0.544 };
   int ldc = 2;
   float C_expected[] = { 0.719, -0.681, 0.003, -0.544 };
   cblas_ssyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssyrk(case 1570)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 111;
   int N = 2;
   int K = 1;
   float alpha = 0;
   float beta = -1;
   float A[] = { 0.788, 0.041 };
   int lda = 2;
   float C[] = { 0.029, 0.365, 0.739, -0.769 };
   int ldc = 2;
   float C_expected[] = { -0.029, -0.365, 0.739, 0.769 };
   cblas_ssyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssyrk(case 1571)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 112;
   int N = 2;
   int K = 1;
   float alpha = -0.3;
   float beta = -1;
   float A[] = { 0.733, 0.678 };
   int lda = 2;
   float C[] = { -0.941, 0.96, 0.07, -0.295 };
   int ldc = 2;
   float C_expected[] = { 0.779813, 0.96, -0.219092, 0.157095 };
   cblas_ssyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssyrk(case 1572)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 112;
   int N = 2;
   int K = 1;
   float alpha = -0.3;
   float beta = -1;
   float A[] = { -0.87, 0.675 };
   int lda = 1;
   float C[] = { -0.602, -0.432, -0.984, 0.384 };
   int ldc = 2;
   float C_expected[] = { 0.37493, 0.608175, -0.984, -0.520687 };
   cblas_ssyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssyrk(case 1573)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 111;
   int N = 2;
   int K = 1;
   double alpha = 0.1;
   double beta = -0.3;
   double A[] = { 0.169, -0.875 };
   int lda = 1;
   double C[] = { 0.159, 0.277, 0.865, 0.346 };
   int ldc = 2;
   double C_expected[] = { -0.0448439, -0.0978875, 0.865, -0.0272375 };
   cblas_dsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsyrk(case 1574)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 111;
   int N = 2;
   int K = 1;
   double alpha = 0.1;
   double beta = -0.3;
   double A[] = { 0.536, -0.725 };
   int lda = 2;
   double C[] = { 0.154, -0.445, -0.841, -0.91 };
   int ldc = 2;
   double C_expected[] = { -0.0174704, -0.445, 0.21344, 0.3255625 };
   cblas_dsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsyrk(case 1575)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 112;
   int N = 2;
   int K = 1;
   double alpha = 0;
   double beta = -1;
   double A[] = { -0.07, 0.8 };
   int lda = 2;
   double C[] = { 0.823, -0.88, -0.136, 0.793 };
   int ldc = 2;
   double C_expected[] = { -0.823, 0.88, -0.136, -0.793 };
   cblas_dsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsyrk(case 1576)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 112;
   int N = 2;
   int K = 1;
   double alpha = 0;
   double beta = -1;
   double A[] = { -0.058, 0.649 };
   int lda = 1;
   double C[] = { -0.187, 0.294, -0.004, -0.933 };
   int ldc = 2;
   double C_expected[] = { 0.187, 0.294, 0.004, 0.933 };
   cblas_dsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsyrk(case 1577)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 111;
   int N = 2;
   int K = 1;
   double alpha = 1;
   double beta = -1;
   double A[] = { 0.263, -0.289 };
   int lda = 1;
   double C[] = { 0.554, -0.679, 0.993, 0.758 };
   int ldc = 2;
   double C_expected[] = { -0.484831, -0.679, -1.069007, -0.674479 };
   cblas_dsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsyrk(case 1578)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 111;
   int N = 2;
   int K = 1;
   double alpha = 1;
   double beta = -1;
   double A[] = { -0.265, -0.837 };
   int lda = 2;
   double C[] = { -0.994, 0.967, -0.34, -0.069 };
   int ldc = 2;
   double C_expected[] = { 1.064225, -0.745195, -0.34, 0.769569 };
   cblas_dsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsyrk(case 1579)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 112;
   int N = 2;
   int K = 1;
   double alpha = -0.3;
   double beta = 1;
   double A[] = { -0.464, 0.394 };
   int lda = 2;
   double C[] = { -0.45, -0.447, 0.649, 0.055 };
   int ldc = 2;
   double C_expected[] = { -0.5145888, -0.447, 0.7038448, 0.0084292 };
   cblas_dsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsyrk(case 1580)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 112;
   int N = 2;
   int K = 1;
   double alpha = -0.3;
   double beta = 1;
   double A[] = { 0.815, 0.168 };
   int lda = 1;
   double C[] = { 0.817, -0.957, -0.395, -0.382 };
   int ldc = 2;
   double C_expected[] = { 0.6177325, -0.998076, -0.395, -0.3904672 };
   cblas_dsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsyrk(case 1581)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 111;
   int N = 2;
   int K = 1;
   float alpha[2] = {0, 0};
   float beta[2] = {-0.3, 0.1};
   float A[] = { 0.447, -0.507, -0.425, 0.701 };
   int lda = 1;
   float C[] = { 0.16, -0.245, 0.922, -0.437, 0.24, 0.008, -0.095, 0.749 };
   int ldc = 2;
   float C_expected[] = { -0.0235, 0.0895, -0.2329, 0.2233, 0.24, 0.008, -0.0464, -0.2342 };
   cblas_csyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyrk(case 1582) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyrk(case 1582) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 111;
   int N = 2;
   int K = 1;
   float alpha[2] = {0, 0};
   float beta[2] = {-0.3, 0.1};
   float A[] = { -0.421, -0.435, -0.914, -0.493 };
   int lda = 2;
   float C[] = { -0.761, -0.38, 0.043, -0.999, 0.779, 0.238, 0.082, 0.394 };
   int ldc = 2;
   float C_expected[] = { 0.2663, 0.0379, 0.043, -0.999, -0.2575, 0.0065, -0.064, -0.11 };
   cblas_csyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyrk(case 1583) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyrk(case 1583) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 112;
   int N = 2;
   int K = 1;
   float alpha[2] = {-1, 0};
   float beta[2] = {-0.3, 0.1};
   float A[] = { 0.827, -0.896, 0.417, 0.865 };
   int lda = 2;
   float C[] = { -0.349, -0.31, 0.972, 0.794, -0.906, -0.595, -0.089, -0.333 };
   int ldc = 2;
   float C_expected[] = { 0.254587, 1.54008, -1.4909, -0.482723, -0.906, -0.595, 0.634336, -0.63041 };
   cblas_csyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyrk(case 1584) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyrk(case 1584) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 112;
   int N = 2;
   int K = 1;
   float alpha[2] = {-1, 0};
   float beta[2] = {-0.3, 0.1};
   float A[] = { 0.607, 0.747, -0.889, 0.333 };
   int lda = 1;
   float C[] = { 0.244, 0.564, 0.009, 0.578, -0.827, 0.558, -0.337, 0.731 };
   int ldc = 2;
   float C_expected[] = { 0.05996, -1.05166, 0.009, 0.578, 0.980674, 0.211852, -0.651432, 0.339074 };
   cblas_csyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyrk(case 1585) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyrk(case 1585) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 113;
   int N = 2;
   int K = 1;
   float alpha[2] = {1, 0};
   float beta[2] = {0, 1};
   float A[] = { 0.784, -0.281, -0.88, 0.479 };
   int lda = 2;
   float C[] = { 0.491, 0.531, 0.805, -0.097, 0.728, 0.674, -0.705, -0.754 };
   int ldc = 2;
   float C_expected[] = { 0.004695, 0.050392, -0.458321, 1.42782, 0.728, 0.674, 1.29896, -1.54804 };
   cblas_csyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyrk(case 1586) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyrk(case 1586) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 113;
   int N = 2;
   int K = 1;
   float alpha[2] = {1, 0};
   float beta[2] = {0, 1};
   float A[] = { 0.272, -0.146, 0.155, 0.038 };
   int lda = 1;
   float C[] = { 0.533, -0.41, -0.904, 0.301, -0.836, 0.57, -0.374, -0.293 };
   int ldc = 2;
   float C_expected[] = { 0.462668, 0.453576, -0.904, 0.301, -0.522292, -0.848294, 0.315581, -0.36222 };
   cblas_csyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyrk(case 1587) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyrk(case 1587) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 111;
   int N = 2;
   int K = 1;
   float alpha[2] = {0, 1};
   float beta[2] = {-1, 0};
   float A[] = { -0.055, -0.127, -0.896, -0.625 };
   int lda = 1;
   float C[] = { -0.619, 0.511, -0.877, 0.557, -0.801, -0.437, -0.922, 0.332 };
   int ldc = 2;
   float C_expected[] = { 0.60503, -0.524104, -0.877, 0.557, 0.652833, 0.406905, -0.198, 0.080191 };
   cblas_csyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyrk(case 1588) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyrk(case 1588) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 111;
   int N = 2;
   int K = 1;
   float alpha[2] = {0, 1};
   float beta[2] = {-1, 0};
   float A[] = { -0.528, 0.759, -0.079, 0.952 };
   int lda = 2;
   float C[] = { 0.775, 0.855, 0.786, 0.525, 0.85, 0.044, 0.658, 0.947 };
   int ldc = 2;
   float C_expected[] = { 0.026504, -1.1523, -0.223383, -1.20586, 0.85, 0.044, -0.507584, -1.84706 };
   cblas_csyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyrk(case 1589) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyrk(case 1589) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 112;
   int N = 2;
   int K = 1;
   float alpha[2] = {1, 0};
   float beta[2] = {1, 0};
   float A[] = { -0.049, -0.687, -0.434, 0.294 };
   int lda = 2;
   float C[] = { 0.937, -0.113, 0.796, 0.293, 0.876, -0.199, -0.757, -0.103 };
   int ldc = 2;
   float C_expected[] = { 0.467432, -0.045674, 0.796, 0.293, 1.09924, 0.084752, -0.65508, -0.358192 };
   cblas_csyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyrk(case 1590) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyrk(case 1590) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 112;
   int N = 2;
   int K = 1;
   float alpha[2] = {1, 0};
   float beta[2] = {1, 0};
   float A[] = { 0.359, -0.364, 0.926, -0.69 };
   int lda = 1;
   float C[] = { 0.306, 0.249, 0.28, 0.229, 0.866, 0.092, 0.886, -0.283 };
   int ldc = 2;
   float C_expected[] = { 0.302385, -0.012352, 0.361274, -0.355774, 0.866, 0.092, 1.26738, -1.56088 };
   cblas_csyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyrk(case 1591) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyrk(case 1591) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 113;
   int N = 2;
   int K = 1;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {0, 0};
   float A[] = { 0.607, 0.555, -0.85, 0.831 };
   int lda = 2;
   float C[] = { 0.069, 0.368, 0.551, -0.912, -0.243, -0.063, -0.924, 0.192 };
   int ldc = 2;
   float C_expected[] = { -0.0855042, -0.196089, 0.551, -0.912, 0.28988, -0.107516, 0.131688, 0.427004 };
   cblas_csyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyrk(case 1592) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyrk(case 1592) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 113;
   int N = 2;
   int K = 1;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {0, 0};
   float A[] = { 0.427, 0.86, -0.136, 0.002 };
   int lda = 1;
   float C[] = { 0.398, -0.47, 0.011, -0.547, -0.106, 0.016, 0.681, 0.246 };
   int ldc = 2;
   float C_expected[] = { 0.0937373, -0.276059, 0.0295482, 0.0288526, -0.106, 0.016, -0.0054932, 0.0020124 };
   cblas_csyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyrk(case 1593) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyrk(case 1593) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 111;
   int N = 2;
   int K = 1;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {1, 0};
   double A[] = { 0.718, 0.023, 0.355, -0.492 };
   int lda = 1;
   double C[] = { -0.637, -0.727, -0.475, -0.776, 0.802, -0.55, -0.837, 0.222 };
   int ldc = 2;
   double C_expected[] = { -0.7948013, -0.6854089, -0.5203527, -0.6458521, 0.802, -0.55, -0.7672563, 0.3151921 };
   cblas_zsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyrk(case 1594) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyrk(case 1594) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 111;
   int N = 2;
   int K = 1;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {1, 0};
   double A[] = { 0.209, 0.139, -0.202, -0.223 };
   int lda = 2;
   double C[] = { -0.695, 0.524, 0.212, -0.88, -0.752, 0.291, 0.684, -0.124 };
   int ldc = 2;
   double C_expected[] = { -0.7081182, 0.5090054, 0.212, -0.88, -0.7411652, 0.3122834, 0.6776683, -0.1519201 };
   cblas_zsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyrk(case 1595) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyrk(case 1595) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 112;
   int N = 2;
   int K = 1;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {1, 0};
   double A[] = { -0.365, -0.624, 0.632, 0.348 };
   int lda = 2;
   double C[] = { 0.877, 0.927, -0.377, 0.967, 0.008, 0.292, -0.779, 0.794 };
   int ldc = 2;
   double C_expected[] = { 0.9082933, 0.7647289, -0.3208028, 1.1220636, 0.008, 0.292, -0.9064832, 0.6898704 };
   cblas_zsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyrk(case 1596) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyrk(case 1596) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 112;
   int N = 2;
   int K = 1;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {1, 0};
   double A[] = { -0.067, -0.586, 0.208, 0.331 };
   int lda = 1;
   double C[] = { 0.584, -0.454, 0.93, 0.782, 0.489, -0.278, 0.081, -0.919 };
   int ldc = 2;
   double C_expected[] = { 0.6778197, -0.5114479, 0.93, 0.782, 0.4493975, -0.2167775, 0.0871195, -0.9669385 };
   cblas_zsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyrk(case 1597) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyrk(case 1597) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 113;
   int N = 2;
   int K = 1;
   double alpha[2] = {0, 0.1};
   double beta[2] = {0, 0.1};
   double A[] = { -0.617, 0.179, -0.626, 0.334 };
   int lda = 2;
   double C[] = { 0.346, -0.903, 0.022, -0.839, -0.715, 0.049, -0.338, 0.149 };
   int ldc = 2;
   double C_expected[] = { 0.1123886, 0.0694648, 0.1157132, 0.0348456, -0.715, 0.049, 0.0269168, -0.005768 };
   cblas_zsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyrk(case 1598) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyrk(case 1598) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 113;
   int N = 2;
   int K = 1;
   double alpha[2] = {0, 0.1};
   double beta[2] = {0, 0.1};
   double A[] = { -0.356, -0.308, 0.493, -0.351 };
   int lda = 1;
   double C[] = { -0.898, -0.905, 0.002, -0.219, 0.881, 0.879, 0.275, -0.351 };
   int ldc = 2;
   double C_expected[] = { 0.0685704, -0.0866128, 0.002, -0.219, -0.0852112, 0.0597384, 0.0697086, 0.0394848 };
   cblas_zsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyrk(case 1599) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyrk(case 1599) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 111;
   int N = 2;
   int K = 1;
   double alpha[2] = {0, 0.1};
   double beta[2] = {1, 0};
   double A[] = { -0.103, -0.951, -0.601, -0.041 };
   int lda = 1;
   double C[] = { -0.918, -0.018, 0.991, -0.789, -0.698, -0.067, 0.956, -0.599 };
   int ldc = 2;
   double C_expected[] = { -0.9375906, -0.1073792, 0.991, -0.789, -0.7555774, -0.0647088, 0.9510718, -0.563048 };
   cblas_zsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyrk(case 1600) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyrk(case 1600) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 111;
   int N = 2;
   int K = 1;
   double alpha[2] = {0, 0.1};
   double beta[2] = {1, 0};
   double A[] = { -0.237, 0.925, -0.904, -0.091 };
   int lda = 2;
   double C[] = { -0.572, 0.915, 0.398, 0.222, 0.016, 0.288, -0.078, -0.507 };
   int ldc = 2;
   double C_expected[] = { -0.528155, 0.8350544, 0.4794633, 0.2518423, 0.016, 0.288, -0.0944528, -0.4261065 };
   cblas_zsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyrk(case 1601) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyrk(case 1601) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 112;
   int N = 2;
   int K = 1;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {0, 0};
   double A[] = { 0.963, -0.23, -0.435, 0.289 };
   int lda = 2;
   double C[] = { 0.282, -0.272, -0.516, -0.594, -0.001, 0.155, -0.39, -0.354 };
   int ldc = 2;
   double C_expected[] = { -0.2180427, 0.2203409, -0.516, -0.594, 0.0678948, -0.1487506, -0.0065682, 0.0859994 };
   cblas_zsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyrk(case 1602) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyrk(case 1602) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 112;
   int N = 2;
   int K = 1;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {0, 0};
   double A[] = { 0.674, 0.1, -0.098, 0.552 };
   int lda = 1;
   double C[] = { 0.089, -0.523, -0.551, 0.618, 0.67, 0.247, 0.975, -0.714 };
   int ldc = 2;
   double C_expected[] = { -0.1467628, 0.0039876, 0.0001508, -0.1207996, 0.67, 0.247, 0.0993492, 0.0029476 };
   cblas_zsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyrk(case 1603) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyrk(case 1603) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 113;
   int N = 2;
   int K = 1;
   double alpha[2] = {0, 1};
   double beta[2] = {1, 0};
   double A[] = { 0.033, -0.864, 0.168, 0.524 };
   int lda = 2;
   double C[] = { 0.788, 0.016, -0.436, 0.749, -0.89, -0.87, 0.421, -0.203 };
   int ldc = 2;
   double C_expected[] = { 0.845024, -0.729407, -0.436, 0.749, -0.76214, -0.41172, 0.244936, -0.449352 };
   cblas_zsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyrk(case 1604) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyrk(case 1604) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 113;
   int N = 2;
   int K = 1;
   double alpha[2] = {0, 1};
   double beta[2] = {1, 0};
   double A[] = { 0.957, -0.079, 0.935, 0.232 };
   int lda = 1;
   double C[] = { -0.744, -0.061, 0.195, -0.574, 0.551, 0.478, -0.337, 0.1 };
   int ldc = 2;
   double C_expected[] = { -0.592794, 0.848608, 0.046841, 0.339123, 0.551, 0.478, -0.77084, 0.920401 };
   cblas_zsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyrk(case 1605) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyrk(case 1605) imag");
     };
   };
  };


}
