#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_trsm () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { 0.773, 0.069, 0.45, 0.189 };
   int lda = 2;
   float B[] = { -0.037, 0.788, 0.015, 0.028, -0.804, -0.357 };
   int ldb = 3;
   float B_expected[] = { 0.0183269, -0.419738, -0.0564036, -0.0444444, 1.27619, 0.566667 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1830)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { -0.13, -0.832, 0.426, 0.195 };
   int lda = 2;
   float B[] = { 0.504, 0.996, 0.872, -0.35, 0.518, -0.8 };
   int ldb = 3;
   float B_expected[] = { -0.06384, -0.428093, -0.06192, 0.105, -0.1554, 0.24 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1831)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { 0.755, -0.053, -0.132, -0.515 };
   int lda = 2;
   float B[] = { -0.735, 0.494, 0.072, -0.882, -0.112, 0.904 };
   int ldb = 3;
   float B_expected[] = { 0.292053, -0.196291, -0.0286093, -0.588643, -0.0149311, 0.533935 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1832)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { -0.88, -0.555, 0.642, 0.751 };
   int lda = 2;
   float B[] = { -0.411, 0.134, 0.657, 0.072, -0.007, -0.34 };
   int ldb = 3;
   float B_expected[] = { 0.1233, -0.0402, -0.1971, -0.100759, 0.0279084, 0.228538 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1833)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { -0.478, 0.938, -0.731, 0.25 };
   int lda = 2;
   float B[] = { -0.859, -0.409, -0.154, -0.54, 0.146, -0.106 };
   int ldb = 2;
   float B_expected[] = { -1.2897, 0.4908, -1.08763, 0.648, -0.102894, 0.1272 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1834)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { 0.953, 0.249, -0.451, -0.781 };
   int lda = 2;
   float B[] = { -0.4, -0.546, 0.839, 0.392, -0.445, -0.818 };
   int ldb = 2;
   float B_expected[] = { 0.193874, 0.1638, -0.304738, -0.1176, 0.244175, 0.2454 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1835)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { 0.831, -0.997, -0.366, 0.307 };
   int lda = 2;
   float B[] = { 0.157, -0.02, 0.57, 0.309, -0.159, 0.266 };
   int ldb = 2;
   float B_expected[] = { -0.0566787, -0.164523, -0.205776, -0.970224, 0.0574007, -0.0735227 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1836)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { -0.842, 0.674, 0.03, 0.628 };
   int lda = 2;
   float B[] = { -0.426, 0.806, 0.299, 0.626, -0.471, 0.208 };
   int ldb = 2;
   float B_expected[] = { 0.1278, -0.327937, -0.0897, -0.127342, 0.1413, -0.157636 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1837)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { -0.095, 0.301, 0.168, 0.934, 0.107, 0.068, 0.384, -0.201, 0.116 };
   int lda = 3;
   float B[] = { 0.534, 0.773, -0.304, -0.402, 0.642, -0.102 };
   int ldb = 3;
   float B_expected[] = { 1.68632, -6.91104, 2.39525, -1.26947, 1.77114, 1.06409 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1838)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { -0.738, -0.353, -0.616, 0.304, 0.403, 0.739, 0.996, 0.329, 0.273 };
   int lda = 3;
   float B[] = { -0.436, 0.074, 0.273, -0.609, 0.858, 0.993 };
   int ldb = 3;
   float B_expected[] = { 0.1308, 0.0239724, -0.0190428, 0.1827, -0.192907, -0.0427986 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1839)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { -0.956, 0.878, 0.156, 0.217, 0.082, -0.869, 0.595, 0.845, 0.064 };
   int lda = 3;
   float B[] = { -0.744, 0.662, -0.31, 0.811, 0.257, 0.98 };
   int ldb = 3;
   float B_expected[] = { -3.27779, -17.3962, 1.45312, 7.92713, 46.3978, -4.59375 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1840)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { 0.313, -0.316, 0.836, 0.359, -0.415, 0.154, -0.948, -0.596, -0.799 };
   int lda = 3;
   float B[] = { 0.29, -0.291, 0.652, 0.614, 0.922, -0.063 };
   int ldb = 3;
   float B_expected[] = { -0.261918, -0.0292776, -0.1956, -0.0710273, -0.265336, 0.0189 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1841)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { -0.634, 0.561, 0.883, -0.136, 0.203, -0.531, 0.733, -0.332, 0.705 };
   int lda = 3;
   float B[] = { 0.133, -0.843, -0.179, 0.94, -0.656, 0.645 };
   int ldb = 2;
   float B_expected[] = { 0.0629338, -0.398896, 0.306695, -1.6564, 0.358145, -0.639766 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1842)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { 0.742, -0.438, 0.991, 0.614, 0.108, -0.125, 0.736, -0.383, 0 };
   int lda = 3;
   float B[] = { -0.792, -0.033, -0.723, 0.885, 0.336, 0.584 };
   int ldb = 2;
   float B_expected[] = { 0.2376, 0.0099, 0.0710136, -0.271579, -0.248475, -0.286501 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1843)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { 0.761, 0.466, 0.907, -0.85, -0.342, -0.058, -0.379, -0.416, 0.599 };
   int lda = 3;
   float B[] = { -0.238, 0.013, 0.473, -0.626, 0.912, -0.003 };
   int ldb = 2;
   float B_expected[] = { 0.336709, 0.329497, 0.492375, -0.549378, -0.456761, 0.0015025 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1844)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { 0.567, -0.532, -0.817, 0.85, -0.135, 0.797, 0.981, -0.75, 0.856 };
   int lda = 3;
   float B[] = { -0.705, 0.326, 0.184, 0.079, -0.173, 0.125 };
   int ldb = 2;
   float B_expected[] = { 0.20253, -0.125146, -0.0965643, 0.0061875, 0.0519, -0.0375 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1845)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { -0.859, 0.563, -0.61, 0.2 };
   int lda = 2;
   float B[] = { -0.241, -0.357, -0.683, -0.718, 0.69, -0.486 };
   int ldb = 3;
   float B_expected[] = { -0.0841676, -0.12468, -0.238533, 1.31393, -0.684026, 1.40047 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1846)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { 0.157, -0.741, 0.844, 0.206 };
   int lda = 2;
   float B[] = { 0.816, -0.692, 0.765, -0.408, 0.404, 0.764 };
   int ldb = 3;
   float B_expected[] = { -0.2448, 0.2076, -0.2295, -0.0589968, 0.0326316, -0.399259 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1847)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { 0.187, 0.354, -0.931, 0.18 };
   int lda = 2;
   float B[] = { -0.215, -0.645, 0.847, 0.014, 0.83, 0.761 };
   int ldb = 3;
   float B_expected[] = { 0.228752, -5.85232, -7.67336, -0.0233333, -1.38333, -1.26833 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1848)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { -0.923, 0.27, -0.319, -0.856 };
   int lda = 2;
   float B[] = { 0.391, 0.01, 0.429, 0.685, 0.332, -0.643 };
   int ldb = 3;
   float B_expected[] = { -0.182855, -0.0347724, -0.0671649, -0.2055, -0.0996, 0.1929 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1849)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { 0.724, 0.201, 0.87, -0.638 };
   int lda = 2;
   float B[] = { -0.533, 0.183, 0.569, 0.85, 0.642, -0.051 };
   int ldb = 2;
   float B_expected[] = { 0.220856, 0.387218, -0.235773, 0.0781772, -0.266022, -0.386739 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1850)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { 0.291, 0.244, 0.931, 0.857 };
   int lda = 2;
   float B[] = { 0.008, -0.478, -0.252, -0.155, 0.419, -0.192 };
   int ldb = 2;
   float B_expected[] = { -0.0024, 0.145634, 0.0756, -0.0238836, -0.1257, 0.174627 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1851)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { -0.634, -0.529, -0.344, 0.375 };
   int lda = 2;
   float B[] = { -0.295, 0.551, 0.832, 0.744, -0.326, 0.111 };
   int ldb = 2;
   float B_expected[] = { 0.228207, -0.4408, 0.890317, -0.5952, -0.0801653, -0.0888 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1852)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { 0.641, 0.989, 0.998, -0.005 };
   int lda = 2;
   float B[] = { -0.168, 0.465, 0.36, 0.356, -0.858, 0.879 };
   int ldb = 2;
   float B_expected[] = { 0.188365, -0.1395, -0.0023748, -0.1068, 0.518199, -0.2637 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1853)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { -0.638, 0.389, 0.997, 0.909, -0.598, -0.43, -0.345, -0.897, 0.119 };
   int lda = 3;
   float B[] = { 0.64, 0.779, -0.129, 0.016, 0.599, -0.668 };
   int ldb = 3;
   float B_expected[] = { 0.904844, 0.156956, 0.32521, 2.08405, -0.910426, 1.68403 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1854)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { 0.289, 0.641, -0.876, -0.503, -0.062, -0.987, 0.1, -0.105, 0.757 };
   int lda = 3;
   float B[] = { -0.285, 0.285, 0.219, -0.986, -0, -0.605 };
   int ldb = 3;
   float B_expected[] = { 0.124319, -0.150346, -0.0657, 0.339965, 0.17914, 0.1815 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1855)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { 0.524, 0.018, 0.292, -0.573, 0.866, 0.749, 0.99, 0.101, 0.871 };
   int lda = 3;
   float B[] = { 0.522, -0.269, -0.142, -0.266, -0.505, -0.55 };
   int ldb = 3;
   float B_expected[] = { -0.298855, -0.104554, 0.400719, 0.15229, 0.275707, -0.0156298 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1856)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { 0.283, 0.62, -0.387, -0.739, -0.599, 0.114, 0.552, 0.083, -0.976 };
   int lda = 3;
   float B[] = { 0.202, 0.169, 0.7, 0.473, 0.86, -0.557 };
   int ldb = 3;
   float B_expected[] = { -0.0606, -0.0954834, -0.168624, -0.1419, -0.362864, 0.275547 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1857)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { -0.185, 0.178, -0.22, -0.645, -0.585, -0.342, -0.594, -0.141, 0.944 };
   int lda = 3;
   float B[] = { 0.22, -0.895, -0.301, -0.683, -0.009, -0.451 };
   int ldb = 2;
   float B_expected[] = { 0.888147, -0.569939, -0.155048, -0.384802, 0.00286017, 0.143326 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1858)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { -0.145, 0.746, 0.541, 0.584, -0.394, 0.371, -0.172, -0.601, 0.542 };
   int lda = 3;
   float B[] = { 0.529, 0.636, 0.668, 0.848, -0.816, -0.925 };
   int ldb = 2;
   float B_expected[] = { -0.0854817, -0.0918985, -0.0532752, -0.0876225, 0.2448, 0.2775 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1859)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { 0.416, -0.526, -0.486, -0.716, 0.361, 0.365, -0.492, 0.544, 0.721 };
   int lda = 3;
   float B[] = { 0.25, 0.746, 0.55, 0.836, -0.024, 0.226 };
   int ldb = 2;
   float B_expected[] = { -0.180288, -0.537981, -0.719755, -1.47861, 0.25283, 0.291864 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1860)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3;
   float A[] = { -0.735, -0.606, -0.124, 0.641, -0.074, -0.053, -0.734, 0.907, 0.558 };
   int lda = 3;
   float B[] = { 0.623, 0.392, -0.808, -0.022, -0.665, -0.616 };
   int ldb = 2;
   float B_expected[] = { -0.1869, -0.1176, 0.129139, -0.0646656, 0.183169, 0.16679 };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1861)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { -0.584, -0.058, -0.964, -0.214 };
   int lda = 2;
   double B[] = { 0.073, -0.734, -0.058, -0.115, 0.513, 0.503 };
   int ldb = 3;
   double B_expected[] = { -0.0178370247087, 0.149492702599, 0.0332751888363, 0.053738317757, -0.239719626168, -0.235046728972 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1862)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { 0.251, -0.8, 0.365, 0.809 };
   int lda = 2;
   double B[] = { -0.632, -0.611, 0.9, 0.063, -0.652, -0.841 };
   int ldb = 3;
   double B_expected[] = { -0.05816, -0.11326, 0.02272, 0.0063, -0.0652, -0.0841 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1863)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { -0.833, 0.934, -0.608, 0.49 };
   int lda = 2;
   double B[] = { 0.336, -0.541, -0.729, -0.382, 0.741, 0.546 };
   int ldb = 3;
   double B_expected[] = { -0.0403361344538, 0.0649459783914, 0.0875150060024, -0.128008917853, 0.231810520126, 0.220018619693 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1864)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { 0.824, 0.907, 0.632, -0.348 };
   int lda = 2;
   double B[] = { 0.351, -0.301, 0.602, 0.873, 0.031, -0.2 };
   int ldb = 3;
   double B_expected[] = { 0.0351, -0.0301, 0.0602, 0.0651168, 0.0221232, -0.0580464 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1865)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { 0.427, 0.193, -0.959, -0.679 };
   int lda = 2;
   double B[] = { -0.646, 0.741, -0.339, 0.049, 0.734, -0.182 };
   int ldb = 2;
   double B_expected[] = { -0.3963857167, -0.10913107511, -0.0955986383061, -0.00721649484536, 0.232096380888, 0.0268041237113 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1866)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { 0.282, 0.766, -0.422, -0.518 };
   int lda = 2;
   double B[] = { 0.269, 0.211, -0.911, -0.685, -0.777, -0.919 };
   int ldb = 2;
   double B_expected[] = { 0.0358042, 0.0211, -0.120007, -0.0685, -0.1164818, -0.0919 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1867)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { -0.877, -0.818, 0.191, 0.468 };
   int lda = 2;
   double B[] = { 0.517, 0.669, 0.337, -0.579, 0.885, -0.677 };
   int ldb = 2;
   double B_expected[] = { -0.0589509692132, 0.039910485435, -0.0384264538198, -0.190882135095, -0.100912200684, -0.321038846495 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1868)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { 0.469, 0.115, 0.284, 0.139 };
   int lda = 2;
   double B[] = { 0.889, -0.002, -0.686, -0.256, 0.028, 0.371 };
   int ldb = 2;
   double B_expected[] = { 0.0889, -0.0104235, -0.0686, -0.017711, 0.0028, 0.036778 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1869)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { -0.218, -0.819, -0.523, 0.042, 0.545, -0.292, 0.283, 0.224, 0.247 };
   int lda = 3;
   double B[] = { 0.677, 0.153, -0.272, -0.226, 0.987, -0.216 };
   int ldb = 3;
   double B_expected[] = { -0.310550458716, -0.438607019611, -1.28619894589, 0.103669724771, 0.336890834105, 0.530329512606 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1870)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { 0.241, 0.561, 0.164, 0.486, 0.891, -0.508, -0.596, -0.074, 0.576 };
   int lda = 3;
   double B[] = { -0.325, 0.382, 0.368, 0.761, -0.349, 0.324 };
   int ldb = 3;
   double B_expected[] = { -0.0325, 0.0564325, 0.07079771, 0.0761, -0.0775921, -0.0194971868 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1871)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { 0.76, 0.58, -0.203, 0.053, 0.792, 0.355, -0.685, 0.449, -0.367 };
   int lda = 3;
   double B[] = { 0.861, -0.44, 0.842, -0.019, -0.382, -0.579 };
   int ldb = 3;
   double B_expected[] = { -0.0986936127734, 0.0745114634079, -0.229427792916, 0.149297547123, -0.137672708006, 0.157765667575 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1872)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { 0.802, 0.298, 0.159, 0.333, 0.515, 0.715, -0.32, -0.217, 0.301 };
   int lda = 3;
   double B[] = { -0.268, 0.1, -0.631, 0.472, 0.796, 0.278 };
   int ldb = 3;
   double B_expected[] = { -0.0457623309, -0.0036927, -0.0631, 0.0275803442, 0.0856326, 0.0278 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1873)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { -0.028, 0.186, -0.435, -0.747, 0.212, 0.257, 0.804, -0.595, 0.64 };
   int lda = 3;
   double B[] = { 0.729, -0.847, -0.577, 0.056, -0.493, 0.619 };
   int ldb = 2;
   double B_expected[] = { -2.60357142857, 3.025, -9.44607479784, 10.685259434, -5.58819230648, 6.23051463001 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1874)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { -0.74, -0.091, 0.484, 0.769, 0.91, 0.817, -0.26, 0.579, 0.393 };
   int lda = 3;
   double B[] = { 0.109, 0.969, -0.668, 0.544, 0.753, 0.796 };
   int ldb = 2;
   double B_expected[] = { 0.0109, 0.0969, -0.0751821, -0.0201161, 0.1216644359, 0.1164412219 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1875)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { 0.123, -0.328, -0.482, 0.083, -0.125, -0.712, -0.757, -0.009, 0.237 };
   int lda = 3;
   double B[] = { -0.18, 0.358, 0.839, -0.725, 0.73, -0.095 };
   int ldb = 2;
   double B_expected[] = { -5.40775366883, 2.28950005146, -2.42566413502, 0.808320675105, 0.308016877637, -0.0400843881857 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1876)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { 0.255, -0.069, -0.137, -0.45, -0.24, 0.221, -0.509, -0.484, -0.131 };
   int lda = 3;
   double B[] = { -0.563, 0.993, 0.508, 0.771, 0.745, 0.233 };
   int ldb = 2;
   double B_expected[] = { -0.0437243505, 0.1074566983, 0.0343355, 0.0719507, 0.0745, 0.0233 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1877)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { -0.772, 0.079, -0.227, 0.998 };
   int lda = 2;
   double B[] = { -0.095, 0.012, -0.988, -0.722, 0.738, 0.05 };
   int ldb = 3;
   double B_expected[] = { -0.123056994819, 0.0155440414508, -1.27979274611, 0.733187878347, -0.740709398071, 0.051206039021 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1878)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { -0.045, 0.059, -0.61, -0.328 };
   int lda = 2;
   double B[] = { 0.302, -0.099, 0.521, 0.487, -0.961, 0.903 };
   int ldb = 3;
   double B_expected[] = { -0.302, 0.099, -0.521, -0.469182, 0.955159, -0.872261 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1879)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { -0.319, 0.642, 0.511, 0.762 };
   int lda = 2;
   double B[] = { 0.883, 0.987, 0.436, -0.783, 0.175, -0.973 };
   int ldb = 3;
   double B_expected[] = { 4.41405227952, 2.72615785879, 3.41221747752, 1.02755905512, -0.229658792651, 1.27690288714 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1880)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { 0.676, 0.038, 0.543, 0.296 };
   int lda = 2;
   double B[] = { 0.804, -0.28, -0.318, 0.382, -0.165, -0.007 };
   int ldb = 3;
   double B_expected[] = { -0.596574, 0.190405, 0.314199, -0.382, 0.165, 0.007 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1881)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { -0.722, -0.355, -0.14, -0.146 };
   int lda = 2;
   double B[] = { -0.44, 0.751, -0.995, 0.625, 0.16, -0.127 };
   int ldb = 2;
   double B_expected[] = { -0.609418282548, 5.72820931203, -1.37811634349, 5.60230334307, 0.221606648199, -1.08236253937 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1882)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { 0.817, -0.619, 0.548, 0.064 };
   int lda = 2;
   double B[] = { -0.756, -0.169, 0.429, -0.789, 0.79, 0.479 };
   int ldb = 2;
   double B_expected[] = { 0.756, -0.245288, -0.429, 1.024092, -0.79, -0.04608 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1883)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { 0.496, -0.734, -0.679, -0.697 };
   int lda = 2;
   double B[] = { -0.483, -0.508, -0.819, 0.237, 0.852, -0.512 };
   int ldb = 2;
   double B_expected[] = { -0.104772180312, -0.728837876614, 2.1543973018, 0.340028694405, -2.80479705651, -0.734576757532 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1884)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { 0.049, 0.079, -0.8, -0.762 };
   int lda = 2;
   double B[] = { 0.426, 0.094, 0.794, -0.098, 0.442, -0.991 };
   int ldb = 2;
   double B_expected[] = { -0.418574, -0.094, -0.801742, 0.098, -0.520289, 0.991 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1885)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { -0.974, 0.848, -0.765, 0.528, -0.693, 0.252, -0.135, -0.507, 0.954 };
   int lda = 3;
   double B[] = { 0.395, 0.791, -0.787, 0.636, 0.271, -0.905 };
   int ldb = 3;
   double B_expected[] = { 1.01254427581, 1.4413950829, 0.824947589099, 0.548697105717, 0.736012415258, 0.948637316562 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1886)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { 0.919, -0.513, -0.38, 0.587, -0.862, 0.598, 0.714, 0.726, 0.491 };
   int lda = 3;
   double B[] = { -0.056, -0.802, -0.962, 0.656, -0.195, -0.679 };
   int ldb = 3;
   double B_expected[] = { 0.537869412, 0.226724, 0.962, -0.506244546, -0.211042, 0.679 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1887)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { -0.287, -0.009, -0.989, -0.062, 0.714, -0.293, -0.875, 0.371, 0.728 };
   int lda = 3;
   double B[] = { -0.14, -0.969, 0.702, -0.317, -0.739, -0.518 };
   int ldb = 3;
   double B_expected[] = { -0.487804878049, 1.31478445037, -2.22062403761, -1.10452961672, 0.939102470256, -1.09460224052 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1888)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { 0.236, 0.605, 0.338, -0.926, 0.362, 0.562, -0.554, 0.076, 0.85 };
   int lda = 3;
   double B[] = { 0.113, 0.604, 0.859, 0.216, -0.6, -0.048 };
   int ldb = 3;
   double B_expected[] = { -0.113, -0.708638, -0.867745512, -0.216, 0.399984, -0.102062784 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1889)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { -0.476, 0.428, -0.214, -0.889, -0.526, -0.704, 0.458, -0.479, 0.077 };
   int lda = 3;
   double B[] = { 0.124, -0.007, 0.452, 0.966, 0.42, 0.369 };
   int ldb = 2;
   double B_expected[] = { -15.8695808776, -16.2060574143, 5.8264777048, 6.20050861686, -5.45454545455, -4.79220779221 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1890)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { 0.748, 0.242, -0.964, 0.422, -0.78, -0.595, -0.926, -0.474, 0.947 };
   int lda = 3;
   double B[] = { 0.242, -0.553, -0.899, -0.714, -0.084, -0.609 };
   int ldb = 2;
   double B_expected[] = { -0.560396352, 0.693808948, 0.938816, 1.002666, 0.084, 0.609 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1891)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { 0.808, -0.074, 0.359, -0.172, -0.934, -0.67, 0.92, -0.617, -0.383 };
   int lda = 3;
   double B[] = { 0.079, 0.978, 0.82, 0.444, -0.597, -0.64 };
   int ldb = 2;
   double B_expected[] = { -0.0977722772277, -1.2103960396, 0.885690737168, 0.571273347892, -3.19977295412, -3.80492250994 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1892)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { 0.786, 0.922, -0.763, 0.498, -0.082, 0.538, 0.742, -0.391, -0.255 };
   int lda = 3;
   double B[] = { 0.911, 0.066, 0.895, 0.255, -0.547, -0.805 };
   int ldb = 2;
   double B_expected[] = { -0.911, -0.066, -0.055058, -0.194148, -0.118471796, 0.859093624 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1893)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 1};
   float A[] = { -0.362, -0.457, -0.347, -0.203, -0.517, 0.462, 0.572, 0.521 };
   int lda = 2;
   float B[] = { 0.118, -0.593, 0.773, 0.053, -0.419, -0.096, 0.846, -0.311, -0.364, 0.161, -0.496, -0.393 };
   int ldb = 3;
   float B_expected[] = { -1.58885, 0.58489, -0.628497, -0.878921, 0.701485, 1.08554, 1.03347, 0.537701, -0.470639, -0.207688, -0.056162, -0.815978 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1894) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1894) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 1};
   float A[] = { 0.845, -0.654, -0.43, -0.834, 0.206, 0.414, 0.761, 0.961 };
   int lda = 2;
   float B[] = { 0.069, 0.005, -0.419, 0.806, 0.857, 0.669, 0.942, 0.657, 0.52, 0.19, -0.609, -0.305 };
   int ldb = 3;
   float B_expected[] = { -1.07314, -0.073878, -1.32138, -0.35386, -0.029944, 0.8495, -0.657, 0.942, -0.19, 0.52, 0.305, -0.609 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1895) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1895) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 1};
   float A[] = { -0.87, 0.218, 0.813, 0.575, -0.848, 0.7, -0.311, 0.374 };
   int lda = 2;
   float B[] = { 0.117, 0.758, -0.189, -0.768, 0.857, -0.269, 0.796, -0.592, -0.499, 0.977, 0.643, 0.282 };
   int ldb = 3;
   float B_expected[] = { 0.851499, 0.0788813, -0.881826, -0.00372192, -0.0586805, -0.999761, -1.37808, -2.51525, 2.45259, 2.57925, 1.0972, 1.8459 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1896) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1896) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 1};
   float A[] = { -0.697, -0.668, 0.746, -0.818, 0.651, 0.275, -0.702, -0.615 };
   int lda = 2;
   float B[] = { -0.876, 0.842, -0.848, 0.901, 0.75, 0.361, -0.702, 0.039, -0.41, 0.541, 0.489, 0.025 };
   int ldb = 3;
   float B_expected[] = { -0.842, -0.876, -0.901, -0.848, -0.361, 0.75, 0.268242, 0.099826, -0.187649, 0.389823, 0.416261, 0.100025 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1897) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1897) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 1};
   float A[] = { -0.394, -0.65, -0.082, -0.632, -0.53, 0.483, 0.149, -0.192 };
   int lda = 2;
   float B[] = { -0.691, 0.732, 0.976, 0.073, 0.607, 0.918, -0.918, 0.67, 0.37, -0.344, -0.114, -0.62 };
   int ldb = 2;
   float B_expected[] = { -1.39367, -3.05481, -3.35679, 2.2248, 4.33836, -1.06673, 1.29393, -4.49373, -1.89826, 2.23995, 1.93461, 1.72783 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1898) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1898) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 1};
   float A[] = { 0.45, 0.972, -0.051, 0.71, -0.127, -0.274, 0.152, 0.789 };
   int lda = 2;
   float B[] = { 0.683, -0.915, -0.773, 0.088, -0.28, 0.17, 0.818, 0.293, -0.551, 0.365, 0.899, 0.257 };
   int ldb = 2;
   float B_expected[] = { 1.11563, 0.560717, -0.088, -0.773, -0.431343, -0.256396, -0.293, 0.818, -0.643965, -0.507245, -0.257, 0.899 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1899) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1899) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 1};
   float A[] = { -0.993, -0.028, -0.547, -0.251, 0.781, -0.315, 0.865, 0.229 };
   int lda = 2;
   float B[] = { 0.578, 0.73, -0.931, 0.288, 0.048, 0.508, -0.168, 0.655, 0.92, -0.26, 0.485, 0.05 };
   int ldb = 2;
   float B_expected[] = { 0.718162, -0.602325, -0.0323644, -1.24023, 0.509813, -0.0627138, -0.410611, 0.0227614, -0.287729, -0.918372, -0.000635986, -0.10338 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1900) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1900) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 1};
   float A[] = { 0.131, -0.494, 0.615, -0.089, -0.22, -0.874, 0.677, 0.074 };
   int lda = 2;
   float B[] = { -0.276, 0.539, 0.647, 0.986, -0.34, 0.983, -0.819, 0.144, 0.361, 0.561, 0.178, -0.433 };
   int ldb = 2;
   float B_expected[] = { -0.539, -0.276, -0.629951, 0.768769, -0.983, -0.34, 0.490805, -0.697387, -0.561, 0.361, 0.745886, -0.093944 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1901) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1901) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 1};
   float A[] = { -0.123, -0.726, -0.011, 0.245, -0.205, 0.77, -0.81, -0.973, 0.354, -0.835, 0.552, 0.396, -0.524, -0.204, -0.814, 0.284, -0.976, -0.835 };
   int lda = 3;
   float B[] = { -0.42, 0.976, -0.845, 0.651, -0.44, -0.862, 0.137, 0.066, -0.63, 0.482, -0.187, 0.724 };
   int ldb = 3;
   float B_expected[] = { 0.783777, -1.21156, 0.66205, -1.40548, 0.886226, 0.0391664, -0.168468, -0.119451, 0.378144, -0.774828, 0.708857, -0.807468 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1902) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1902) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 1};
   float A[] = { 0.921, 0.167, -0.41, 0.578, -0.372, 0.106, 0.551, 0.668, 0.295, 0.855, -0.167, 0.976, -0.782, -0.777, 0.278, -0.98, 0.038, -0.832 };
   int lda = 3;
   float B[] = { 0.459, 0.06, 0.387, 0.871, -0.366, 0.926, 0.236, -0.889, 0.619, 0.319, -0.709, 0.884 };
   int ldb = 3;
   float B_expected[] = { -0.06, 0.459, -0.630298, 0.60987, -0.409693, 0.528127, 0.889, 0.236, 0.181898, 0.201918, -0.300827, -0.859254 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1903) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1903) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 1};
   float A[] = { -0.691, -0.817, 0.954, -0.969, -0.574, -0.026, 0.992, 0.529, 0.135, -0.413, -0.314, -0.859, -0.284, -0.849, 0.781, 0.534, -0.018, 0.282 };
   int lda = 3;
   float B[] = { -0.028, -0.429, 0.066, -0.854, -0.316, 0.514, -0.465, -0.857, 0.286, 0.415, -0.486, 0.538 };
   int ldb = 3;
   float B_expected[] = { 6.83575, 2.7232, 3.79999, 5.15624, -1.00015, 1.88653, 5.42614, 2.69261, 2.30584, 3.85628, -1.59513, 2.00962 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1904) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1904) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 1};
   float A[] = { -0.839, -0.318, 0.175, 0.72, -0.683, 0.395, -0.279, 0.151, -0.71, 0.445, 0.533, -0.38, -0.749, -0.833, 0.871, -0.426, 0.195, 0.889 };
   int lda = 3;
   float B[] = { 0.804, -0.346, 0.234, 0.782, 0.033, 0.581, 0.981, -0.68, 0.919, -0.758, 0.152, -0.503 };
   int ldb = 3;
   float B_expected[] = { -0.20395, 0.376748, -0.290007, -0.042249, -0.581, 0.033, 1.15245, 1.75457, 0.255135, 1.00089, 0.503, 0.152 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1905) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1905) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 1};
   float A[] = { -0.365, -0.662, 0.188, -0.571, 0.082, 0.192, -0.833, -0.958, 0.159, -0.203, 0.481, 0.08, -0.954, 0.681, -0.015, 0.146, -0.352, -0.068 };
   int lda = 3;
   float B[] = { 0.779, -0.691, -0.516, 0.148, 0.721, 0.217, -0.976, -0.963, 0.532, -0.366, 0.176, 0.4 };
   int ldb = 2;
   float B_expected[] = { -1.34375, 0.302916, 0.692272, 0.158126, -2.93098, -5.71682, 3.87247, 3.8052, 3.25028, -6.53201, -2.34332, 2.30748 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1906) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1906) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 1};
   float A[] = { 0.779, 0.065, -0.616, -0.245, 0.823, 0.689, 0.06, -0.164, 0.768, -0.727, 0.897, -0.556, -0.875, 0.862, 0.863, -0.085, 0.171, 0.063 };
   int lda = 3;
   float B[] = { -0.621, 0.428, 0.096, 0.711, 0.416, -0.684, 0.806, 0.491, 0.037, -0.776, -0.312, 0.391 };
   int ldb = 2;
   float B_expected[] = { -0.428, -0.621, -0.711, 0.096, 0.811524, 0.383068, -0.464084, 0.683636, -0.866708, -0.399047, -0.587978, -0.244543 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1907) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1907) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 1};
   float A[] = { -0.349, 0.224, -0.323, 0.457, 0.081, 0.443, 0.809, 0.037, -0.543, 0.554, 0.779, 0.632, -0.852, -0.148, -0.649, -0.78, 0.469, -0.515 };
   int lda = 3;
   float B[] = { 0.162, 0.754, -0.978, -0.097, 0.986, 0.943, 0.676, 0.718, 0.204, 0.264, -0.124, -0.73 };
   int ldb = 2;
   float B_expected[] = { 0.0811068, 1.92921, -3.74716, 1.18561, 1.80842, -0.638944, 0.528341, 1.20828, -0.471728, -0.083028, 0.837267, 0.654994 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1908) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1908) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 1};
   float A[] = { -0.469, -0.164, -0.792, -0.454, 0.206, 0.785, 0.504, -0.561, 0.205, 0.463, -0.8, 0.803, 0.283, 0.131, 0.576, -0.431, 0.297, -0.415 };
   int lda = 3;
   float B[] = { -0.364, 0.853, 0.056, -0.78, 0.05, 0.223, -0.166, -0.097, 0.24, 0.721, 0.023, 0.508 };
   int ldb = 2;
   float B_expected[] = { -1.3696, 0.527133, 0.554099, 0.524136, -0.60708, 0.820963, -0.290931, 0.260324, -0.721, 0.24, -0.508, 0.023 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1909) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1909) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {1, 0};
   float A[] = { 0.917, 0.367, -0.115, -0.321, -0.811, -0.563, 0.78, -0.742 };
   int lda = 2;
   float B[] = { 0.797, 0.166, 0.737, -0.685, 0.677, -0.04, -0.652, 0.327, 0.094, -0.656, 0.496, -0.646 };
   int ldb = 3;
   float B_expected[] = { 0.811592, -0.143789, 0.435059, -0.92112, 0.621302, -0.292277, -0.710488, 0.0561583, 0.694329, -0.137285, 0.752465, 0.100199 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1910) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1910) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {1, 0};
   float A[] = { -0.384, 0.783, -0.086, -0.649, -0.574, 0.216, -0.809, -0.608 };
   int lda = 2;
   float B[] = { 0.067, -0.183, -0.524, 0.77, 0.169, 0.769, -0.982, -0.522, -0.051, -0.129, 0.595, 0.56 };
   int ldb = 3;
   float B_expected[] = { 0.067, -0.183, -0.524, 0.77, 0.169, 0.769, -0.857471, -0.494255, -0.595794, -0.402856, 0.110453, 0.735815 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1911) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1911) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {1, 0};
   float A[] = { 0.228, -0.644, 0.731, 0.458, 0.051, -0.725, 0.731, 0.537 };
   int lda = 2;
   float B[] = { -0.588, 0.01, -0.009, -0.374, 0.422, 0.758, -0.428, 0.263, 0.659, 0.171, -0.239, 0.968 };
   int ldb = 3;
   float B_expected[] = { -0.232749, -1.39168, -0.124158, 0.287962, -1.55821, 0.0298572, -0.208619, 0.513035, 0.697138, -0.278198, 0.419466, 1.01607 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1912) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1912) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {1, 0};
   float A[] = { -0.321, 0.761, 0.809, -0.017, -0.009, -0.975, 0.057, 0.396 };
   int lda = 2;
   float B[] = { 0.377, 0.776, -0.686, -0.561, 0.29, 0.601, 0.755, 0.518, 0.313, -0.394, 0.945, 0.395 };
   int ldb = 3;
   float B_expected[] = { -0.121255, 1.51679, -0.299033, -0.259371, -0.08662, 1.52593, 0.755, 0.518, 0.313, -0.394, 0.945, 0.395 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1913) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1913) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {1, 0};
   float A[] = { 0.186, 0.818, -0.142, -0.376, 0.332, 0.746, 0.413, -0.151 };
   int lda = 2;
   float B[] = { -0.374, -0.787, 0.221, -0.104, 0.74, -0.548, 0.88, -0.66, 0.65, 0.046, -0.839, -0.783 };
   int ldb = 2;
   float B_expected[] = { -1.01366, 0.226724, 1.10152, 1.79962, -0.441403, -1.00501, 0.588898, 0.222456, 0.225271, -0.743398, -2.5862, -2.65075 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1914) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1914) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {1, 0};
   float A[] = { -0.574, 0.018, -0.584, -0.184, 0.41, 0.075, 0.92, 0.022 };
   int lda = 2;
   float B[] = { 0.524, -0.234, 0.198, 0.079, -0.449, -0.433, -0.14, -0.201, -0.242, -0.368, -0.298, 0.693 };
   int ldb = 2;
   float B_expected[] = { 0.524, -0.234, -0.03439, 0.13564, -0.449, -0.433, 0.011615, 0.010205, -0.242, -0.368, -0.22638, 0.86203 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1915) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1915) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {1, 0};
   float A[] = { 0.422, -0.38, 0.919, -0.229, -0.849, -0.19, 0.02, -0.181 };
   int lda = 2;
   float B[] = { 0.971, -0.339, 0.203, 0.083, 0.461, -0.623, 0.334, 0.653, 0.694, 0.42, 0.239, -0.061 };
   int ldb = 2;
   float B_expected[] = { 3.06394, -0.745692, -0.330599, 1.15808, 8.0252, -0.902398, -3.36278, 2.21688, 0.70369, -0.872941, 0.477097, 1.26772 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1916) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1916) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {1, 0};
   float A[] = { 0.564, 0.483, -0.635, 0.84, 0.238, 0.35, 0.96, 0.397 };
   int lda = 2;
   float B[] = { 0.963, -0.513, 0.989, 0.404, -0.352, 0.924, 0.052, -0.059, -0.771, 0.341, -0.566, -0.844 };
   int ldb = 2;
   float B_expected[] = { 1.93037, -1.08722, 0.989, 0.404, -0.36854, 0.842855, 0.052, -0.059, -1.83937, 0.2805, -0.566, -0.844 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1917) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1917) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {1, 0};
   float A[] = { 0.182, 0.699, 0.303, -0.273, -0.363, 0.02, 0.991, -0.206, -0.347, 0.269, -0.384, 0.797, 0.392, -0.966, 0.347, 0.87, 0.016, -0.097 };
   int lda = 3;
   float B[] = { 0.587, 0.875, -0.848, 0.154, -0.887, -0.709, 0.824, -0.895, 0.159, 0.933, -0.011, -0.393 };
   int ldb = 3;
   float B_expected[] = { -14.9753, 2.31554, 0.613295, 24.1527, 5.64728, -10.0758, -3.79783, -5.34545, -5.38045, 2.99977, 3.92602, -0.760993 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1918) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1918) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {1, 0};
   float A[] = { 0.637, -0.57, 0.322, -0.303, 0.618, 0.261, 0.654, -0.238, 0.66, -0.485, 0.223, -0.196, -0.252, 0.929, -0.012, 0.965, 0.783, 0.489 };
   int lda = 3;
   float B[] = { 0.894, 0.93, 0.648, 0.914, 0.7, -0.138, 0.63, -0.173, -0.671, -0.327, -0.922, 0.816 };
   int ldb = 3;
   float B_expected[] = { -0.0695574, 0.64143, 0.518948, 1.08197, 0.7, -0.138, 1.8231, -0.404044, -0.62533, -0.68968, -0.922, 0.816 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1919) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1919) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {1, 0};
   float A[] = { 0.274, 0.721, -0.445, 0.14, 0.023, -0.945, 0.859, -0.522, -0.227, 0.722, 0.165, 0.969, -0.212, -0.816, 0.908, -0.652, -0.208, -0.229 };
   int lda = 3;
   float B[] = { 0.011, -0.818, 0.067, -0.191, -0.911, 0.84, -0.162, -0.951, -0.502, -0.21, 0.492, 0.767 };
   int ldb = 3;
   float B_expected[] = { -0.986296, -0.390076, -0.910328, -1.26205, -3.05033, 0.930902, -1.22716, -0.241667, -1.07925, -0.600129, -2.84941, 5.27338 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1920) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1920) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {1, 0};
   float A[] = { -0.186, 0.118, -0.545, 0.784, 0.057, 0.39, 0.77, -0.518, -0.97, 0.271, 0.488, 0.637, -0.482, -0.993, -0.797, -0.945, 0.257, 0.3 };
   int lda = 3;
   float B[] = { -0.783, 0.649, 0.698, 0.046, -0.153, 0.473, -0.996, -0.211, 0.84, 0.201, -0.457, 0.918 };
   int ldb = 3;
   float B_expected[] = { -0.783, 0.649, 0.964728, -0.859324, 0.406086, 0.235086, -0.996, -0.211, 1.71622, -0.152458, 0.78435, 1.32759 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1921) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1921) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {1, 0};
   float A[] = { -0.681, -0.342, -0.195, -0.053, 0.016, -0.191, 0.989, -0.718, -0.59, 0.646, -0.41, -0.809, -0.359, -0.783, -0.902, 0.917, -0.703, 0.795 };
   int lda = 3;
   float B[] = { -0.27, 0.037, 0.349, 0.36, -0.293, 0.128, -0.481, -0.834, -0.815, -0.6, 0.728, 0.122 };
   int ldb = 2;
   float B_expected[] = { -0.69977, -2.39368, 2.17354, 1.74016, 0.260417, -1.25151, 0.175881, 1.93577, 0.085191, 0.949825, -0.368302, -0.590043 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1922) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1922) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {1, 0};
   float A[] = { -0.132, 0.33, 0.357, 0.32, 0.833, -0.111, -0.192, -0.643, -0.622, -0.663, -0.58, 0.423, -0.874, 0.86, -0.281, -0.992, 0.055, 0.137 };
   int lda = 3;
   float B[] = { 0.104, -0.906, -0.712, 0.103, -0.474, -0.591, 0.073, -0.906, -0.261, -0.391, 0.881, -0.345 };
   int ldb = 2;
   float B_expected[] = { 0.126148, -1.31009, -0.0285057, -0.554776, -0.159469, -0.959783, 0.662801, -0.128993, -0.261, -0.391, 0.881, -0.345 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1923) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1923) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {1, 0};
   float A[] = { -0.353, -0.581, -0.648, 0.894, 0.825, -0.23, -0.529, 0.213, 0.568, 0.296, 0.372, 0.442, 0.515, -0.409, 0.222, -0.246, -0.524, 0.318 };
   int lda = 3;
   float B[] = { -0.467, 0.632, 0.672, 0.777, -0.609, 0.511, -0.991, 0.311, -0.617, -0.732, -0.585, 0.152 };
   int ldb = 2;
   float B_expected[] = { -0.437806, -1.06979, -1.49004, 0.251317, -2.40924, 1.62379, -1.09482, 3.75003, -1.80514, -2.07012, -4.8059, -0.418185 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1924) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1924) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {1, 0};
   float A[] = { -0.201, -0.918, -0.514, -0.889, -0.948, 0.34, 0.818, 0.557, 0.341, 0.484, 0.235, 0.561, 0.874, -0.342, -0.411, -0.975, -0.85, -0.621 };
   int lda = 3;
   float B[] = { -0.389, -0.252, 0.322, -0.763, -0.839, -0.744, -0.946, -0.312, 0.051, -0.686, -0.626, -0.043 };
   int ldb = 2;
   float B_expected[] = { -0.389, -0.252, 0.322, -0.763, -0.814918, -1.21935, -0.102185, -0.417924, -0.896001, -0.04892, -0.790606, -0.720266 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1925) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1925) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3, 0.1};
   float A[] = { 0.571, 0.12, 0.665, 0.425, -0.977, -0.772, -0.944, -0.154 };
   int lda = 2;
   float B[] = { -0.357, -0.213, 0.57, 0.134, 0.089, 0.046, 0.027, 0.825, -0.127, 0.658, -0.332, 0.247 };
   int ldb = 3;
   float B_expected[] = { 0.205417, 0.092557, -0.315204, -0.0368205, -0.0507703, -0.0192512, 0.238158, 0.270895, -0.257649, 0.296502, -0.140106, 0.100105 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1926) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1926) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3, 0.1};
   float A[] = { 0.051, 0.966, 0.04, -0.765, 0.276, -0.798, 0.766, -0.37 };
   int lda = 2;
   float B[] = { 0.532, 0.59, 0.305, 0.443, 0.036, 0.655, -0.145, -0.864, -0.483, -0.45, -0.327, -0.365 };
   int ldb = 3;
   float B_expected[] = { -0.2186, -0.1238, -0.1358, -0.1024, -0.0763, -0.1929, 0.043937, 0.416881, 0.116996, 0.194683, -0.0099165, 0.142885 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1927) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1927) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3, 0.1};
   float A[] = { -0.163, -0.238, -0.032, 0.494, 0.863, 0.96, 0.669, 0.415 };
   int lda = 2;
   float B[] = { -0.724, -0.682, 0.034, 0.352, 0.42, 0.253, 0.186, -0.061, 0.278, -0.764, -0.484, 0.051 };
   int ldb = 3;
   float B_expected[] = { -0.532386, -1.09223, -1.1606, 1.43429, 1.04476, 0.724237, -0.0783541, 0.00655162, -0.179639, 0.27272, 0.193877, 0.0250509 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1928) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1928) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3, 0.1};
   float A[] = { 0.533, 0.575, 0.808, -0.631, 0.185, 0.296, -0.757, -0.279 };
   int lda = 2;
   float B[] = { -0.744, -0.881, -0.594, 0.629, -0.924, 0.017, -0.089, -0.052, 0.959, -0.486, 0.39, -0.378 };
   int ldb = 3;
   float B_expected[] = { 0.303415, 0.198103, 0.0879903, -0.363588, 0.245042, -0.149137, 0.0319, 0.0067, -0.2391, 0.2417, -0.0792, 0.1524 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1929) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1929) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3, 0.1};
   float A[] = { -0.001, -0.948, -0.97, -0.285, -0.664, -0.977, -0.746, 0.192 };
   int lda = 2;
   float B[] = { 0.997, -0.852, 0.87, -0.955, 0.007, -0.071, -0.263, -0.077, -0.856, 0.228, -0.81, 0.476 };
   int ldb = 2;
   float B_expected[] = { 0.375027, 0.225237, -0.432345, -0.0987217, 0.0232012, -0.00529874, -0.112225, 0.0682749, -0.162707, -0.246664, 0.267117, 0.237712 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1930) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1930) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3, 0.1};
   float A[] = { 0.804, 0.476, -0.898, -0.966, 0.51, -0.346, 0.622, -0.749 };
   int lda = 2;
   float B[] = { -0.964, 0.453, 0.799, -0.949, -0.055, 0.803, 0.99, -0.162, 0.913, -0.081, -0.057, 0.014 };
   int ldb = 2;
   float B_expected[] = { 0.2439, -0.2323, -0.349565, 0.398684, -0.0638, -0.2464, -0.333516, 0.295339, -0.2658, 0.1156, 0.191256, 0.0231108 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1931) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1931) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3, 0.1};
   float A[] = { -0.144, -0.611, -0.054, 0.618, 0.213, 0.49, -0.465, -0.488 };
   int lda = 2;
   float B[] = { -0.225, -0.663, 0.073, -0.379, -0.297, 0.822, -0.038, -0.935, -0.81, 0.885, -0.065, 0.412 };
   int ldb = 2;
   float B_expected[] = { 0.287563, -0.439427, 0.113582, -0.141015, -0.375321, -0.339988, 0.189826, -0.395838, -0.655583, 0.0702722, -0.117522, 0.15645 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1932) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1932) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3, 0.1};
   float A[] = { -0.367, 0.502, -0.309, 0.404, 0.531, -0.188, 0.181, 0.583 };
   int lda = 2;
   float B[] = { 0.861, -0.648, 0.906, -0.402, 0.455, 0.412, 0.34, -0.248, 0.107, 0.507, 0.088, -0.593 };
   int ldb = 2;
   float B_expected[] = { -0.350389, 0.252194, -0.2316, 0.2112, -0.245348, -0.0757932, -0.0772, 0.1084, -0.148061, -0.0704181, 0.0329, 0.1867 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1933) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1933) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3, 0.1};
   float A[] = { -0.476, 0.202, -0.66, 0.774, -0.943, -0.99, -0.035, 0.901, -0.742, -0.085, -0.335, -0.591, 0.799, 0.515, 0.753, 0.76, -0.042, -0.011 };
   int lda = 3;
   float B[] = { 0.025, -0.976, -0.44, 0.741, -0.126, 0.527, 0.743, 0.216, 0.661, -0.071, 0.564, -0.093 };
   int ldb = 3;
   float B_expected[] = { -6.73789, 0.501263, -2.62173, -2.22684, -0.664138, 3.89034, 4.11106, 5.79368, -1.20958, 3.39994, 4.05469, -0.945199 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1934) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1934) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3, 0.1};
   float A[] = { -0.999, 0.418, 0.687, 0.6, 0.106, -0.737, -0.165, 0.263, 0.998, -0.092, 0.555, -0.671, -0.162, -0.814, 0.317, 0.582, 0.302, -0.48 };
   int lda = 3;
   float B[] = { 0.699, 0.128, 0.296, -0.021, 0.654, 0.14, 0.008, 0.94, -0.963, 0.333, -0.481, -0.917 };
   int ldb = 3;
   float B_expected[] = { -0.312717, 0.0986958, 0.0456624, 0.163957, -0.2102, 0.0234, 0.143952, 0.0170999, 0.276937, -0.480541, 0.236, 0.227 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1935) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1935) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3, 0.1};
   float A[] = { -0.932, 0.532, -0.763, -0.029, -0.524, -0.938, 0.007, -0.445, -0.659, 0.709, -0.581, 0.825, -0.904, -0.453, 0.119, 0.964, -0.649, 0.48 };
   int lda = 3;
   float B[] = { -0.571, 0.138, 0.038, -0.175, 0.737, 0.567, -0.569, 0.062, 0.522, -0.625, 0.156, 0.799 };
   int ldb = 3;
   float B_expected[] = { -0.0819591, 0.15247, -0.121808, -0.00810757, 0.287388, -0.154159, -0.0982488, 0.13709, -0.190946, -0.223188, 0.0729118, 0.274542 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1936) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1936) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3, 0.1};
   float A[] = { 0.527, 0.272, 0.431, 0.642, -0.239, -0.254, -0.231, 0.766, 0.85, -0.09, 0.679, -0.898, 0.192, -0.651, -0.869, 0.859, 0.68, 0.03 };
   int lda = 3;
   float B[] = { 0.867, 0.816, -0.643, 0.509, -0.594, -0.833, -0.174, 0.51, 0.676, 0.115, 0.261, -0.409 };
   int ldb = 3;
   float B_expected[] = { -0.3417, -0.1581, 0.184172, -0.515263, 0.82684, 0.153742, 0.0012, -0.1704, -0.0834964, -0.0053432, -0.216529, 0.104369 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1937) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1937) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3, 0.1};
   float A[] = { 0.606, -0.757, 0.503, -0.649, -0.269, -0.484, 0.626, -0.107, -0.867, -0.047, -0.779, 0.675, 0.249, 0.645, -0.755, 0.242, 0.941, 0.189 };
   int lda = 3;
   float B[] = { -0.402, 0.252, -0.214, 0.745, 0.342, -0.98, -0.096, 0.38, -0.543, 0.605, 0.63, -0.059 };
   int ldb = 2;
   float B_expected[] = { 0.349049, -0.0955741, -0.472341, -0.259287, -0.176304, -0.239347, 0.191174, 0.170679, 0.152979, -0.219859, -0.203592, 0.0448683 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1938) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1938) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3, 0.1};
   float A[] = { -0.797, -0.288, 0.943, -0.821, -0.565, 0.73, -0.146, -0.967, 0.473, -0.095, 0.877, 0.178, -0.159, 0.021, -0.988, 0.296, 0.279, -0.513 };
   int lda = 3;
   float B[] = { -0.455, 0.859, -0.21, 0.702, -0.591, -0.235, 0.519, 0.279, -0.444, 0.816, -0.507, 0.893 };
   int ldb = 2;
   float B_expected[] = { -0.136371, -0.712172, -0.311667, -0.302476, 0.337384, -0.259056, -0.027248, -0.327988, 0.0516, -0.2892, 0.0628, -0.3186 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1939) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1939) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3, 0.1};
   float A[] = { 0.791, 0.19, -0.549, 0.994, -0.822, 0.679, -0.586, 0.042, -0.159, -0.86, 0.065, 0.943, -0.545, 0.403, 0.199, 0.76, 0.159, 0.715 };
   int lda = 3;
   float B[] = { -0.336, 0.317, 0.502, 0.543, 0.027, 0.802, 0.391, 0.716, -0.154, 0.436, 0.738, -0.029 };
   int ldb = 2;
   float B_expected[] = { 0.119543, -0.133991, -0.212552, -0.193533, -0.239565, -0.0842153, -0.531028, 0.229828, 0.61223, 0.265016, 0.850081, -0.810046 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1940) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1940) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3, 0.1};
   float A[] = { -0.182, -0.821, -0.756, -0.479, -0.191, -0.989, -0.466, 0.018, 0.85, 0.516, -0.826, 0.209, -0.321, -0.988, -0.936, -0.745, -0.57, -0.362 };
   int lda = 3;
   float B[] = { -0.501, 0.915, -0.928, 0.722, -0.542, -0.828, -0.875, -0.981, 0.425, 0.347, -0.929, -0.596 };
   int ldb = 2;
   float B_expected[] = { 0.0588, -0.3246, 0.2062, -0.3094, 0.134369, -0.0793628, 0.368285, -0.125876, -0.344423, -0.219222, 0.402199, -0.204129 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1941) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1941) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { 0.976, -0.41, -0.313, -0.779, -0.164, 0.571, 0.056, -0.526 };
   int lda = 2;
   double B[] = { -0.177, 0.837, 0.391, -0.853, -0.633, 0.693, -0.392, -0.356, -0.708, 0.926, -0.093, -0.337 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1942) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1942) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { -0.383, 0.141, 0.889, -0.007, -0.148, -0.068, 0.481, 0.675 };
   int lda = 2;
   double B[] = { 0.469, 0.735, -0.47, -0.164, 0.994, -0.483, -0.354, 0.357, 0.51, 0.523, 0.934, -0.592 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1943) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1943) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { -0.089, -0.391, -0.317, -0.349, 0.618, -0.541, -0.84, 0.31 };
   int lda = 2;
   double B[] = { 0.931, -0.257, -0.048, 0.633, -0.32, -0.576, -0.682, 0.953, -0.412, 0.408, -0.809, 0.092 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1944) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1944) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { 0.599, -0.01, -0.045, 0.567, 0.827, -0.969, -0.729, 0.538 };
   int lda = 2;
   double B[] = { 0.971, -0.626, -0.77, -0.882, 0.434, 0.269, -0.456, 0.497, 0.289, 0.957, 0.447, -0.921 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1945) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1945) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { 0.441, -0.501, 0.607, -0.4, -0.976, -0.523, -0.136, -0.492 };
   int lda = 2;
   double B[] = { 0.639, 0.872, -0.436, 0.518, 0.164, -0.04, 0.489, 0.201, 0.723, -0.958, 0.934, -0.549 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1946) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1946) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { -0.603, -0.475, 0.598, -0.666, -0.733, 0.04, 0.491, -0.592 };
   int lda = 2;
   double B[] = { 0.71, -0.827, 0.947, -0.364, 0.235, 0.294, 0.298, -0.401, -0.193, -0.008, 0.122, -0.47 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1947) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1947) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { -0.73, -0.823, 0.636, -0.965, 0.886, -0.236, 0.501, -0.301 };
   int lda = 2;
   double B[] = { 0.259, 0.701, -0.033, 0.616, -0.646, -0.177, -0.886, 0.589, -0.736, -0.303, -0.995, 0.982 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1948) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1948) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { 0.829, -0.889, 0.382, 0.083, 0.006, -0.76, -0.338, -0.601 };
   int lda = 2;
   double B[] = { 0.006, 0.381, 0.241, 0.096, -0.672, 0.664, 0.952, -0.376, -0.803, 0.344, -0.09, -0.175 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1949) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1949) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { 0.879, -0.511, -0.814, -0.94, 0.91, 0.761, 0.223, 0.03, -0.689, -0.739, -0.814, 0.463, 0.389, 0.615, -0.175, 0.129, -0.904, 0.102 };
   int lda = 3;
   double B[] = { 0.383, 0.328, 0.589, -0.29, 0.912, 0.327, 0.629, 0.883, -0.578, -0.708, 0.168, -0.982 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1950) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1950) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { -0.91, -0.182, 0.333, 0.193, 0.14, 0.538, 0.161, -0.034, -0.614, -0.154, 0.881, 0.842, 0.183, -0.229, 0.099, 0.062, -0.121, 0.179 };
   int lda = 3;
   double B[] = { -0.138, 0.109, -0.87, -0.161, 0.917, 0.443, 0.798, 0.677, -0.574, 0.327, -0.626, 0.446 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1951) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1951) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { -0.491, -0.021, -0.833, 0.921, -0.71, 0.282, 0.638, 0.223, -0.434, 0.921, -0.949, 0.457, -0.665, -0.844, -0.633, -0.874, -0.73, 0.637 };
   int lda = 3;
   double B[] = { -0.047, 0.714, 0.678, 0.756, 0.003, 0.359, 0.507, -0.197, -0.726, 0.873, -0.118, -0.996 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1952) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1952) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { 0.372, 0.354, -0.537, 0.948, -0.348, 0.808, 0.573, -0.797, 0.818, 0.701, -0.749, -0.801, -0.959, -0.781, 0.727, -0.189, 0.244, 0.414 };
   int lda = 3;
   double B[] = { 0.852, -0.714, 0.455, 0.171, -0.128, 0.554, 0.342, -0.203, 0.669, 0.619, -0.76, 0.759 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1953) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1953) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { 0.1, -0.975, 0.885, -0.608, -0.303, 0.87, -0.763, 0.409, 0.501, 0.522, -0.176, 0.679, -0.681, -0.815, -0.878, 0.86, 0.348, -0.65 };
   int lda = 3;
   double B[] = { -0.245, 0.954, -0.465, -0.931, 0.327, 0.288, -0.067, 0.252, 0.124, -0.073, -0.731, 0.176 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1954) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1954) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { 0.572, 0.045, -0.465, 0.113, 0.996, -0.597, 0.712, 0.945, 0.053, -0.436, 0.36, 0.035, -0.489, -0.012, 0.23, 0.22, 0.068, -0.586 };
   int lda = 3;
   double B[] = { -0.543, -0.809, -0.641, -0.744, 0.507, -0.742, -0.279, -0.835, -0.097, -0.968, 0.984, -0.813 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1955) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1955) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { 0.837, 0.576, -0.396, 0.013, -0.567, 0.59, 0.513, 0.824, 0.045, 0.486, 0.386, 0.766, 0.222, 0.042, 0.091, -0.008, 0.43, 0.102 };
   int lda = 3;
   double B[] = { 0.16, -0.958, -0.125, 0.833, 0.344, 0.213, 0.2, -0.689, 0.81, 0.415, -0.198, 0.001 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1956) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1956) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { -0.351, 0.7, -0.495, 0.448, -0.229, 0.925, -0.269, 0.251, -0.783, -0.223, 0.582, 0.373, -0.095, -0.383, -0.087, -0.043, -0.315, -0.999 };
   int lda = 3;
   double B[] = { -0.067, -0.104, 0.92, -0.333, 0.367, 0.995, 0.86, 0.425, 0.12, -0.756, 0.441, -0.214 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1957) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1957) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.572, -0.073, 0.878, -0.688, -0.615, -0.213, -0.643, 0.809 };
   int lda = 2;
   double B[] = { -0.973, -0.481, 0.071, -0.71, -0.669, 0.717, -0.09, -0.304, -0.427, 0.625, 0.539, -0.565 };
   int ldb = 3;
   double B_expected[] = { 0.574560994608, 0.155494672389, 0.0371747871512, 0.389534544514, 0.283820482207, -0.45678514825, 0.591891359193, 0.214411302729, -0.27258111691, 0.507180331171, 0.645135319443, -0.46315922005 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1958) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1958) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.022, 0.475, 0.444, 0.252, -0.871, 0.867, -0.093, 0.264 };
   int lda = 2;
   double B[] = { 0.696, 0.259, 0.494, 0.162, -0.9, 0.143, 0.436, 0.487, -0.733, 0.138, -0.618, 0.572 };
   int ldb = 3;
   double B_expected[] = { -0.2347, -0.0081, -0.1644, 0.0008, 0.2557, -0.1329, -0.0773344, -0.0397592, 0.2792952, -0.0736264, -0.0188216, -0.2388288 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1959) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1959) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.918, -0.459, 0.547, 0.887, 0.4, -0.497, 0.49, -0.313 };
   int lda = 2;
   double B[] = { 0.028, 0.482, -0.59, -0.533, -0.594, 0.544, -0.717, -0.524, 0.07, -0.839, 0.538, -0.548 };
   int ldb = 3;
   double B_expected[] = { -0.258092239243, -0.278373561582, 0.128448307703, -0.0949352940165, 0.35005709854, -0.355276452021, 0.308556833073, 0.371588344391, -0.148348709879, 0.433197660833, -0.356526626221, 0.217565644883 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1960) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1960) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.86, -0.532, -0.396, -0.116, 0.766, -0.818, -0.335, 0.271 };
   int lda = 2;
   double B[] = { -0.029, -0.754, -0.566, -0.108, 0.904, -0.038, 0.07, -0.476, -0.48, 0.961, 0.864, -0.593 };
   int ldb = 3;
   double B_expected[] = { -0.058812, 0.130312, 0.419002, 0.272588, -0.330474, -0.264172, 0.0266, 0.1498, 0.0479, -0.3363, -0.1999, 0.2643 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1961) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1961) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.043, 0.25, -0.831, 0.609, -0.896, 0.886, 0.653, 0.065 };
   int lda = 2;
   double B[] = { 0.548, 0.076, 0.429, 0.873, -0.559, -0.329, -0.326, -0.174, 0.633, 0.489, 0.317, -0.896 };
   int ldb = 2;
   double B_expected[] = { 0.239257797324, 0.64684765886, 0.889006221152, 0.139062311692, 0.0322336011438, -0.807944179397, -0.977615509726, -1.02501063893, -0.164440783851, 0.983483814822, 1.28991055447, 1.90436729944 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1962) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1962) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.268, -0.062, -0.017, 0.326, 0.561, -0.203, -0.665, 0.338 };
   int lda = 2;
   double B[] = { -0.46, 0.954, 0.823, 0.945, -0.825, 0.882, -0.214, -0.095, -0.935, -0.245, 0.902, 0.904 };
   int ldb = 2;
   double B_expected[] = { 0.0426, -0.3322, -0.297862, -0.006188, 0.1593, -0.3471, 0.054794, 0.234161, 0.305, -0.02, -0.528045, -0.107865 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1963) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1963) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.075, 0.178, -0.321, -0.056, -0.124, -0.483, 0.685, -0.052 };
   int lda = 2;
   double B[] = { -0.47, -0.363, 0.766, -0.961, -0.391, -0.691, 0.42, -0.339, 0.45, -0.975, 0.991, -0.198 };
   int ldb = 2;
   double B_expected[] = { 0.874038948948, -0.779868445448, -0.234271045009, 0.514916650598, 0.810533012472, -1.05664738101, -0.149515922946, 0.198430908039, 2.17245126703, 0.115946317124, -0.420252834642, 0.199484456348 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1964) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1964) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.206, -0.461, -0.681, 0.358, 0.21, -0.318, 0.082, -0.097 };
   int lda = 2;
   double B[] = { 0.576, -0.249, 0.718, 0.424, 0.728, -0.464, 0.774, 0.541, -0.112, 0.803, 0.275, -0.638 };
   int ldb = 2;
   double B_expected[] = { -0.343295, 0.186865, -0.2578, -0.0554, -0.3973645, 0.2566785, -0.2863, -0.0849, 0.0189315, -0.0963345, -0.0187, 0.2189 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1965) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1965) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.117, 0.983, -0.929, -0.69, -0.144, 0.28, 0.658, 0.304, -0.657, 0.543, -0.051, -0.98, -0.846, -0.484, 0.052, 0.691, 0.613, -0.178 };
   int lda = 3;
   double B[] = { -0.688, 0.453, -0.63, 0.067, 0.193, 0.359, -0.792, 0.307, -0.501, -0.616, -0.595, 0.817 };
   int ldb = 3;
   double B_expected[] = { -0.566587593051, 0.340892661842, -0.458137993587, -0.0857620879204, -0.102500656517, -0.173972458173, -1.32599192297, -0.284341349955, -0.284178293736, -0.823318590512, 0.278700120014, -0.415972885216 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1966) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1966) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.773, -0.614, 0.782, -0.728, -0.727, -0.715, 0.858, -0.065, 0.922, 0.178, 0.588, 0.215, -0.92, -0.443, -0.583, -0.244, 0.996, -0.539 };
   int lda = 3;
   double B[] = { 0.159, 0.669, -0.692, 0.808, -0.146, 0.489, -0.385, -0.646, 0.704, -0.968, 0.551, -0.281 };
   int ldb = 3;
   double B_expected[] = { 0.0796383322, -0.0678193334, 0.0951193, -0.2156591, -0.0051, -0.1613, -0.2408434996, -0.0853028168, -0.0037554, 0.3083308, -0.1372, 0.1394 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1967) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1967) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.864, -0.382, -0.715, 0.227, -0.973, -0.709, -0.247, -0.601, 0.467, -0.133, 0.988, 0.937, -0.272, -0.334, 0.719, 0.992, 0.203, -0.646 };
   int lda = 3;
   double B[] = { 0.285, -0.409, -0.347, -0.925, -0.616, 0.422, 0.631, -0.954, -0.053, -0.255, -0.749, -0.979 };
   int ldb = 3;
   double B_expected[] = { -0.0215414266825, -0.165475896999, 0.469240391843, 0.538308411392, 1.71185240759, 0.063655952267, -0.0586080545035, -0.378370049976, 0.536158413721, 0.02961076215, 0.67769157898, -0.0939027988826 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1968) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1968) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.373, -0.335, -0.106, 0.542, -0.504, 0.574, -0.318, 0.043, -0.801, -0.331, 0.699, 0.776, -0.56, 0.131, 0.742, -0.692, -0.614, -0.874 };
   int lda = 3;
   double B[] = { -0.823, 0.929, -0.55, 0.172, -0.44, 0.067, 0.99, -0.013, 0.513, -0.438, -0.591, -0.302 };
   int ldb = 3;
   double B_expected[] = { 0.154, -0.361, 0.181249, -0.22802, 0.187552082, 0.008181148, -0.2957, 0.1029, -0.1997079, 0.2281373, 0.0457001502, -0.1796150434 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1969) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1969) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.912, 0.523, 0.314, -0.205, -0.895, 0.033, 0.157, -0.936, -0.582, 0.104, -0.868, 0.851, -0.131, 0.836, 0.993, 0.319, -0.684, -0.035 };
   int lda = 3;
   double B[] = { 0.07, -0.556, 0.018, -0.245, -0.405, 0.77, 0.888, 0.01, -0.81, -0.42, 0.66, -0.387 };
   int ldb = 2;
   double B_expected[] = { -0.132542904863, 0.151203976135, 0.45996395874, -0.700981460432, -0.771115355304, 0.0234040392321, 1.04091400336, -0.314874142966, -0.418936175202, -0.0443526810935, 0.218699329114, -0.27741882532 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1970) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1970) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.904, 0.983, -0.777, 0.503, 0.061, -0.442, 0.797, 0.415, -0.49, -0.466, 0.386, -0.147, -0.793, -0.381, -0.481, 0.33, 0.69, 0.35 };
   int lda = 3;
   double B[] = { 0.152, 0.832, 0.687, -0.287, -0.571, -0.187, -0.456, 0.631, 0.976, 0.833, -0.527, -0.188 };
   int ldb = 2;
   double B_expected[] = { -0.3155234788, -0.5211211034, -0.2870272698, 0.3910522396, -0.0411631, 0.0498567, 0.1600099, -0.2914973, -0.3761, -0.1523, 0.1769, 0.0037 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1971) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1971) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.527, 0.434, 0.025, 0.505, 0.724, 0.961, -0.071, 0.675, -0.334, 0.259, 0.167, 0.898, 0.116, 0.723, 0.086, 0.042, -0.483, -0.862 };
   int lda = 3;
   double B[] = { -0.874, 0.252, 0.924, 0.251, 0.559, -0.619, -0.131, -0.286, 0.09, -0.111, 0.062, -0.973 };
   int ldb = 2;
   double B_expected[] = { 0.116195543731, -0.404988360492, -0.325886265381, 0.300824742268, 0.86553022636, 0.0931927221532, -0.0931167995431, -0.760087414797, 0.774460770553, -0.204189465459, -0.501996021978, -0.354684266966 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1972) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1972) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.383, -0.184, 0.14, 0.131, -0.494, -0.025, -0.396, -0.183, 0.519, 0.806, -0.737, 0.764, -0.03, 0.622, -0.826, 0.605, 0.638, 0.935 };
   int lda = 3;
   double B[] = { 0.975, -0.816, -0.996, -0.038, -0.316, -0.31, -0.003, -0.974, 0.364, -0.217, 0.909, -0.656 };
   int ldb = 2;
   double B_expected[] = { -0.2109, 0.3423, 0.3026, -0.0882, 0.2001673, 0.0411059, 0.0443818, 0.2646074, -0.0213138923, 0.1426909311, 0.1794588402, 0.4128021586 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1973) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1973) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.433, -0.405, -0.008, 0.13, 0.377, -0.664, 0.421, -0.779 };
   int lda = 2;
   double B[] = { 0.022, -0.326, -0.905, 0.323, -0.722, 0.282, -0.877, -0.793, -0.906, -0.999, -0.607, -0.979 };
   int ldb = 3;
   double B_expected[] = { 0.0831887207906, -0.153137570623, -0.510564586332, -0.0447544052299, -0.412732352054, -0.0239182507667, 0.35364638809, -0.274824473121, 0.341954849059, -0.294570686181, 0.328230337479, -0.181800438645 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1974) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1974) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.007, 0.289, 0.434, 0.931, 0.776, -0.861, 0.83, -0.753 };
   int lda = 2;
   double B[] = { 0.775, -0.299, -0.45, 0.923, 0.251, 0.934, 0.388, -0.958, -0.732, 0.263, -0.5, 0.097 };
   int ldb = 3;
   double B_expected[] = { -0.2026, 0.1672, 0.0427, -0.3219, -0.1687, -0.2551, -0.0883348, 0.0650146, 0.4744571, 0.0273583, 0.4510139, -0.1254463 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1975) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1975) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.768, -0.738, 0.236, 0.721, 0.691, -0.963, -0.36, -0.376 };
   int lda = 2;
   double B[] = { -0.822, 0.174, 0.799, 0.8, -0.985, -0.169, 0.652, -0.529, -0.51, -0.506, -0.542, -0.786 };
   int ldb = 3;
   double B_expected[] = { -0.212429545832, 0.508667487335, 0.591670151369, 0.238559438419, 0.40264717438, -0.154881488703, 0.500259801606, -0.0994508738781, -0.130621162022, -0.416426547, -0.0684577231932, -0.575944733113 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1976) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1976) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.294, 0.843, 0.52, 0.53, 0.392, 0.293, 0.209, 0.497 };
   int lda = 2;
   double B[] = { 0.765, -0.547, 0.451, -0.581, 0.166, 0.834, -0.541, 0.278, -0.832, 0.66, -0.718, -0.664 };
   int ldb = 3;
   double B_expected[] = { -0.1872365, 0.3339085, -0.0667796, 0.3834252, -0.2809938, -0.2009734, 0.1345, -0.1375, 0.1836, -0.2812, 0.2818, 0.1274 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1977) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1977) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.642, 0.513, 0.428, 0.273, -0.612, 0.531, -0.664, 0.801 };
   int lda = 2;
   double B[] = { 0.429, -0.049, -0.661, 0.36, -0.247, 0.523, -0.227, 0.459, -0.902, 0.328, 0.37, -0.225 };
   int ldb = 2;
   double B_expected[] = { -0.161443909893, -0.0392846195877, 0.158306491417, 0.236544282705, 0.158671944063, -0.1560767799, 0.00300493937503, 0.254905467713, 0.369328020399, 0.00134777953987, -0.306971508873, -0.0836654236493 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1978) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1978) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.229, -0.461, -0.279, 0.674, -0.797, -0.286, 0.397, 0.329 };
   int lda = 2;
   double B[] = { 0.402, 0.728, 0.824, -0.691, -0.362, 0.437, 0.192, 0.788, -0.259, 0.599, 0.79, 0.076 };
   int ldb = 2;
   double B_expected[] = { -0.1934, -0.1782, -0.383205, 0.202987, 0.0649, -0.1673, -0.1325225, -0.3690995, 0.0178, -0.2056, -0.289215, -0.112754 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1979) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1979) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.942, -0.832, 0.595, -0.092, 0.01, 0.001, 0.944, 0.256 };
   int lda = 2;
   double B[] = { 0.73, 0.488, -0.363, -0.01, -0.112, 0.169, -0.268, -0.13, -0.657, 0.573, 0.91, 0.632 };
   int ldb = 2;
   double B_expected[] = { 0.158268746344, 0.226988691038, 0.117355164571, -0.00345029435376, -0.0289643723553, 0.0722018494696, 0.0888981803586, 0.0370317099277, -0.233113998714, -0.101765761072, -0.305361921327, -0.187259165106 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1980) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1980) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.508, 0.053, -0.516, 0.785, -0.451, -0.53, 0.551, 0.235 };
   int lda = 2;
   double B[] = { -0.09, 0.46, 0.948, 0.918, -0.337, 0.012, -0.786, -0.676, 0.906, -0.38, -0.566, 0.645 };
   int ldb = 2;
   double B_expected[] = { -0.0713482, -0.5355066, -0.3762, -0.1806, 0.1589574, 0.2649562, 0.3034, 0.1242, 0.0168633, 0.1582089, 0.1053, -0.2501 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1981) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1981) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.499, -0.268, 0.234, 0.032, -0.158, 0.684, -0.878, 0.613, 0.968, 0.812, 0.013, 0.34, -0.485, -0.565, 0.316, 0.286, -0.459, 0.637 };
   int lda = 3;
   double B[] = { -0.964, 0.804, 0.197, 0.141, 0.942, 0.474, 0.741, -0.441, -0.738, -0.703, -0.27, 0.98 };
   int ldb = 3;
   double B_expected[] = { 0.561582612433, -0.70128258354, -0.0253749021391, 0.0631927226609, 0.295313488523, -0.305260767297, -0.0937671252683, 0.884164549696, 0.000683977216651, 0.260184505619, 0.344358828778, 0.221445372699 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1982) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1982) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.453, 0.917, 0.131, 0.361, 0.087, 0.441, -0.439, 0.439, 0.777, 0.131, 0.535, 0.646, 0.508, 0.746, -0.347, -0.911, -0.874, -0.525 };
   int lda = 3;
   double B[] = { -0.739, -0.776, -0.049, 0.548, -0.39, -0.856, -0.757, 0.307, -0.533, -0.342, 0.431, 0.618 };
   int ldb = 3;
   double B_expected[] = { 0.2794424312, 0.1451980676, -0.2891898, -0.1549434, 0.2026, 0.2178, 0.2242026328, -0.0997909546, 0.3882643, 0.0019799, -0.1911, -0.1423 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1983) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1983) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.835, -0.775, -0.384, -0.128, -0.41, -0.511, -0.282, -0.341, -0.856, -0.662, 0.721, -0.939, 0.175, -0.899, 0.832, -0.519, 0.652, -0.318 };
   int lda = 3;
   double B[] = { -0.654, 0.105, -0.39, 0.645, 0.867, 0.045, -0.842, -0.896, -0.249, 0.419, 0.575, 0.561 };
   int ldb = 3;
   double B_expected[] = { -0.177337134492, -0.0485464421929, -0.0947130836909, 0.143712701441, -0.0502556531648, 0.286334558029, -0.109929498786, -0.323108217437, -0.0362323282558, 0.21056630482, -0.514117706819, 0.0792536824901 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1984) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1984) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.136, 0.272, 0.676, 0.673, -0.659, 0.668, 0.991, -0.569, -0.489, 0.581, -0.232, -0.249, -0.396, -0.832, 0.763, -0.092, 0.117, 0.108 };
   int lda = 3;
   double B[] = { 0.721, -0.141, -0.604, 0.318, 0.387, 0.73, -0.549, 0.302, 0.101, 0.721, -0.064, 0.673 };
   int ldb = 3;
   double B_expected[] = { -0.2022, 0.1144, 0.4148738, -0.1541186, -0.5047180206, 0.1126569022, 0.1345, -0.1455, -0.318479, -0.13854, 0.114359797, -0.242815912 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1985) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1985) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.578, 0.601, -0.43, -0.187, -0.934, 0.635, 0.157, -0.561, -0.964, 0.025, 0.435, -0.674, -0.575, 0.275, 0.609, 0.228, -0.202, -0.267 };
   int lda = 3;
   double B[] = { 0.505, -0.347, 0.213, -0.392, -0.465, -0.918, -0.737, -0.974, -0.051, 0.97, 0.066, 0.604 };
   int ldb = 2;
   double B_expected[] = { -0.206165616299, -0.811510964363, -0.328765954464, -0.593889594613, -0.410790365608, 0.365230809488, -0.377900693873, 0.166778025696, -0.558066070138, 0.728199798382, -0.271362172482, 0.505674752215 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1986) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1986) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.311, -0.737, -0.738, -0.214, -0.387, -0.043, -0.168, 0.563, 0.165, 0.007, -0.121, 0.408, -0.75, -0.641, -0.997, -0.347, 0.523, -0.922 };
   int lda = 3;
   double B[] = { 0.46, 0.376, -0.623, -0.092, 0.233, 0.981, -0.435, -0.493, 0.405, 0.855, -0.391, 0.572 };
   int ldb = 2;
   double B_expected[] = { -0.311417159, -0.418726217, 0.2053384662, -0.1587052684, -0.449331, -0.414523, 0.1666068, -0.1265226, -0.207, -0.216, 0.0601, -0.2107 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1987) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1987) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.288, -0.421, 0.451, 0.234, 0.67, -0.483, 0.273, 0.131, 0.005, 0.091, -0.706, -0.191, 0.285, -0.434, 0.648, -0.556, -0.886, 0.798 };
   int lda = 3;
   double B[] = { 0.359, -0.682, -0.618, 0.479, 0.463, 0.468, -0.43, 0.058, -0.361, -0.058, -0.028, -0.729 };
   int ldb = 2;
   double B_expected[] = { 0.432870841901, -0.202296442916, -0.484714722217, 0.00498299287046, -1.27917612947, -3.59551100448, 2.13407463306, 3.62604336509, 2.50059207751, 0.44116664838, -3.08374361183, -0.156015309482 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1988) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1988) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.288, 0.961, -0.571, -0.341, -0.443, 0.116, -0.928, 0.157, 0.035, 0.822, 0.733, -0.15, 0.851, -0.634, -0.769, -0.709, 0.346, -0.943 };
   int lda = 3;
   double B[] = { -0.708, 0.945, -0.144, 0.505, 0.827, -0.467, 0.883, 0.194, -0.607, -0.332, 0.716, -0.117 };
   int ldb = 2;
   double B_expected[] = { 0.1179, -0.3543, -0.0073, -0.1659, -0.2548954, -0.0197092, -0.3450402, -0.0621396, 0.4925104482, -0.0516973464, 0.0565040266, 0.1296638568 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1989) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1989) imag");
     };
   };
  };


}
