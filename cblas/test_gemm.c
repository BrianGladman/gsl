#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_gemm () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int transA = 111;
   int transB = 111;
   int M = 1;
   int N = 2;
   int K = 4;
   float alpha = 1;
   float beta = 0;
   float A[] = { 0.199, 0.237, 0.456, 0.377 };
   int lda = 4;
   float B[] = { 0.842, -0.734, 0.323, -0.957, -0.303, -0.873, -0.871, -0.819 };
   int ldb = 2;
   float C[] = { 0.498, -0.925 };
   int ldc = 2;
   float C_expected[] = { -0.222426, -1.07973 };
   cblas_sgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "sgemm(case 1466)");
     }
   };
  };


  {
   int order = 102;
   int transA = 111;
   int transB = 111;
   int M = 1;
   int N = 2;
   int K = 4;
   float alpha = 1;
   float beta = 0;
   float A[] = { -0.83, 0.922, -0.228, -0.003 };
   int lda = 1;
   float B[] = { 0.072, 0.345, 0.944, -0.39, -0.577, 0.656, -0.693, -0.453 };
   int ldb = 4;
   float C[] = { 0.583, 0.522 };
   int ldc = 1;
   float C_expected[] = { 0.044268, 1.24311 };
   cblas_sgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "sgemm(case 1467)");
     }
   };
  };


  {
   int order = 101;
   int transA = 111;
   int transB = 112;
   int M = 1;
   int N = 2;
   int K = 4;
   float alpha = 0.1;
   float beta = 0.1;
   float A[] = { -0.838, 0.622, -0.494, 0.304 };
   int lda = 4;
   float B[] = { 0.147, 0.134, 0.169, 0.734, -0.7, 0.541, -0.794, -0.256 };
   int ldb = 4;
   float C[] = { -0.632, -0.559 };
   int ldc = 2;
   float C_expected[] = { -0.0532188, 0.0678514 };
   cblas_sgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "sgemm(case 1468)");
     }
   };
  };


  {
   int order = 102;
   int transA = 111;
   int transB = 112;
   int M = 1;
   int N = 2;
   int K = 4;
   float alpha = 0.1;
   float beta = 0.1;
   float A[] = { -0.937, 0.635, 0.596, -0.51 };
   int lda = 1;
   float B[] = { -0.688, -0.265, 0.049, 0.133, -0.918, -0.147, 0.977, -0.21 };
   int ldb = 2;
   float C[] = { 0.844, 0.999 };
   int ldc = 1;
   float C_expected[] = { 0.0474373, 0.135125 };
   cblas_sgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "sgemm(case 1469)");
     }
   };
  };


  {
   int order = 101;
   int transA = 112;
   int transB = 111;
   int M = 1;
   int N = 2;
   int K = 4;
   float alpha = -0.3;
   float beta = 0.1;
   float A[] = { -0.165, 0.638, 0.346, -0.697 };
   int lda = 1;
   float B[] = { 0.499, -0.73, 0.262, 0.759, 0.664, 0.997, -0.702, -0.839 };
   int ldb = 2;
   float C[] = { 0.17, 0.425 };
   int ldc = 2;
   float C_expected[] = { -0.224158, -0.417831 };
   cblas_sgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "sgemm(case 1470)");
     }
   };
  };


  {
   int order = 102;
   int transA = 112;
   int transB = 111;
   int M = 1;
   int N = 2;
   int K = 4;
   float alpha = -0.3;
   float beta = 0.1;
   float A[] = { -0.603, -0.714, -0.893, 0.046 };
   int lda = 4;
   float B[] = { 0.859, -0.694, -0.868, -0.98, -0.103, 0.567, -0.277, -0.734 };
   int ldb = 4;
   float C[] = { 0.517, -0.622 };
   int ldc = 1;
   float C_expected[] = { -0.160575, -0.0234604 };
   cblas_sgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "sgemm(case 1471)");
     }
   };
  };


  {
   int order = 101;
   int transA = 112;
   int transB = 112;
   int M = 1;
   int N = 2;
   int K = 4;
   float alpha = 0.1;
   float beta = 1;
   float A[] = { -0.087, -0.047, -0.051, -0.615 };
   int lda = 1;
   float B[] = { -0.722, -0.077, 0.563, 0.501, 0.855, 0.605, 0.556, -0.627 };
   int ldb = 4;
   float C[] = { -0.181, -0.89 };
   int ldc = 2;
   float C_expected[] = { -0.208039, -0.864557 };
   cblas_sgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "sgemm(case 1472)");
     }
   };
  };


  {
   int order = 102;
   int transA = 112;
   int transB = 112;
   int M = 1;
   int N = 2;
   int K = 4;
   float alpha = 0.1;
   float beta = 1;
   float A[] = { -0.753, -0.074, -0.247, -0.19 };
   int lda = 4;
   float B[] = { 0.061, 0.743, 0.22, -0.682, 0.733, 0.417, 0.772, 0.665 };
   int ldb = 2;
   float C[] = { -0.253, 0.972 };
   int ldc = 1;
   float C_expected[] = { -0.291994, 0.898164 };
   cblas_sgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "sgemm(case 1473)");
     }
   };
  };


  {
   int order = 101;
   int transA = 111;
   int transB = 111;
   int M = 1;
   int N = 2;
   int K = 4;
   double alpha = 0;
   double beta = 0;
   double A[] = { 0.017, 0.191, 0.863, -0.97 };
   int lda = 4;
   double B[] = { -0.207, -0.916, -0.278, 0.403, 0.885, 0.409, -0.772, -0.27 };
   int ldb = 2;
   double C[] = { -0.274, -0.858 };
   int ldc = 2;
   double C_expected[] = { 0, 0 };
   cblas_dgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dgemm(case 1474)");
     }
   };
  };


  {
   int order = 102;
   int transA = 111;
   int transB = 111;
   int M = 1;
   int N = 2;
   int K = 4;
   double alpha = 0;
   double beta = 0;
   double A[] = { 0.571, 0.081, 0.109, 0.988 };
   int lda = 1;
   double B[] = { -0.048, -0.753, -0.8, -0.89, -0.535, -0.017, -0.018, -0.544 };
   int ldb = 4;
   double C[] = { -0.876, -0.792 };
   int ldc = 1;
   double C_expected[] = { 0, 0 };
   cblas_dgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dgemm(case 1475)");
     }
   };
  };


  {
   int order = 101;
   int transA = 111;
   int transB = 112;
   int M = 1;
   int N = 2;
   int K = 4;
   double alpha = -0.3;
   double beta = 1;
   double A[] = { 0.939, 0.705, 0.977, 0.4 };
   int lda = 4;
   double B[] = { -0.089, -0.822, 0.937, 0.159, 0.789, -0.413, -0.172, 0.88 };
   int ldb = 4;
   double C[] = { -0.619, 0.063 };
   int ldc = 2;
   double C_expected[] = { -0.7137904, -0.1270986 };
   cblas_dgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dgemm(case 1476)");
     }
   };
  };


  {
   int order = 102;
   int transA = 111;
   int transB = 112;
   int M = 1;
   int N = 2;
   int K = 4;
   double alpha = -0.3;
   double beta = 1;
   double A[] = { -0.795, 0.81, 0.388, 0.09 };
   int lda = 1;
   double B[] = { -0.847, 0.031, -0.938, 0.09, -0.286, -0.478, -0.981, 0.881 };
   int ldb = 2;
   double C[] = { -0.242, -0.02 };
   int ldc = 1;
   double C_expected[] = { -0.1562981, -0.0026243 };
   cblas_dgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dgemm(case 1477)");
     }
   };
  };


  {
   int order = 101;
   int transA = 112;
   int transB = 111;
   int M = 1;
   int N = 2;
   int K = 4;
   double alpha = -1;
   double beta = 0;
   double A[] = { -0.556, 0.532, 0.746, 0.673 };
   int lda = 1;
   double B[] = { -0.525, 0.967, 0.687, -0.024, 0.527, 0.485, 0.109, -0.46 };
   int ldb = 2;
   double C[] = { -0.495, 0.859 };
   int ldc = 2;
   double C_expected[] = { -1.123883, 0.49819 };
   cblas_dgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dgemm(case 1478)");
     }
   };
  };


  {
   int order = 102;
   int transA = 112;
   int transB = 111;
   int M = 1;
   int N = 2;
   int K = 4;
   double alpha = -1;
   double beta = 0;
   double A[] = { -0.358, 0.224, -0.941, 0.513 };
   int lda = 4;
   double B[] = { -0.201, -0.159, -0.586, -0.016, -0.324, 0.411, 0.115, -0.229 };
   int ldb = 4;
   double C[] = { 0.558, 0.596 };
   int ldc = 1;
   double C_expected[] = { -0.57956, 0.017636 };
   cblas_dgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dgemm(case 1479)");
     }
   };
  };


  {
   int order = 101;
   int transA = 112;
   int transB = 112;
   int M = 1;
   int N = 2;
   int K = 4;
   double alpha = -0.3;
   double beta = 1;
   double A[] = { -0.586, 0.809, 0.709, -0.524 };
   int lda = 1;
   double B[] = { 0.768, 0.7, 0.619, -0.478, -0.129, -0.778, -0.432, 0.454 };
   int ldb = 4;
   double C[] = { 0.042, 0.252 };
   int ldc = 2;
   double C_expected[] = { -0.1996785, 0.5813976 };
   cblas_dgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dgemm(case 1480)");
     }
   };
  };


  {
   int order = 102;
   int transA = 112;
   int transB = 112;
   int M = 1;
   int N = 2;
   int K = 4;
   double alpha = -0.3;
   double beta = 1;
   double A[] = { -0.164, 0.522, 0.948, -0.624 };
   int lda = 4;
   double B[] = { -0.142, 0.778, 0.359, 0.622, -0.637, -0.757, -0.282, -0.805 };
   int ldb = 2;
   double C[] = { -0.09, 0.183 };
   int ldc = 1;
   double C_expected[] = { -0.0248334, 0.1884672 };
   cblas_dgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dgemm(case 1481)");
     }
   };
  };


  {
   int order = 101;
   int transA = 111;
   int transB = 111;
   int M = 1;
   int N = 2;
   int K = 4;
   float alpha[2] = {0, 1};
   float beta[2] = {0, 0};
   float A[] = { -0.082, -0.281, -0.096, 0.913, 0.974, -0.706, -0.773, 0.522 };
   int lda = 4;
   float B[] = { 0.745, -0.664, 0.352, -0.733, 0.304, -0.555, -0.493, -0.089, 0.188, 0.631, 0.235, 0.152, -0.299, -0.731, -0.686, -0.332 };
   int ldb = 2;
   float C[] = { -0.179, -0.284, -0.996, -0.414 };
   int ldc = 2;
   float C_expected[] = { -1.06679, 1.47116, 0.599689, 0.933532 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1482) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1482) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 111;
   int transB = 111;
   int M = 1;
   int N = 2;
   int K = 4;
   float alpha[2] = {0, 1};
   float beta[2] = {0, 0};
   float A[] = { 0.044, -0.33, 0.279, 0.712, -0.363, -0.788, -0.768, -0.551 };
   int lda = 1;
   float B[] = { 0.138, 0.927, -0.178, -0.864, 0.888, 0.844, -0.199, 0.706, -0.034, 0.483, 0.499, 0.664, 0.648, 0.324, 0.97, 0.609 };
   int ldb = 4;
   float C[] = { -0.129, 0.842, 0.214, -0.626 };
   int ldc = 1;
   float C_expected[] = { 1.81122, 1.76205, 1.0574, -0.564966 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1483) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1483) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 111;
   int transB = 112;
   int M = 1;
   int N = 2;
   int K = 4;
   float alpha[2] = {0, 0};
   float beta[2] = {-1, 0};
   float A[] = { 0.812, -0.471, 0.241, 0.795, 0.439, 0.131, -0.636, 0.531 };
   int lda = 4;
   float B[] = { 0.062, 0.807, 0.873, 0.372, 0.239, 0.804, 0.537, -0.954, -0.396, 0.838, 0.081, 0.15, 0.489, -0.438, 0.165, 0.429 };
   int ldb = 4;
   float C[] = { 0.868, 0.329, -0.509, 0.724 };
   int ldc = 2;
   float C_expected[] = { -0.868, -0.329, 0.509, -0.724 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1484) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1484) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 111;
   int transB = 112;
   int M = 1;
   int N = 2;
   int K = 4;
   float alpha[2] = {0, 0};
   float beta[2] = {-1, 0};
   float A[] = { 0.832, 0.198, 0.794, -0.522, -0.319, 0.578, 0.332, 0.746 };
   int lda = 1;
   float B[] = { -0.361, 0.187, -0.163, -0.781, 0.536, 0.888, -0.969, 0.899, 0.961, -0.583, 0.753, 0.29, -0.997, 0.729, -0.352, -0.2 };
   int ldb = 2;
   float C[] = { 0.864, 0.735, -0.074, -0.228 };
   int ldc = 1;
   float C_expected[] = { -0.864, -0.735, 0.074, 0.228 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1485) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1485) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 111;
   int transB = 113;
   int M = 1;
   int N = 2;
   int K = 4;
   float alpha[2] = {0, 1};
   float beta[2] = {0, 0.1};
   float A[] = { 0.149, 0.187, 0.263, -0.715, -0.882, -0.907, 0.87, -0.527 };
   int lda = 4;
   float B[] = { -0.915, -0.249, -0.986, -0.799, -0.136, 0.712, 0.964, 0.799, -0.569, 0.686, 0.603, 0.758, 0.161, -0.698, -0.263, -0.256 };
   int ldb = 4;
   float C[] = { 0.622, -0.824, -0.482, -0.161 };
   int ldc = 2;
   float C_expected[] = { -0.246901, 0.083044, 1.25556, 0.009106 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1486) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1486) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 111;
   int transB = 113;
   int M = 1;
   int N = 2;
   int K = 4;
   float alpha[2] = {0, 1};
   float beta[2] = {0, 0.1};
   float A[] = { 0.963, -0.943, -0.734, -0.253, 0.832, 0.545, -0.815, -0.434 };
   int lda = 1;
   float B[] = { 0.23, -0.211, 0.906, 0.232, -0.339, 0.597, -0.919, 0.793, 0.535, 0.526, 0.119, 0.053, 0.751, 0.044, 0.752, -0.469 };
   int ldb = 2;
   float C[] = { 0.483, -0.266, -0.224, -0.692 };
   int ldc = 1;
   float C_expected[] = { -0.047537, 0.667177, 1.02025, 0.823778 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1487) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1487) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 112;
   int transB = 111;
   int M = 1;
   int N = 2;
   int K = 4;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {-1, 0};
   float A[] = { -0.657, -0.497, -0.293, -0.168, -0.943, -0.181, 0.569, 0.91 };
   int lda = 1;
   float B[] = { -0.047, 0.796, -0.913, 0.998, 0.365, 0.467, -0.627, -0.523, 0.885, 0.234, -0.494, 0.071, -0.361, -0.154, -0.055, -0.32 };
   int ldb = 2;
   float C[] = { 0.956, 0.268, 0.152, 0.717 };
   int ldc = 2;
   float C_expected[] = { -0.668685, 0.134477, -0.715786, -0.478065 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1488) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1488) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 112;
   int transB = 111;
   int M = 1;
   int N = 2;
   int K = 4;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {-1, 0};
   float A[] = { 0.394, -0.482, 0.631, -0.833, 0.221, 0.672, 0.2, 0.967 };
   int lda = 4;
   float B[] = { 0.708, 0.695, 0.111, -0.912, 0.376, 0.606, -0.997, -0.741, 0.349, 0.543, 0.372, -0.563, 0.129, -0.295, -0.672, -0.95 };
   int ldb = 4;
   float C[] = { 0.436, 0.752, 0.074, 0.209 };
   int ldc = 1;
   float C_expected[] = { -0.325083, -0.301952, -0.283022, 0.339919 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1489) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1489) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 112;
   int transB = 112;
   int M = 1;
   int N = 2;
   int K = 4;
   float alpha[2] = {1, 0};
   float beta[2] = {-0.3, 0.1};
   float A[] = { 0.827, -0.862, 0.373, -0.265, -0.9, 0.892, -0.319, 0.151 };
   int lda = 1;
   float B[] = { 0.603, 0.816, -0.511, 0.831, -0.36, -0.954, -0.978, 0.485, 0.675, 0.186, 0.463, 0.144, 0.851, -0.458, 0.766, -0.213 };
   int ldb = 4;
   float C[] = { -0.335, 0.333, -0.4, 0.422 };
   int ldc = 2;
   float C_expected[] = { 2.7126, 0.702111, 0.437661, 0.691294 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1490) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1490) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 112;
   int transB = 112;
   int M = 1;
   int N = 2;
   int K = 4;
   float alpha[2] = {1, 0};
   float beta[2] = {-0.3, 0.1};
   float A[] = { 0.966, 0.476, -0.013, -0.655, 0.773, -0.543, -0.231, -0.353 };
   int lda = 4;
   float B[] = { -0.684, 0.144, 0.018, -0.77, -0.688, 0.909, -0.094, -0.938, -0.757, 0.574, -0.479, 0.473, 0, 0.064, -0.168, 0.858 };
   int ldb = 2;
   float C[] = { -0.912, 0.54, 0.756, 0.024 };
   int ldc = 1;
   float C_expected[] = { -0.156236, 0.839112, -0.230206, -0.106256 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1491) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1491) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 112;
   int transB = 113;
   int M = 1;
   int N = 2;
   int K = 4;
   float alpha[2] = {0, 0};
   float beta[2] = {1, 0};
   float A[] = { 0.66, -0.113, -0.663, -0.856, 0.614, -0.344, -0.964, -0.532 };
   int lda = 1;
   float B[] = { -0.606, -0.965, -0.279, -0.312, 0.63, 0.967, 0.041, -0.557, 0.663, 0.619, -0.134, 0.261, -0.388, 0.525, 0.222, 0.538 };
   int ldb = 4;
   float C[] = { 0.114, -0.376, -0.851, -0.682 };
   int ldc = 2;
   float C_expected[] = { 0.114, -0.376, -0.851, -0.682 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1492) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1492) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 112;
   int transB = 113;
   int M = 1;
   int N = 2;
   int K = 4;
   float alpha[2] = {0, 0};
   float beta[2] = {1, 0};
   float A[] = { 0.212, -0.752, 0.679, 0.49, -0.029, -0.488, 0.567, 0.374 };
   int lda = 4;
   float B[] = { -0.914, 0.734, -0.845, 0.059, -0.297, 0.152, -0.417, -0.669, 0.831, -0.544, 0.022, 0.102, -0.379, -0.357, -0.394, -0.588 };
   int ldb = 2;
   float C[] = { -0.584, 0.373, 0.235, 0.521 };
   int ldc = 1;
   float C_expected[] = { -0.584, 0.373, 0.235, 0.521 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1493) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1493) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 113;
   int transB = 111;
   int M = 1;
   int N = 2;
   int K = 4;
   float alpha[2] = {0, 0.1};
   float beta[2] = {-1, 0};
   float A[] = { 0.135, 0.128, 0.909, -0.963, 0.299, -0.944, 0.944, 0.942 };
   int lda = 1;
   float B[] = { 0.924, -0.317, -0.992, -0.854, -0.435, 0.102, 0.126, 0.862, 0.952, 0.68, 0.545, 0.168, 0.752, 0.549, 0.687, -0.76 };
   int ldb = 2;
   float C[] = { -0.369, -0.33, 0.849, -0.632 };
   int ldc = 2;
   float C_expected[] = { 0.326537, 0.37603, -0.86067, 0.529817 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1494) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1494) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 113;
   int transB = 111;
   int M = 1;
   int N = 2;
   int K = 4;
   float alpha[2] = {0, 0.1};
   float beta[2] = {-1, 0};
   float A[] = { 0.061, -0.271, -0.043, -0.023, 0.694, 0.333, 0.733, -0.967 };
   int lda = 4;
   float B[] = { 0.088, -0.607, 0.589, 0.375, -0.897, -0.954, -0.216, -0.195, -0.865, -0.511, -0.219, 0.535, 0.976, 0.582, 0.464, -0.041 };
   int ldb = 4;
   float C[] = { 0.533, -0.63, 0.405, 0.667 };
   int ldc = 1;
   float C_expected[] = { -0.459906, 0.552595, -0.425391, -0.533626 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1495) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1495) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 113;
   int transB = 112;
   int M = 1;
   int N = 2;
   int K = 4;
   float alpha[2] = {0, 0};
   float beta[2] = {1, 0};
   float A[] = { -0.676, -0.116, 0.707, -0.256, -0.893, -0.966, 0.159, -0.246 };
   int lda = 1;
   float B[] = { 0.059, 0.281, -0.93, -0.263, 0.583, -0.11, 0.639, -0.96, -0.878, 0.984, 0.058, 0.977, -0.567, 0.561, -0.048, -0.798 };
   int ldb = 4;
   float C[] = { 0.362, -0.808, 0.428, -0.112 };
   int ldc = 2;
   float C_expected[] = { 0.362, -0.808, 0.428, -0.112 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1496) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1496) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 113;
   int transB = 112;
   int M = 1;
   int N = 2;
   int K = 4;
   float alpha[2] = {0, 0};
   float beta[2] = {1, 0};
   float A[] = { -0.915, 0.439, 0.171, -0.019, 0.843, 0.944, -0.581, 0.856 };
   int lda = 4;
   float B[] = { -0.284, 0.207, -0.27, 0.832, 0.894, -0.626, -0.305, -0.006, 0.562, -0.744, -0.533, 0.126, -0.375, -0.333, 0.275, 0.748 };
   int ldb = 2;
   float C[] = { -0.763, -0.829, 0.708, -0.613 };
   int ldc = 1;
   float C_expected[] = { -0.763, -0.829, 0.708, -0.613 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1497) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1497) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 113;
   int transB = 113;
   int M = 1;
   int N = 2;
   int K = 4;
   float alpha[2] = {0, 0.1};
   float beta[2] = {0, 1};
   float A[] = { 0.496, -0.9, 0.825, -0.678, 0.41, -0.585, -0.264, 0.308 };
   int lda = 1;
   float B[] = { 0.907, 0.972, -0.724, 0.745, -0.601, 0.589, 0.759, -0.521, -0.161, -0.321, 0.341, -0.981, -0.378, -0.671, -0.314, -0.878 };
   int ldb = 4;
   float C[] = { -0.293, 0.07, 0.087, -0.542 };
   int ldc = 2;
   float C_expected[] = { 0.10357, -0.163927, 0.444626, -0.0076744 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1498) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1498) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 113;
   int transB = 113;
   int M = 1;
   int N = 2;
   int K = 4;
   float alpha[2] = {0, 0.1};
   float beta[2] = {0, 1};
   float A[] = { -0.225, -0.629, -0.939, -0.836, -0.841, -0.794, 0.836, -0.65 };
   int lda = 4;
   float B[] = { 0.869, -0.453, 0.8, -0.947, 0.545, 0.716, -0.507, -0.228, 0.722, 0.372, 0.77, 0.317, -0.153, -0.524, -0.465, -0.684 };
   int ldb = 2;
   float C[] = { -0.896, 0.91, -0.973, -0.269 };
   int ldc = 1;
   float C_expected[] = { -1.18974, -1.0134, 0.189027, -1.14494 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1499) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1499) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 111;
   int transB = 111;
   int M = 1;
   int N = 2;
   int K = 4;
   double alpha[2] = {1, 0};
   double beta[2] = {-1, 0};
   double A[] = { -0.33, 0.457, 0.428, -0.19, 0.86, -0.53, 0.058, -0.942 };
   int lda = 4;
   double B[] = { 0.434, 0.653, -0.124, 0.191, -0.112, -0.84, -0.72, 0.075, -0.503, -0.109, 0.3, -0.898, 0.489, 0.384, 0.993, -0.804 };
   int ldb = 2;
   double C[] = { -0.792, -0.155, -0.608, -0.243 };
   int ldc = 2;
   double C_expected[] = { 0.042563, -0.465908, -0.649991, -1.621116 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1500) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1500) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 111;
   int transB = 111;
   int M = 1;
   int N = 2;
   int K = 4;
   double alpha[2] = {1, 0};
   double beta[2] = {-1, 0};
   double A[] = { 0.726, -0.438, -0.23, -0.054, -0.019, 0.902, -0.883, -0.235 };
   int lda = 1;
   double B[] = { 0.159, -0.18, 0.386, -0.167, 0.971, -0.072, 0.87, -0.839, 0.474, 0.956, -0.235, 0.332, 0.826, -0.056, -0.941, 0.01 };
   int ldb = 4;
   double C[] = { -0.799, 0.973, -0.549, -0.177 };
   int ldc = 1;
   double C_expected[] = { -0.181084, 0.257841, 2.251901, 1.558195 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1501) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1501) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 111;
   int transB = 112;
   int M = 1;
   int N = 2;
   int K = 4;
   double alpha[2] = {0, 0.1};
   double beta[2] = {1, 0};
   double A[] = { 0.109, 0.892, -0.723, 0.793, 0.109, -0.419, -0.534, 0.448 };
   int lda = 4;
   double B[] = { -0.875, -0.31, -0.027, 0.067, 0.274, -0.126, -0.548, 0.497, 0.681, 0.388, 0.909, 0.889, 0.982, -0.074, -0.788, 0.233 };
   int ldb = 4;
   double C[] = { 0.503, 0.067, 0.239, 0.876 };
   int ldc = 2;
   double C_expected[] = { 0.6553584, 0.0864583, 0.2559136, 0.7518389 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1502) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1502) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 111;
   int transB = 112;
   int M = 1;
   int N = 2;
   int K = 4;
   double alpha[2] = {0, 0.1};
   double beta[2] = {1, 0};
   double A[] = { 0.334, 0.192, -0.992, -0.168, 0.154, -0.75, -0.797, -0.76 };
   int lda = 1;
   double B[] = { -0.82, 0.147, -0.237, 0.68, 0.317, 0.257, -0.406, -0.802, 0.058, 0.012, -0.832, 0.949, -0.263, -0.085, -0.064, 0.492 };
   int ldb = 2;
   double C[] = { 0.079, -0.602, -0.392, 0.316 };
   int ldc = 1;
   double C_expected[] = { 0.0980569, -0.6430449, -0.539207, 0.4226848 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1503) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1503) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 111;
   int transB = 113;
   int M = 1;
   int N = 2;
   int K = 4;
   double alpha[2] = {0, 0};
   double beta[2] = {-1, 0};
   double A[] = { -0.305, -0.698, -0.072, -0.383, 0.364, -0.656, 0.819, 0.194 };
   int lda = 4;
   double B[] = { 0.682, 0.498, -0.389, 0.923, -0.853, -0.558, -0.722, -0.085, -0.27, 0.026, -0.107, -0.036, 0.644, -0.327, -0.894, 0.34 };
   int ldb = 4;
   double C[] = { 0.981, -0.336, -0.377, -0.41 };
   int ldc = 2;
   double C_expected[] = { -0.981, 0.336, 0.377, 0.41 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1504) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1504) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 111;
   int transB = 113;
   int M = 1;
   int N = 2;
   int K = 4;
   double alpha[2] = {0, 0};
   double beta[2] = {-1, 0};
   double A[] = { -0.306, -0.709, -0.196, 0.285, 0.873, -0.802, 0.715, -0.179 };
   int lda = 1;
   double B[] = { 0.028, 0.109, 0.87, -0.446, 0.735, 0.731, 0.021, -0.186, 0.541, 0.97, -0.333, 0.002, -0.089, -0.01, 0.331, 0.851 };
   int ldb = 2;
   double C[] = { 0.902, -0.584, -0.695, -0.607 };
   int ldc = 1;
   double C_expected[] = { -0.902, 0.584, 0.695, 0.607 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1505) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1505) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 112;
   int transB = 111;
   int M = 1;
   int N = 2;
   int K = 4;
   double alpha[2] = {-1, 0};
   double beta[2] = {1, 0};
   double A[] = { 0.517, -0.136, 0.72, -0.237, 0.121, -0.66, 0.005, 0.759 };
   int lda = 1;
   double B[] = { -0.606, 0.049, 0.807, -0.236, -0.258, -0.412, 0.75, -0.659, 0.993, -0.029, -0.968, 0.707, -0.362, -0.005, 0.096, -0.241 };
   int ldb = 2;
   double C[] = { 0.63, 0.922, 0.025, -0.535 };
   int ldc = 2;
   double C_expected[] = { 1.117044, 1.983417, -1.276831, -0.447092 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1506) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1506) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 112;
   int transB = 111;
   int M = 1;
   int N = 2;
   int K = 4;
   double alpha[2] = {-1, 0};
   double beta[2] = {1, 0};
   double A[] = { 0.064, 0.371, -0.01, -0.262, 0.143, -0.081, 0.1, -0.062 };
   int lda = 4;
   double B[] = { -0.749, 0.289, -0.239, -0.226, 0.284, 0.668, 0.305, 0.075, -0.36, 0.166, -0.416, 0.234, -0.267, 0.525, 0.116, -0.561 };
   int ldb = 4;
   double C[] = { 0.671, 0.763, 0.444, -0.246 };
   int ldc = 1;
   double C_expected[] = { 0.753107, 0.896395, 0.481996, -0.263126 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1507) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1507) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 112;
   int transB = 112;
   int M = 1;
   int N = 2;
   int K = 4;
   double alpha[2] = {1, 0};
   double beta[2] = {-0.3, 0.1};
   double A[] = { -0.956, -0.751, 0.671, -0.633, 0.648, -0.042, 0.948, 0.826 };
   int lda = 1;
   double B[] = { 0.921, 0.506, -0.609, 0.817, -0.686, 0.991, 0.616, -0.482, -0.02, -0.34, 0.559, 0.976, 0.431, 0.385, -0.164, -0.778 };
   int ldb = 4;
   double C[] = { 0.074, -0.01, 0.165, 0.166 };
   int ldc = 2;
   double C_expected[] = { 0.166046, 0.491557, 1.473191, -0.033821 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1508) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1508) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 112;
   int transB = 112;
   int M = 1;
   int N = 2;
   int K = 4;
   double alpha[2] = {1, 0};
   double beta[2] = {-0.3, 0.1};
   double A[] = { -0.698, -0.062, 0.023, 0.704, 0.443, -0.46, 0.541, 0.296 };
   int lda = 4;
   double B[] = { 0.787, -0.199, 0.835, -0.276, -0.515, 0.467, -0.76, -0.483, 0.015, -0.394, -0.748, 0.02, 0.573, 0.3, -0.088, -0.238 };
   int ldb = 2;
   double C[] = { 0.935, -0.655, -0.797, 0.071 };
   int ldc = 1;
   double C_expected[] = { -1.070679, 0.178755, -0.344714, -0.308137 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1509) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1509) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 112;
   int transB = 113;
   int M = 1;
   int N = 2;
   int K = 4;
   double alpha[2] = {0, 0.1};
   double beta[2] = {-0.3, 0.1};
   double A[] = { -0.202, -0.219, 0.741, 0.527, 0.054, 0.16, -0.359, 0.338 };
   int lda = 1;
   double B[] = { -0.872, 0.995, 0.722, 0.618, -0.27, 0.939, -0.743, 0.547, -0.864, 0.376, -0.997, -0.63, 0.887, -0.454, 0.436, -0.039 };
   int ldb = 4;
   double C[] = { -0.684, 0.463, -0.386, -0.524 };
   int ldc = 2;
   double C_expected[] = { 0.1423153, -0.066679, 0.1175618, 0.0012949 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1510) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1510) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 112;
   int transB = 113;
   int M = 1;
   int N = 2;
   int K = 4;
   double alpha[2] = {0, 0.1};
   double beta[2] = {-0.3, 0.1};
   double A[] = { -0.855, -0.173, -0.679, 0.824, 0.469, 0.786, 0.757, -0.109 };
   int lda = 4;
   double B[] = { 0.483, -0.888, -0.757, 0.551, -0.81, 0.23, -0.078, 0.725, -0.592, 0.394, 0.884, 0.802, -0.813, -0.016, -0.853, 0.783 };
   int ldb = 2;
   double C[] = { 0.181, -0.368, -0.864, -0.784 };
   int ldc = 1;
   double C_expected[] = { 0.1728438, 0.1183508, 0.2526999, 0.3004174 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1511) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1511) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 113;
   int transB = 111;
   int M = 1;
   int N = 2;
   int K = 4;
   double alpha[2] = {-1, 0};
   double beta[2] = {-0.3, 0.1};
   double A[] = { 0.446, -0.65, -0.724, 0.014, 0.792, -0.695, -0.81, -0.358 };
   int lda = 1;
   double B[] = { -0.08, 0.216, 0.689, 0.699, 0.073, -0.346, 0.821, -0.668, -0.798, 0.869, 0.451, -0.061, -0.41, 0.316, 0.104, -0.514 };
   int ldb = 2;
   double C[] = { -0.476, 0.211, -0.912, -0.243 };
   int ldc = 2;
   double C_expected[] = { 1.372475, -0.135616, 0.549353, -1.968747 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1512) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1512) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 113;
   int transB = 111;
   int M = 1;
   int N = 2;
   int K = 4;
   double alpha[2] = {-1, 0};
   double beta[2] = {-0.3, 0.1};
   double A[] = { 0.669, 0.046, -0.094, 0.666, 0.23, 0.448, -0.795, -0.142 };
   int lda = 4;
   double B[] = { 0.037, -0.154, -0.739, 0.905, 0.793, -0.53, -0.34, 0.428, 0.072, -0.263, -0.603, -0.905, 0.681, -0.083, -0.511, -0.337 };
   int ldb = 4;
   double C[] = { 0.247, 0.575, -0.836, -0.883 };
   int ldc = 1;
   double C_expected[] = { -0.975939, 0.415528, 0.275533, 0.002716 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1513) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1513) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 113;
   int transB = 112;
   int M = 1;
   int N = 2;
   int K = 4;
   double alpha[2] = {0, 0};
   double beta[2] = {-1, 0};
   double A[] = { 0.369, 0.506, 0.217, -0.739, -0.395, 0.16, -0.329, -0.954 };
   int lda = 1;
   double B[] = { -0.622, -0.945, 0.416, -0.884, 0.797, -0.74, 0.519, -0.789, -0.348, 0.563, -0.398, -0.956, 0.227, 0.84, -0.079, 0.847 };
   int ldb = 4;
   double C[] = { 0.833, 0.761, 0.074, -0.448 };
   int ldc = 2;
   double C_expected[] = { -0.833, -0.761, -0.074, 0.448 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1514) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1514) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 113;
   int transB = 112;
   int M = 1;
   int N = 2;
   int K = 4;
   double alpha[2] = {0, 0};
   double beta[2] = {-1, 0};
   double A[] = { -0.141, 0.275, 0.717, 0.775, -0.701, -0.689, -0.883, -0.077 };
   int lda = 4;
   double B[] = { -0.526, -0.437, 0.133, -0.209, -0.83, 0.328, 0.916, -0.337, 0.762, -0.664, -0.566, 0.955, 0.168, 0.488, -0.172, -0.535 };
   int ldb = 2;
   double C[] = { -0.88, 0.945, 0.416, 0.99 };
   int ldc = 1;
   double C_expected[] = { 0.88, -0.945, -0.416, -0.99 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1515) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1515) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 113;
   int transB = 113;
   int M = 1;
   int N = 2;
   int K = 4;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {0, 0.1};
   double A[] = { -0.534, -0.013, -0.258, -0.31, -0.211, -0.883, -0.89, -0.499 };
   int lda = 1;
   double B[] = { -0.185, -0.798, -0.34, 0.716, 0.035, 0.968, -0.26, 0.784, -0.889, -0.344, -0.685, -0.647, -0.764, 0.03, 0.626, -0.989 };
   int ldb = 4;
   double C[] = { -0.793, -0.551, 0.182, 0.838 };
   int ldc = 2;
   double C_expected[] = { -0.5507177, -0.0286821, 0.2222276, 0.5197398 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1516) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1516) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 113;
   int transB = 113;
   int M = 1;
   int N = 2;
   int K = 4;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {0, 0.1};
   double A[] = { 0.575, -0.128, -0.702, 0.758, 0.383, -0.914, 0.157, 0.368 };
   int lda = 4;
   double B[] = { 0.572, -0.841, 0.223, -0.334, -0.823, -0.84, 0.671, -0.871, 0.241, 0.927, -0.344, 0.281, -0.034, -0.104, 0.587, -0.329 };
   int ldb = 2;
   double C[] = { -0.612, 0.167, 0.647, 0.447 };
   int ldc = 1;
   double C_expected[] = { -0.7876717, 0.0341179, -0.0800018, 0.5717566 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1517) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1517) imag");
     };
   };
  };


}
