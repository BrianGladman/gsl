#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

#include "tests.h"

void
test_syr2k (void) {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 2;
   float alpha = -0.3;
   float beta = 0.1;
   float A[] = { -0.635, 0.805 };
   int lda = 2;
   float B[] = { 0.773, 0.375 };
   int ldb = 2;
   float C[] = { 0.616 };
   int ldc = 1;
   float C_expected[] = { 0.174988 };
   cblas_ssyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssyr2k(case 1622)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 2;
   float alpha = -0.3;
   float beta = 0.1;
   float A[] = { -0.396, -0.131 };
   int lda = 2;
   float B[] = { -0.603, -0.288 };
   int ldb = 2;
   float C[] = { -0.434 };
   int ldc = 1;
   float C_expected[] = { -0.20931 };
   cblas_ssyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssyr2k(case 1623)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 2;
   float alpha = -0.3;
   float beta = 0.1;
   float A[] = { -0.085, -0.444 };
   int lda = 1;
   float B[] = { 0.936, 0.752 };
   int ldb = 1;
   float C[] = { -0.64 };
   int ldc = 1;
   float C_expected[] = { 0.184069 };
   cblas_ssyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssyr2k(case 1624)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 2;
   float alpha = -0.3;
   float beta = 0.1;
   float A[] = { 0.655, 0.45 };
   int lda = 1;
   float B[] = { 0.16, -0.747 };
   int ldb = 1;
   float C[] = { 0.576 };
   int ldc = 1;
   float C_expected[] = { 0.19641 };
   cblas_ssyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssyr2k(case 1625)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 112;
   int N = 1;
   int K = 2;
   float alpha = -0.3;
   float beta = 0;
   float A[] = { 0.259, -0.334 };
   int lda = 1;
   float B[] = { -0.911, -0.426 };
   int ldb = 1;
   float C[] = { 0.432 };
   int ldc = 1;
   float C_expected[] = { 0.056199 };
   cblas_ssyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssyr2k(case 1626)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 112;
   int N = 1;
   int K = 2;
   float alpha = -0.3;
   float beta = 0;
   float A[] = { -0.765, 0.7 };
   int lda = 1;
   float B[] = { 0.487, 0.768 };
   int ldb = 1;
   float C[] = { 0.836 };
   int ldc = 1;
   float C_expected[] = { -0.099027 };
   cblas_ssyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssyr2k(case 1627)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 112;
   int N = 1;
   int K = 2;
   float alpha = -0.3;
   float beta = 0;
   float A[] = { -0.584, 0.056 };
   int lda = 2;
   float B[] = { 0.928, -0.101 };
   int ldb = 2;
   float C[] = { -0.529 };
   int ldc = 1;
   float C_expected[] = { 0.328565 };
   cblas_ssyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssyr2k(case 1628)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 112;
   int N = 1;
   int K = 2;
   float alpha = -0.3;
   float beta = 0;
   float A[] = { 0.25, 0.8 };
   int lda = 2;
   float B[] = { 0.489, -0.642 };
   int ldb = 2;
   float C[] = { 0.322 };
   int ldc = 1;
   float C_expected[] = { 0.23481 };
   cblas_ssyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssyr2k(case 1629)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 2;
   double alpha = 0;
   double beta = 0;
   double A[] = { 0.591, 0.21 };
   int lda = 2;
   double B[] = { -0.718, -0.579 };
   int ldb = 2;
   double C[] = { -0.856 };
   int ldc = 1;
   double C_expected[] = { -0 };
   cblas_dsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsyr2k(case 1630)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 2;
   double alpha = 0;
   double beta = 0;
   double A[] = { -0.971, -0.824 };
   int lda = 2;
   double B[] = { -0.227, 0.457 };
   int ldb = 2;
   double C[] = { 0.521 };
   int ldc = 1;
   double C_expected[] = { 0 };
   cblas_dsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsyr2k(case 1631)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 2;
   double alpha = 0;
   double beta = 0;
   double A[] = { -0.274, 0.583 };
   int lda = 1;
   double B[] = { 0.668, -0.83 };
   int ldb = 1;
   double C[] = { 0.907 };
   int ldc = 1;
   double C_expected[] = { 0 };
   cblas_dsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsyr2k(case 1632)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 2;
   double alpha = 0;
   double beta = 0;
   double A[] = { -0.512, -0.436 };
   int lda = 1;
   double B[] = { -0.443, -0.259 };
   int ldb = 1;
   double C[] = { -0.667 };
   int ldc = 1;
   double C_expected[] = { 0 };
   cblas_dsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsyr2k(case 1633)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 112;
   int N = 1;
   int K = 2;
   double alpha = -0.3;
   double beta = 0;
   double A[] = { 0.741, -0.341 };
   int lda = 1;
   double B[] = { 0.743, -0.315 };
   int ldb = 1;
   double C[] = { -0.776 };
   int ldc = 1;
   double C_expected[] = { -0.3947868 };
   cblas_dsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsyr2k(case 1634)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 112;
   int N = 1;
   int K = 2;
   double alpha = -0.3;
   double beta = 0;
   double A[] = { 0.03, 0.175 };
   int lda = 1;
   double B[] = { -0.832, 0.291 };
   int ldb = 1;
   double C[] = { 0.281 };
   int ldc = 1;
   double C_expected[] = { -0.015579 };
   cblas_dsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsyr2k(case 1635)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 112;
   int N = 1;
   int K = 2;
   double alpha = -0.3;
   double beta = 0;
   double A[] = { 0.476, 0.496 };
   int lda = 2;
   double B[] = { -0.626, -0.159 };
   int ldb = 2;
   double C[] = { -0.964 };
   int ldc = 1;
   double C_expected[] = { 0.226104 };
   cblas_dsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsyr2k(case 1636)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 112;
   int N = 1;
   int K = 2;
   double alpha = -0.3;
   double beta = 0;
   double A[] = { -0.489, 0.611 };
   int lda = 2;
   double B[] = { -0.285, -0.673 };
   int ldb = 2;
   double C[] = { -0.11 };
   int ldc = 1;
   double C_expected[] = { 0.1631028 };
   cblas_dsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsyr2k(case 1637)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 2;
   float alpha[2] = {1, 0};
   float beta[2] = {0, 0.1};
   float A[] = { 0.796, 0.872, -0.919, 0.748 };
   int lda = 2;
   float B[] = { -0.945, 0.915, -0.252, -0.276 };
   int ldb = 2;
   float C[] = { 0.07, -0.957 };
   int ldc = 1;
   float C_expected[] = { -2.12843, -0.054104 };
   cblas_csyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyr2k(case 1638) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyr2k(case 1638) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 2;
   float alpha[2] = {1, 0};
   float beta[2] = {0, 0.1};
   float A[] = { 0.984, 0.526, 0.284, 0.806 };
   int lda = 2;
   float B[] = { -0.509, -0.178, 0.188, -0.221 };
   int ldb = 2;
   float C[] = { -0.388, 0.795 };
   int ldc = 1;
   float C_expected[] = { -0.43092, -0.747044 };
   cblas_csyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyr2k(case 1639) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyr2k(case 1639) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 2;
   float alpha[2] = {1, 0};
   float beta[2] = {0, 0.1};
   float A[] = { -0.16, 0.628, -0.06, -0.645 };
   int lda = 1;
   float B[] = { 0.846, 0.545, 0.032, 0.493 };
   int ldb = 1;
   float C[] = { -0.041, -0.621 };
   int ldc = 1;
   float C_expected[] = { -0.26101, 0.783636 };
   cblas_csyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyr2k(case 1640) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyr2k(case 1640) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 2;
   float alpha[2] = {1, 0};
   float beta[2] = {0, 0.1};
   float A[] = { -0.478, -0.556, 0.519, 0.177 };
   int lda = 1;
   float B[] = { -0.946, 0.423, -0.859, 0.736 };
   int ldb = 1;
   float C[] = { -0.54, -0.035 };
   int ldc = 1;
   float C_expected[] = { 0.226066, 1.05345 };
   cblas_csyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyr2k(case 1641) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyr2k(case 1641) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 112;
   int N = 1;
   int K = 2;
   float alpha[2] = {0, 0};
   float beta[2] = {1, 0};
   float A[] = { -0.582, 0.09, -0.176, 0.784 };
   int lda = 1;
   float B[] = { 0.687, -0.859, 0.945, 0.756 };
   int ldb = 1;
   float C[] = { -0.663, -0.186 };
   int ldc = 1;
   float C_expected[] = { -0.663, -0.186 };
   cblas_csyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyr2k(case 1642) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyr2k(case 1642) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 112;
   int N = 1;
   int K = 2;
   float alpha[2] = {0, 0};
   float beta[2] = {1, 0};
   float A[] = { 0.231, -0.452, -0.112, -0.837 };
   int lda = 1;
   float B[] = { -0.258, 0.464, -0.224, 0.893 };
   int ldb = 1;
   float C[] = { -0.448, 0.046 };
   int ldc = 1;
   float C_expected[] = { -0.448, 0.046 };
   cblas_csyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyr2k(case 1643) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyr2k(case 1643) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 112;
   int N = 1;
   int K = 2;
   float alpha[2] = {0, 0};
   float beta[2] = {1, 0};
   float A[] = { 0.115, 0.178, -0.193, -0.491 };
   int lda = 2;
   float B[] = { 0.545, -0.665, 0.979, -0.4 };
   int ldb = 2;
   float C[] = { 0.522, 0.712 };
   int ldc = 1;
   float C_expected[] = { 0.522, 0.712 };
   cblas_csyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyr2k(case 1644) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyr2k(case 1644) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 112;
   int N = 1;
   int K = 2;
   float alpha[2] = {0, 0};
   float beta[2] = {1, 0};
   float A[] = { -0.725, -0.808, -0.244, 0.145 };
   int lda = 2;
   float B[] = { 0.447, -0.413, -0.226, -0.585 };
   int ldb = 2;
   float C[] = { -0.531, 0.227 };
   int ldc = 1;
   float C_expected[] = { -0.531, 0.227 };
   cblas_csyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyr2k(case 1645) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyr2k(case 1645) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 2;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {-0.3, 0.1};
   double A[] = { 0.337, -0.737, -0.993, 0.69 };
   int lda = 2;
   double B[] = { -0.39, -0.836, -0.32, 0.368 };
   int ldb = 2;
   double C[] = { 0.844, -0.763 };
   int ldc = 1;
   double C_expected[] = { 0.3494384, 0.5248712 };
   cblas_zsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyr2k(case 1646) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyr2k(case 1646) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 2;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {-0.3, 0.1};
   double A[] = { 0.386, -0.465, 0.719, -0.378 };
   int lda = 2;
   double B[] = { 0.099, -0.879, 0.864, 0.141 };
   int ldb = 2;
   double C[] = { -0.599, -0.47 };
   int ldc = 1;
   double C_expected[] = { 0.1664126, 0.5082238 };
   cblas_zsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyr2k(case 1647) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyr2k(case 1647) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 2;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {-0.3, 0.1};
   double A[] = { 0.914, 0.128, -0.262, -0.26 };
   int lda = 1;
   double B[] = { 0.431, 0.276, 0.75, 0.904 };
   int ldb = 1;
   double C[] = { 0.287, 0.537 };
   int ldc = 1;
   double C_expected[] = { -0.3532044, 0.0216788 };
   cblas_zsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyr2k(case 1648) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyr2k(case 1648) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 2;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {-0.3, 0.1};
   double A[] = { 0.618, 0.72, 0.392, -0.737 };
   int lda = 1;
   double B[] = { 0.783, 0.531, 0.375, 0.203 };
   int ldb = 1;
   double C[] = { 0.058, -0.116 };
   int ldc = 1;
   double C_expected[] = { -0.3837348, -0.2968344 };
   cblas_zsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyr2k(case 1649) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyr2k(case 1649) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 112;
   int N = 1;
   int K = 2;
   double alpha[2] = {0, 1};
   double beta[2] = {1, 0};
   double A[] = { -0.372, -0.735, -0.711, 0.051 };
   int lda = 1;
   double B[] = { 0.257, 0.097, 0.338, -0.484 };
   int ldb = 1;
   double C[] = { -0.142, -0.197 };
   int ldc = 1;
   double C_expected[] = { -0.414766, -0.676886 };
   cblas_zsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyr2k(case 1650) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyr2k(case 1650) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 112;
   int N = 1;
   int K = 2;
   double alpha[2] = {0, 1};
   double beta[2] = {1, 0};
   double A[] = { 0.1, -0.878, 0.28, -0.381 };
   int lda = 1;
   double B[] = { -0.208, 0.309, -0.276, 0.123 };
   int ldb = 1;
   double C[] = { 0.483, -0.541 };
   int ldc = 1;
   double C_expected[] = { -0.22324, -0.10083 };
   cblas_zsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyr2k(case 1651) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyr2k(case 1651) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 112;
   int N = 1;
   int K = 2;
   double alpha[2] = {0, 1};
   double beta[2] = {1, 0};
   double A[] = { -0.918, 0.515, -0.985, 0.067 };
   int lda = 2;
   double B[] = { -0.034, 0.408, 0.66, -0.945 };
   int ldb = 2;
   double C[] = { -0.063, -0.018 };
   int ldc = 1;
   double C_expected[] = { -1.228982, -1.549386 };
   cblas_zsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyr2k(case 1652) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyr2k(case 1652) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 112;
   int N = 1;
   int K = 2;
   double alpha[2] = {0, 1};
   double beta[2] = {1, 0};
   double A[] = { 0.443, -0.009, -0.245, -0.008 };
   int lda = 2;
   double B[] = { 0.495, -0.239, -0.973, -0.032 };
   int ldb = 2;
   double C[] = { -0.85, -0.799 };
   int ldc = 1;
   double C_expected[] = { -0.660584, 0.111526 };
   cblas_zsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyr2k(case 1653) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyr2k(case 1653) imag");
     };
   };
  };


}
