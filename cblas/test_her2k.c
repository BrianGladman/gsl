#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_her2k () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 1;
   float alpha[2] = {-0.3, 0.1};
   float beta = 0.1;
   float A[] = { 0.17, 0.425 };
   int lda = 1;
   float B[] = { -0.165, 0.638 };
   int ldb = 1;
   float C[] = { 0.776, -0.536 };
   int ldc = 1;
   float C_expected[] = { -0.032543, 0 };
   cblas_cher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cher2k(case 1376) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cher2k(case 1376) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 1;
   float alpha[2] = {-0.3, 0.1};
   float beta = 0.1;
   float A[] = { 0.499, -0.73 };
   int lda = 1;
   float B[] = { 0.262, 0.759 };
   int ldb = 1;
   float C[] = { 0.346, -0.697 };
   int ldc = 1;
   float C_expected[] = { 0.4025994, 0 };
   cblas_cher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cher2k(case 1377) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cher2k(case 1377) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 1;
   float alpha[2] = {-0.3, 0.1};
   float beta = 0.1;
   float A[] = { -0.702, -0.839 };
   int lda = 1;
   float B[] = { 0.517, -0.622 };
   int ldb = 1;
   float C[] = { 0.664, 0.997 };
   int ldc = 1;
   float C_expected[] = { 0.145127, 0 };
   cblas_cher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cher2k(case 1378) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cher2k(case 1378) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 1;
   float alpha[2] = {-0.3, 0.1};
   float beta = 0.1;
   float A[] = { -0.893, 0.046 };
   int lda = 1;
   float B[] = { 0.859, -0.694 };
   int ldb = 1;
   float C[] = { -0.603, -0.714 };
   int ldc = 1;
   float C_expected[] = { 0.5351522, 0 };
   cblas_cher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cher2k(case 1379) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cher2k(case 1379) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 113;
   int N = 1;
   int K = 1;
   float alpha[2] = {0, 0};
   float beta = 0;
   float A[] = { -0.277, -0.734 };
   int lda = 1;
   float B[] = { -0.5, -0.101 };
   int ldb = 1;
   float C[] = { -0.103, 0.567 };
   int ldc = 1;
   float C_expected[] = { 0.000000000000e+00, 0.000000000000e+00 };
   cblas_cher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cher2k(case 1380) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cher2k(case 1380) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 113;
   int N = 1;
   int K = 1;
   float alpha[2] = {0, 0};
   float beta = 0;
   float A[] = { -0.087, -0.047 };
   int lda = 1;
   float B[] = { -0.051, -0.615 };
   int ldb = 1;
   float C[] = { -0.181, -0.89 };
   int ldc = 1;
   float C_expected[] = { 0.000000000000e+00, 0.000000000000e+00 };
   cblas_cher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cher2k(case 1381) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cher2k(case 1381) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 113;
   int N = 1;
   int K = 1;
   float alpha[2] = {0, 0};
   float beta = 0;
   float A[] = { 0.563, 0.501 };
   int lda = 1;
   float B[] = { 0.855, 0.605 };
   int ldb = 1;
   float C[] = { -0.722, -0.077 };
   int ldc = 1;
   float C_expected[] = { 0.000000000000e+00, 0.000000000000e+00 };
   cblas_cher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cher2k(case 1382) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cher2k(case 1382) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 113;
   int N = 1;
   int K = 1;
   float alpha[2] = {0, 0};
   float beta = 0;
   float A[] = { -0.253, 0.972 };
   int lda = 1;
   float B[] = { -0.753, -0.074 };
   int ldb = 1;
   float C[] = { 0.556, -0.627 };
   int ldc = 1;
   float C_expected[] = { 0.000000000000e+00, 0.000000000000e+00 };
   cblas_cher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cher2k(case 1383) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cher2k(case 1383) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 1;
   double alpha[2] = {1, 0};
   double beta = 1;
   double A[] = { 0.22, -0.682 };
   int lda = 1;
   double B[] = { 0.733, 0.417 };
   int ldb = 1;
   double C[] = { 0.061, 0.743 };
   int ldc = 1;
   double C_expected[] = { -0.185268, 0 };
   cblas_zher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zher2k(case 1384) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zher2k(case 1384) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 1;
   double alpha[2] = {1, 0};
   double beta = 1;
   double A[] = { -0.619, -0.76 };
   int lda = 1;
   double B[] = { -0.274, -0.858 };
   int ldb = 1;
   double C[] = { 0.772, 0.665 };
   int ldc = 1;
   double C_expected[] = { 2.415372, 0 };
   cblas_zher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zher2k(case 1385) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zher2k(case 1385) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 1;
   double alpha[2] = {1, 0};
   double beta = 1;
   double A[] = { 0.863, -0.97 };
   int lda = 1;
   double B[] = { -0.207, -0.916 };
   int ldb = 1;
   double C[] = { 0.017, 0.191 };
   int ldc = 1;
   double C_expected[] = { 1.436758, 0 };
   cblas_zher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zher2k(case 1386) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zher2k(case 1386) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 1;
   double alpha[2] = {1, 0};
   double beta = 1;
   double A[] = { 0.885, 0.409 };
   int lda = 1;
   double B[] = { -0.772, -0.27 };
   int ldb = 1;
   double C[] = { -0.278, 0.403 };
   int ldc = 1;
   double C_expected[] = { -1.8653, 0 };
   cblas_zher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zher2k(case 1387) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zher2k(case 1387) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 113;
   int N = 1;
   int K = 1;
   double alpha[2] = {0, 0};
   double beta = 0;
   double A[] = { 0.109, 0.988 };
   int lda = 1;
   double B[] = { -0.048, -0.753 };
   int ldb = 1;
   double C[] = { 0.571, 0.081 };
   int ldc = 1;
   double C_expected[] = { 0.000000000000e+00, 0.000000000000e+00 };
   cblas_zher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zher2k(case 1388) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zher2k(case 1388) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 113;
   int N = 1;
   int K = 1;
   double alpha[2] = {0, 0};
   double beta = 0;
   double A[] = { -0.535, -0.017 };
   int lda = 1;
   double B[] = { -0.018, -0.544 };
   int ldb = 1;
   double C[] = { -0.8, -0.89 };
   int ldc = 1;
   double C_expected[] = { 0.000000000000e+00, 0.000000000000e+00 };
   cblas_zher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zher2k(case 1389) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zher2k(case 1389) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 113;
   int N = 1;
   int K = 1;
   double alpha[2] = {0, 0};
   double beta = 0;
   double A[] = { -0.619, 0.063 };
   int lda = 1;
   double B[] = { 0.939, 0.705 };
   int ldb = 1;
   double C[] = { 0.829, -0.136 };
   int ldc = 1;
   double C_expected[] = { 0.000000000000e+00, 0.000000000000e+00 };
   cblas_zher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zher2k(case 1390) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zher2k(case 1390) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 113;
   int N = 1;
   int K = 1;
   double alpha[2] = {0, 0};
   double beta = 0;
   double A[] = { -0.089, -0.822 };
   int lda = 1;
   double B[] = { 0.937, 0.159 };
   int ldb = 1;
   double C[] = { 0.977, 0.4 };
   int ldc = 1;
   double C_expected[] = { 0.000000000000e+00, 0.000000000000e+00 };
   cblas_zher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zher2k(case 1391) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zher2k(case 1391) imag");
     };
   };
  };


}
