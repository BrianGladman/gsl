#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_hemm () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int M = 1;
   int N = 2;
   float alpha[2] = {0, 0.1};
   float beta[2] = {0, 0.1};
   float A[] = { -0.126, 0.079 };
   int lda = 1;
   float B[] = { -0.954, -0.059, 0.296, -0.988 };
   int ldb = 2;
   float C[] = { -0.859, -0.731, 0.737, 0.593 };
   int ldc = 2;
   float C_expected[] = { 0.0723566, -0.0738796, -0.0717488, 0.0699704 };
   cblas_chemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "chemm(case 1550) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "chemm(case 1550) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int M = 1;
   int N = 2;
   float alpha[2] = {0, 0.1};
   float beta[2] = {0, 0.1};
   float A[] = { 0.652, 0.584 };
   int lda = 1;
   float B[] = { -0.983, -0.734, -0.422, -0.825 };
   int ldb = 1;
   float C[] = { 0.387, 0.341, -0.734, 0.632 };
   int ldc = 1;
   float C_expected[] = { 0.0137568, -0.0253916, -0.00941, -0.100914 };
   cblas_chemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "chemm(case 1551) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "chemm(case 1551) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int M = 1;
   int N = 2;
   float alpha[2] = {0, 1};
   float beta[2] = {-1, 0};
   float A[] = { 0.78, 0.885, 0.507, 0.765, 0.911, -0.461, 0.707, 0.508 };
   int lda = 2;
   float B[] = { -0.905, 0.633, 0.85, -0.943 };
   int ldb = 2;
   float C[] = { 0.045, -0.237, 0.078, -0.252 };
   int ldc = 2;
   float C_expected[] = { 0.589611, -0.759345, 0.960095, -0.09013 };
   cblas_chemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "chemm(case 1552) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "chemm(case 1552) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int M = 1;
   int N = 2;
   float alpha[2] = {0, 1};
   float beta[2] = {-1, 0};
   float A[] = { 0.947, 0.939, -0.267, -0.819, -0.827, -0.937, 0.991, 0.838 };
   int lda = 2;
   float B[] = { 0.871, -0.988, -0.232, -0.434 };
   int ldb = 1;
   float C[] = { -0.261, 0.927, -0.351, -0.203 };
   int ldc = 1;
   float C_expected[] = { 1.0551, 0.496359, 0.780145, -1.67298 };
   cblas_chemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "chemm(case 1553) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "chemm(case 1553) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int M = 1;
   int N = 2;
   float alpha[2] = {-1, 0};
   float beta[2] = {0, 0};
   float A[] = { -0.593, -0.9 };
   int lda = 1;
   float B[] = { -0.861, 0.747, -0.984, 0.595 };
   int ldb = 2;
   float C[] = { -0.589, -0.671, -0.011, -0.417 };
   int ldc = 2;
   float C_expected[] = { -0.510573, 0.442971, -0.583512, 0.352835 };
   cblas_chemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "chemm(case 1554) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "chemm(case 1554) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int M = 1;
   int N = 2;
   float alpha[2] = {-1, 0};
   float beta[2] = {0, 0};
   float A[] = { -0.79, 0.132 };
   int lda = 1;
   float B[] = { -0.243, -0.12, 0.633, -0.556 };
   int ldb = 1;
   float C[] = { -0.658, -0.74, -0.47, 0.481 };
   int ldc = 1;
   float C_expected[] = { -0.19197, -0.0948, 0.50007, -0.43924 };
   cblas_chemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "chemm(case 1555) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "chemm(case 1555) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int M = 1;
   int N = 2;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {0, 1};
   float A[] = { -0.114, -0.515, -0.513, -0.527, -0.995, 0.986, 0.229, -0.076 };
   int lda = 2;
   float B[] = { 0.084, 0.522, 0.61, 0.694 };
   int ldb = 2;
   float C[] = { 0.802, 0.136, -0.161, -0.364 };
   int ldc = 2;
   float C_expected[] = { 0.269101, 0.716492, 0.237088, 0.0290666 };
   cblas_chemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "chemm(case 1556) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "chemm(case 1556) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int M = 1;
   int N = 2;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {0, 1};
   float A[] = { 0.798, -0.324, -0.693, -0.893, -0.223, 0.749, 0.102, -0.357 };
   int lda = 2;
   float B[] = { -0.572, -0.569, -0.391, -0.938 };
   int ldb = 1;
   float C[] = { 0.152, -0.834, -0.633, -0.473 };
   int ldc = 1;
   float C_expected[] = { 1.08642, -0.113853, 0.234826, -0.48289 };
   cblas_chemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "chemm(case 1557) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "chemm(case 1557) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int M = 1;
   int N = 2;
   double alpha[2] = {0, 0.1};
   double beta[2] = {0, 0.1};
   double A[] = { -0.359, 0.089 };
   int lda = 1;
   double B[] = { -0.451, -0.337, -0.901, -0.871 };
   int ldb = 2;
   double C[] = { 0.729, 0.631, 0.364, 0.246 };
   int ldc = 2;
   double C_expected[] = { -0.0751983, 0.0890909, -0.0558689, 0.0687459 };
   cblas_zhemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zhemm(case 1558) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zhemm(case 1558) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int M = 1;
   int N = 2;
   double alpha[2] = {0, 0.1};
   double beta[2] = {0, 0.1};
   double A[] = { 0.044, -0.496 };
   int lda = 1;
   double B[] = { -0.674, 0.281, 0.366, 0.888 };
   int ldb = 1;
   double C[] = { -0.9, 0.919, 0.857, -0.049 };
   int ldc = 1;
   double C_expected[] = { -0.0931364, -0.0929656, 0.0009928, 0.0873104 };
   cblas_zhemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zhemm(case 1559) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zhemm(case 1559) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int M = 1;
   int N = 2;
   double alpha[2] = {0, 0};
   double beta[2] = {0, 0.1};
   double A[] = { -0.314, 0.115, 0.114, 0.878, 0.961, -0.224, 0.973, 0.771 };
   int lda = 2;
   double B[] = { 0.5, -0.016, -0.5, 0.149 };
   int ldb = 2;
   double C[] = { -0.054, 0.064, 0.02, 0.245 };
   int ldc = 2;
   double C_expected[] = { -0.0064, -0.0054, -0.0245, 0.002 };
   cblas_zhemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zhemm(case 1560) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zhemm(case 1560) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int M = 1;
   int N = 2;
   double alpha[2] = {0, 0};
   double beta[2] = {0, 0.1};
   double A[] = { 0.186, 0.578, 0.797, -0.957, -0.539, -0.969, -0.21, 0.354 };
   int lda = 2;
   double B[] = { 0.641, -0.968, 0.15, -0.569 };
   int ldb = 1;
   double C[] = { -0.556, -0.9, 0.197, 0.31 };
   int ldc = 1;
   double C_expected[] = { 0.09, -0.0556, -0.031, 0.0197 };
   cblas_zhemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zhemm(case 1561) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zhemm(case 1561) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int M = 1;
   int N = 2;
   double alpha[2] = {1, 0};
   double beta[2] = {1, 0};
   double A[] = { 0.323, 0.641 };
   int lda = 1;
   double B[] = { -0.188, 0.091, -0.235, 0.523 };
   int ldb = 2;
   double C[] = { 0.919, 0.806, 0.823, -0.94 };
   int ldc = 2;
   double C_expected[] = { 0.858276, 0.835393, 0.747095, -0.771071 };
   cblas_zhemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zhemm(case 1562) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zhemm(case 1562) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int M = 1;
   int N = 2;
   double alpha[2] = {1, 0};
   double beta[2] = {1, 0};
   double A[] = { -0.688, 0.915 };
   int lda = 1;
   double B[] = { 0.914, -0.204, 0.205, -0.476 };
   int ldb = 1;
   double C[] = { 0.27, -0.628, -0.079, 0.507 };
   int ldc = 1;
   double C_expected[] = { -0.358832, -0.487648, -0.22004, 0.834488 };
   cblas_zhemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zhemm(case 1563) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zhemm(case 1563) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int M = 1;
   int N = 2;
   double alpha[2] = {0, 1};
   double beta[2] = {0, 0.1};
   double A[] = { 0.681, 0.574, -0.425, -0.64, 0.792, 0.661, -0.009, 0.005 };
   int lda = 2;
   double B[] = { -0.221, 0.554, -0.465, -0.95 };
   int ldb = 2;
   double C[] = { 0.331, -0.958, -0.826, -0.972 };
   int ldc = 2;
   double C_expected[] = { 0.778291, 0.142269, -0.496199, 0.112747 };
   cblas_zhemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zhemm(case 1564) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zhemm(case 1564) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int M = 1;
   int N = 2;
   double alpha[2] = {0, 1};
   double beta[2] = {0, 0.1};
   double A[] = { 0.959, 0.34, -0.23, 0.064, 0.516, -0.275, 0.714, 0.899 };
   int lda = 2;
   double B[] = { -0.502, -0.987, -0.134, 0.215 };
   int ldb = 1;
   double C[] = { 0.929, 0.181, -0.16, -0.921 };
   int ldc = 1;
   double C_expected[] = { 0.986459, -0.371458, -0.320548, -0.059384 };
   cblas_zhemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zhemm(case 1565) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zhemm(case 1565) imag");
     };
   };
  };


}
