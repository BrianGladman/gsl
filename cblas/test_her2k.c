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
   int K = 2;
   float alpha[2] = {0, 0.1};
   float beta = 0.1;
   float A[] = { 0.531, 0.721, -0.848, 0.826 };
   int lda = 2;
   float B[] = { -0.711, -0.2, -0.92, -0.676 };
   int ldb = 2;
   float C[] = { -0.447, 0.701 };
   int ldc = 1;
   float C_expected[] = { 0.30322, 0 };
   cblas_cher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cher2k(case 1654) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cher2k(case 1654) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 2;
   float alpha[2] = {0, 0.1};
   float beta = 0.1;
   float A[] = { 0.68, 0.079, 0.837, -0.814 };
   int lda = 2;
   float B[] = { -0.986, 0.024, 0.584, -0.248 };
   int ldb = 2;
   float C[] = { 0.477, -0.551 };
   int ldc = 1;
   float C_expected[] = { 0.120103, 0 };
   cblas_cher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cher2k(case 1655) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cher2k(case 1655) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 2;
   float alpha[2] = {0, 0.1};
   float beta = 0.1;
   float A[] = { 0.354, -0.63, -0.85, 0.426 };
   int lda = 1;
   float B[] = { 0.787, -0.228, -0.568, 0.83 };
   int ldb = 1;
   float C[] = { 0.428, -0.388 };
   int ldc = 1;
   float C_expected[] = { 0.0331132, 0 };
   cblas_cher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cher2k(case 1656) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cher2k(case 1656) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 2;
   float alpha[2] = {0, 0.1};
   float beta = 0.1;
   float A[] = { -0.49, 0.224, -0.606, 0.46 };
   int lda = 1;
   float B[] = { -0.191, -0.815, 0.464, 0.066 };
   int ldb = 1;
   float C[] = { 0.302, 0.023 };
   int ldc = 1;
   float C_expected[] = { 0.0679396, 0 };
   cblas_cher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cher2k(case 1657) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cher2k(case 1657) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 113;
   int N = 1;
   int K = 2;
   float alpha[2] = {-1, 0};
   float beta = 0;
   float A[] = { 0.943, 0.075, 0.15, -0.141 };
   int lda = 1;
   float B[] = { -0.962, 0.422, -0.592, -0.789 };
   int ldb = 1;
   float C[] = { 0.728, 0.601 };
   int ldc = 1;
   float C_expected[] = { 1.70613, 0 };
   cblas_cher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cher2k(case 1658) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cher2k(case 1658) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 113;
   int N = 1;
   int K = 2;
   float alpha[2] = {-1, 0};
   float beta = 0;
   float A[] = { -0.93, -0.386, 0.565, 0.141 };
   int lda = 1;
   float B[] = { -0.801, 0.022, 0.558, -0.932 };
   int ldb = 1;
   float C[] = { 0.068, 0.501 };
   int ldc = 1;
   float C_expected[] = { -1.84059, 0 };
   cblas_cher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cher2k(case 1659) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cher2k(case 1659) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 113;
   int N = 1;
   int K = 2;
   float alpha[2] = {-1, 0};
   float beta = 0;
   float A[] = { -0.383, 0.124, 0.458, -0.221 };
   int lda = 2;
   float B[] = { -0.107, 0.199, 0.18, 0.122 };
   int ldb = 2;
   float C[] = { 0.896, -0.874 };
   int ldc = 1;
   float C_expected[] = { -0.24227, 0 };
   cblas_cher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cher2k(case 1660) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cher2k(case 1660) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 113;
   int N = 1;
   int K = 2;
   float alpha[2] = {-1, 0};
   float beta = 0;
   float A[] = { 0.131, 0.692, 0.533, -0.672 };
   int lda = 2;
   float B[] = { -0.435, -0.453, 0.195, -0.579 };
   int ldb = 2;
   float C[] = { -0.547, 0.736 };
   int ldc = 1;
   float C_expected[] = { -0.245124, 0 };
   cblas_cher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cher2k(case 1661) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cher2k(case 1661) imag");
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
   double beta = 0.1;
   double A[] = { 0.972, -0.353, 0.712, -0.53 };
   int lda = 2;
   double B[] = { 0.787, -0.379, 0.889, 0.901 };
   int ldb = 2;
   double C[] = { 0.002, 0.266 };
   int ldc = 1;
   double C_expected[] = { -0.4278924, 0 };
   cblas_zher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zher2k(case 1662) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zher2k(case 1662) imag");
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
   double beta = 0.1;
   double A[] = { -0.36, 0.192, 0.539, 0.198 };
   int lda = 2;
   double B[] = { -0.673, 0.781, 0.792, 0.335 };
   int ldb = 2;
   double C[] = { 0.719, -0.339 };
   int ldc = 1;
   double C_expected[] = { -0.485009, 0 };
   cblas_zher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zher2k(case 1663) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zher2k(case 1663) imag");
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
   double beta = 0.1;
   double A[] = { -0.143, 0.456, 0.677, -0.522 };
   int lda = 1;
   double B[] = { 0.851, 0.196, 0.586, 0.64 };
   int ldb = 1;
   double C[] = { 0.617, 0.118 };
   int ldc = 1;
   double C_expected[] = { 0.1081226, 0 };
   cblas_zher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zher2k(case 1664) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zher2k(case 1664) imag");
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
   double beta = 0.1;
   double A[] = { 0.801, 0.91, 0.376, -0.006 };
   int lda = 1;
   double B[] = { -0.613, -0.758, -0.966, 0.194 };
   int ldb = 1;
   double C[] = { -0.723, -0.765 };
   int ldc = 1;
   double C_expected[] = { 0.8583678, 0 };
   cblas_zher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zher2k(case 1665) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zher2k(case 1665) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 113;
   int N = 1;
   int K = 2;
   double alpha[2] = {-0.3, 0.1};
   double beta = 0.1;
   double A[] = { -0.359, 0.913, 0.738, -0.227 };
   int lda = 1;
   double B[] = { 0.787, 0.745, 0.036, -0.606 };
   int ldb = 1;
   double C[] = { -0.652, -0.281 };
   int ldc = 1;
   double C_expected[] = { -0.1172608, 0 };
   cblas_zher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zher2k(case 1666) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zher2k(case 1666) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 113;
   int N = 1;
   int K = 2;
   double alpha[2] = {-0.3, 0.1};
   double beta = 0.1;
   double A[] = { -0.933, 0.598, 0.952, 0.25 };
   int lda = 1;
   double B[] = { -0.508, -0.461, -0.727, 0.162 };
   int ldb = 1;
   double C[] = { 0.215, 0.943 };
   int ldc = 1;
   double C_expected[] = { 0.0795166, 0 };
   cblas_zher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zher2k(case 1667) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zher2k(case 1667) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 113;
   int N = 1;
   int K = 2;
   double alpha[2] = {-0.3, 0.1};
   double beta = 0.1;
   double A[] = { -0.735, 0.372, -0.251, -0.168 };
   int lda = 2;
   double B[] = { 0.217, 0.863, -0.179, -0.057 };
   int ldb = 2;
   double C[] = { 0.579, -0.305 };
   int ldc = 1;
   double C_expected[] = { 0.0744312, 0 };
   cblas_zher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zher2k(case 1668) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zher2k(case 1668) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 113;
   int N = 1;
   int K = 2;
   double alpha[2] = {-0.3, 0.1};
   double beta = 0.1;
   double A[] = { -0.587, -0.994, -0.625, 0.681 };
   int lda = 2;
   double B[] = { -0.577, -0.014, -0.434, 0.204 };
   int ldb = 2;
   double C[] = { 0.256, 0.093 };
   int ldc = 1;
   double C_expected[] = { -0.3526202, 0 };
   cblas_zher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zher2k(case 1669) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zher2k(case 1669) imag");
     };
   };
  };


}
