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
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.824, -0.161 };
   int lda = 1;
   float B[] = { 0.622, -0.482 };
   int ldb = 1;
   float B_expected[] = { -0.070550732944, -0.0617006456262 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1440) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1440) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.882, 0.87 };
   int lda = 1;
   float B[] = { 0.149, 0.263 };
   int ldb = 1;
   float B_expected[] = { -0.0263, 0.0149 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1441) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1441) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.907, -0.527 };
   int lda = 1;
   float B[] = { 0.187, -0.715 };
   int ldb = 1;
   float B_expected[] = { -0.0678906702969, 0.018829529489 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1442) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1442) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.136, 0.964 };
   int lda = 1;
   float B[] = { -0.915, -0.986 };
   int ldb = 1;
   float B_expected[] = { 0.0986, -0.0915 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1443) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1443) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { 0.161, -0.263 };
   int lda = 1;
   float B[] = { -0.569, 0.603 };
   int ldb = 1;
   float B_expected[] = { 0.055278157535, -0.263117047008 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1444) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1444) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { 0.712, 0.799 };
   int lda = 1;
   float B[] = { -0.249, -0.799 };
   int ldb = 1;
   float B_expected[] = { 0.0799, -0.0249 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1445) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1445) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.698, -0.256 };
   int lda = 1;
   float B[] = { 0.686, 0.758 };
   int ldb = 1;
   float B_expected[] = { 0.0639483301371, -0.12173463111 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1446) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1446) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.266, -0.692 };
   int lda = 1;
   float B[] = { 0.483, -0.224 };
   int ldb = 1;
   float B_expected[] = { 0.0224, 0.0483 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1447) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1447) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { 0.832, -0.815 };
   int lda = 1;
   float B[] = { 0.963, -0.734 };
   int ldb = 1;
   float B_expected[] = { -0.012839185255, 0.103168346174 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1448) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1448) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { 0.545, -0.434 };
   int lda = 1;
   float B[] = { -0.943, -0.253 };
   int ldb = 1;
   float B_expected[] = { 0.0253, -0.0943 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1449) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1449) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.339, -0.919 };
   int lda = 1;
   float B[] = { 0.23, 0.906 };
   int ldb = 1;
   float B_expected[] = { 0.00998080214116, -0.0949037084594 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1450) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1450) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { 0.751, 0.752 };
   int lda = 1;
   float B[] = { 0.535, 0.119 };
   int ldb = 1;
   float B_expected[] = { -0.0119, 0.0535 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1451) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1451) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { 0.597, 0.793 };
   int lda = 1;
   float B[] = { -0.211, 0.232 };
   int ldb = 1;
   float B_expected[] = { -0.0310402960443, 0.00588769642063 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1452) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1452) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { 0.044, -0.469 };
   int lda = 1;
   float B[] = { 0.526, 0.053 };
   int ldb = 1;
   float B_expected[] = { -0.0053, 0.0526 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1453) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1453) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { 0.956, 0.152 };
   int lda = 1;
   float B[] = { 0.985, 0.634 };
   int ldb = 1;
   float B_expected[] = { -0.0487048578503, 0.110777341416 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1454) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1454) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.657, -0.293 };
   int lda = 1;
   float B[] = { 0.268, 0.717 };
   int ldb = 1;
   float B_expected[] = { -0.0717, 0.0268 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1455) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1455) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0};
   float A[] = { -0.168, -0.181 };
   int lda = 1;
   float B[] = { 0.569, -0.497 };
   int ldb = 1;
   float B_expected[] = { -0.000000000000e+00, 0.000000000000e+00 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1456) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1456) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0};
   float A[] = { -0.913, 0.365 };
   int lda = 1;
   float B[] = { 0.91, -0.047 };
   int ldb = 1;
   float B_expected[] = { 0.000000000000e+00, 0.000000000000e+00 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1457) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1457) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0};
   float A[] = { -0.494, -0.361 };
   int lda = 1;
   float B[] = { -0.627, 0.885 };
   int ldb = 1;
   float B_expected[] = { 0.000000000000e+00, 0.000000000000e+00 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1458) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1458) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0};
   float A[] = { 0.998, 0.467 };
   int lda = 1;
   float B[] = { -0.055, 0.796 };
   int ldb = 1;
   float B_expected[] = { -0.000000000000e+00, 0.000000000000e+00 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1459) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1459) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0};
   float A[] = { 0.071, -0.154 };
   int lda = 1;
   float B[] = { -0.523, 0.234 };
   int ldb = 1;
   float B_expected[] = { -0.000000000000e+00, 0.000000000000e+00 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1460) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1460) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0};
   float A[] = { 0.074, 0.752 };
   int lda = 1;
   float B[] = { -0.32, 0.436 };
   int ldb = 1;
   float B_expected[] = { -0.000000000000e+00, 0.000000000000e+00 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1461) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1461) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0};
   float A[] = { 0.631, 0.221 };
   int lda = 1;
   float B[] = { 0.209, 0.394 };
   int ldb = 1;
   float B_expected[] = { 0.000000000000e+00, 0.000000000000e+00 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1462) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1462) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0};
   float A[] = { -0.833, 0.672 };
   int lda = 1;
   float B[] = { 0.2, -0.482 };
   int ldb = 1;
   float B_expected[] = { 0.000000000000e+00, 0.000000000000e+00 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1463) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1463) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0};
   float A[] = { 0.111, 0.376 };
   int lda = 1;
   float B[] = { 0.967, 0.708 };
   int ldb = 1;
   float B_expected[] = { 0.000000000000e+00, 0.000000000000e+00 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1464) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1464) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0};
   float A[] = { 0.372, 0.129 };
   int lda = 1;
   float B[] = { -0.997, 0.349 };
   int ldb = 1;
   float B_expected[] = { -0.000000000000e+00, 0.000000000000e+00 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1465) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1465) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0};
   float A[] = { -0.912, 0.606 };
   int lda = 1;
   float B[] = { -0.672, 0.695 };
   int ldb = 1;
   float B_expected[] = { 0.000000000000e+00, 0.000000000000e+00 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1466) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1466) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0};
   float A[] = { -0.563, -0.295 };
   int lda = 1;
   float B[] = { -0.741, 0.543 };
   int ldb = 1;
   float B_expected[] = { -0.000000000000e+00, 0.000000000000e+00 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1467) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1467) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0};
   float A[] = { 0.754, -0.335 };
   int lda = 1;
   float B[] = { -0.95, -0.261 };
   int ldb = 1;
   float B_expected[] = { -0.000000000000e+00, 0.000000000000e+00 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1468) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1468) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0};
   float A[] = { 0.422, 0.827 };
   int lda = 1;
   float B[] = { -0.4, 0.333 };
   int ldb = 1;
   float B_expected[] = { -0.000000000000e+00, 0.000000000000e+00 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1469) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1469) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0};
   float A[] = { -0.319, -0.862 };
   int lda = 1;
   float B[] = { 0.373, -0.9 };
   int ldb = 1;
   float B_expected[] = { -0.000000000000e+00, 0.000000000000e+00 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1470) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1470) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0};
   float A[] = { 0.151, 0.603 };
   int lda = 1;
   float B[] = { -0.265, 0.892 };
   int ldb = 1;
   float B_expected[] = { -0.000000000000e+00, 0.000000000000e+00 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1471) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1471) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 113;
   int diag = 131;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { 0.675, 0.463 };
   int lda = 1;
   float B[] = { -0.36, -0.978 };
   int ldb = 1;
   float B_expected[] = { 0.123408567838, 0.0313158028281 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1472) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1472) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 113;
   int diag = 132;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { 0.816, 0.831 };
   int lda = 1;
   float B[] = { 0.851, 0.766 };
   int ldb = 1;
   float B_expected[] = { -0.0766, 0.0851 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1473) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1473) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 113;
   int diag = 131;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { 0.186, 0.144 };
   int lda = 1;
   float B[] = { -0.954, 0.485 };
   int ldb = 1;
   float B_expected[] = { 0.0852418130557, -0.446909564086 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1474) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1474) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 113;
   int diag = 132;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.912, 0.756 };
   int lda = 1;
   float B[] = { -0.458, -0.213 };
   int ldb = 1;
   float B_expected[] = { 0.0213, -0.0458 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1475) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1475) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 113;
   int diag = 131;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { 0.966, -0.013 };
   int lda = 1;
   float B[] = { 0.54, 0.024 };
   int ldb = 1;
   float B_expected[] = { -0.00173187260601, 0.0559239278922 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1476) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1476) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 113;
   int diag = 132;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { 0.476, -0.655 };
   int lda = 1;
   float B[] = { 0.773, -0.231 };
   int ldb = 1;
   float B_expected[] = { 0.0231, 0.0773 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1477) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1477) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 113;
   int diag = 131;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.684, 0.018 };
   int lda = 1;
   float B[] = { -0.543, -0.353 };
   int ldb = 1;
   float B_expected[] = { -0.0494848135333, 0.0806881968474 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1478) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1478) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 113;
   int diag = 132;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.757, -0.479 };
   int lda = 1;
   float B[] = { -0.688, -0.094 };
   int ldb = 1;
   float B_expected[] = { 0.0094, -0.0688 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1479) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1479) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 113;
   int diag = 131;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { 0.144, -0.77 };
   int lda = 1;
   float B[] = { 0, -0.168 };
   int ldb = 1;
   float B_expected[] = { 0.00394240233624, -0.0210809013813 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1480) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1480) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 113;
   int diag = 132;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { 0.574, 0.473 };
   int lda = 1;
   float B[] = { 0.909, -0.938 };
   int ldb = 1;
   float B_expected[] = { 0.0938, 0.0909 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1481) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1481) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 113;
   int diag = 131;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.836, -0.096 };
   int lda = 1;
   float B[] = { 0.064, 0.858 };
   int ldb = 1;
   float B_expected[] = { 0.102163499559, 0.00407619133696 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1482) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1482) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 113;
   int diag = 132;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.376, -0.682 };
   int lda = 1;
   float B[] = { 0.114, -0.851 };
   int ldb = 1;
   float B_expected[] = { 0.0851, 0.0114 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1483) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1483) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 113;
   int diag = 131;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { 0.614, -0.964 };
   int lda = 1;
   float B[] = { 0.66, -0.663 };
   int ldb = 1;
   float B_expected[] = { 0.079868972634, -0.017905031953 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1484) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1484) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 113;
   int diag = 132;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.344, -0.532 };
   int lda = 1;
   float B[] = { -0.113, -0.856 };
   int ldb = 1;
   float B_expected[] = { 0.0856, -0.0113 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1485) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1485) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 113;
   int diag = 131;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { 0.63, 0.041 };
   int lda = 1;
   float B[] = { -0.606, -0.279 };
   int ldb = 1;
   float B_expected[] = { 0.0503325547379, -0.0929148654853 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1486) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1486) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 113;
   int diag = 132;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.388, 0.222 };
   int lda = 1;
   float B[] = { 0.663, -0.134 };
   int ldb = 1;
   float B_expected[] = { 0.0134, 0.0663 };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1487) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1487) imag");
     };
   };
  };


}
