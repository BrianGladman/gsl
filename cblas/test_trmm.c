#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_trmm () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 1;
   int N = 1;
   float alpha[2] = {-0.3, 0.1};
   float A[] = { 0.88, -0.242 };
   int lda = 1;
   float B[] = { -0.413, -0.172 };
   int ldb = 1;
   float B_expected[] = { 0.1266606, -0.0250822 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1392) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1392) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { 0.81, 0.388 };
   int lda = 1;
   float B[] = { -0.02, -0.795 };
   int ldb = 1;
   float B_expected[] = { 0.0855, 0.2365 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1393) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1393) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { 0.031, -0.938 };
   int lda = 1;
   float B[] = { 0.09, -0.847 };
   int ldb = 1;
   float B_expected[] = { 0.2485765, -0.0459665 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1394) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1394) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { -0.478, -0.981 };
   int lda = 1;
   float B[] = { 0.09, -0.286 };
   int ldb = 1;
   float B_expected[] = { 0.0016, 0.0948 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1395) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1395) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { -0.99, -0.495 };
   int lda = 1;
   float B[] = { 0.881, 0.378 };
   int ldb = 1;
   float B_expected[] = { 0.2865555, 0.1745865 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1396) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1396) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { 0.532, 0.746 };
   int lda = 1;
   float B[] = { 0.859, -0.556 };
   int ldb = 1;
   float B_expected[] = { -0.2021, 0.2527 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1397) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1397) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { 0.967, 0.687 };
   int lda = 1;
   float B[] = { 0.673, -0.525 };
   int ldb = 1;
   float B_expected[] = { -0.2989074, 0.1147438 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1398) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1398) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { 0.485, 0.109 };
   int lda = 1;
   float B[] = { -0.024, 0.527 };
   int ldb = 1;
   float B_expected[] = { -0.0455, -0.1605 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1399) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1399) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { 0.596, -0.358 };
   int lda = 1;
   float B[] = { -0.46, 0.558 };
   int ldb = 1;
   float B_expected[] = { -0.027406, -0.156614 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1400) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1400) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { 0.513, -0.201 };
   int lda = 1;
   float B[] = { 0.224, -0.941 };
   int ldb = 1;
   float B_expected[] = { 0.0269, 0.3047 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1401) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1401) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { -0.016, -0.324 };
   int lda = 1;
   float B[] = { -0.159, -0.586 };
   int ldb = 1;
   float B_expected[] = { 0.0501068, -0.0369996 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1402) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1402) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { -0.229, 0.901 };
   int lda = 1;
   float B[] = { 0.411, 0.115 };
   int ldb = 1;
   float B_expected[] = { -0.1348, 0.0066 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1403) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1403) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { 0.252, -0.586 };
   int lda = 1;
   float B[] = { 0.112, 0.042 };
   int ldb = 1;
   float B_expected[] = { -0.010346, 0.021798 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1404) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1404) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { -0.524, 0.768 };
   int lda = 1;
   float B[] = { 0.809, 0.709 };
   int ldb = 1;
   float B_expected[] = { -0.3136, -0.1318 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1405) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1405) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { -0.478, -0.129 };
   int lda = 1;
   float B[] = { 0.7, 0.619 };
   int ldb = 1;
   float B_expected[] = { 0.1150429, 0.0903797 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1406) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1406) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { 0.454, -0.09 };
   int lda = 1;
   float B[] = { -0.778, -0.432 };
   int ldb = 1;
   float B_expected[] = { 0.2766, 0.0518 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1407) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1407) imag");
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
   float alpha[2] = {0, 1};
   float A[] = { 0.948, -0.624 };
   int lda = 1;
   float B[] = { -0.164, 0.522 };
   int ldb = 1;
   float B_expected[] = { -0.597192, 0.170256 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1408) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1408) imag");
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
   float alpha[2] = {0, 1};
   float A[] = { 0.359, 0.622 };
   int lda = 1;
   float B[] = { -0.142, 0.778 };
   int ldb = 1;
   float B_expected[] = { -0.778, -0.142 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1409) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1409) imag");
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
   float alpha[2] = {0, 1};
   float A[] = { -0.282, -0.805 };
   int lda = 1;
   float B[] = { -0.637, -0.757 };
   int ldb = 1;
   float B_expected[] = { -0.726259, -0.429751 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1410) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1410) imag");
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
   float alpha[2] = {0, 1};
   float A[] = { -0.179, -0.996 };
   int lda = 1;
   float B[] = { 0.216, -0.972 };
   int ldb = 1;
   float B_expected[] = { 0.972, 0.216 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1411) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1411) imag");
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
   float alpha[2] = {0, 1};
   float A[] = { -0.082, -0.096 };
   int lda = 1;
   float B[] = { -0.284, -0.414 };
   int ldb = 1;
   float B_expected[] = { -0.061212, -0.016456 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1412) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1412) imag");
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
   float alpha[2] = {0, 1};
   float A[] = { -0.281, 0.913 };
   int lda = 1;
   float B[] = { 0.974, -0.773 };
   int ldb = 1;
   float B_expected[] = { 0.773, 0.974 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1413) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1413) imag");
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
   float alpha[2] = {0, 1};
   float A[] = { 0.745, 0.352 };
   int lda = 1;
   float B[] = { -0.706, 0.522 };
   int ldb = 1;
   float B_expected[] = { -0.140378, -0.709714 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1414) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1414) imag");
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
   float alpha[2] = {0, 1};
   float A[] = { 0.188, 0.235 };
   int lda = 1;
   float B[] = { 0.304, -0.493 };
   int ldb = 1;
   float B_expected[] = { 0.493, 0.304 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1415) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1415) imag");
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
   float alpha[2] = {0, 1};
   float A[] = { -0.664, -0.733 };
   int lda = 1;
   float B[] = { -0.299, -0.686 };
   int ldb = 1;
   float B_expected[] = { -0.674671, -0.304302 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1416) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1416) imag");
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
   float alpha[2] = {0, 1};
   float A[] = { 0.631, 0.152 };
   int lda = 1;
   float B[] = { -0.555, -0.089 };
   int ldb = 1;
   float B_expected[] = { 0.089, -0.555 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1417) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1417) imag");
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
   float alpha[2] = {0, 1};
   float A[] = { -0.129, 0.214 };
   int lda = 1;
   float B[] = { -0.731, -0.332 };
   int ldb = 1;
   float B_expected[] = { 0.113606, 0.165347 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1418) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1418) imag");
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
   float alpha[2] = {0, 1};
   float A[] = { 0.044, 0.279 };
   int lda = 1;
   float B[] = { 0.842, -0.626 };
   int ldb = 1;
   float B_expected[] = { 0.626, 0.842 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1419) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1419) imag");
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
   float alpha[2] = {0, 1};
   float A[] = { -0.33, 0.712 };
   int lda = 1;
   float B[] = { -0.363, -0.768 };
   int ldb = 1;
   float B_expected[] = { 0.005016, 0.666606 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1420) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1420) imag");
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
   float alpha[2] = {0, 1};
   float A[] = { 0.138, -0.178 };
   int lda = 1;
   float B[] = { -0.788, -0.551 };
   int ldb = 1;
   float B_expected[] = { 0.551, -0.788 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1421) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1421) imag");
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
   float alpha[2] = {0, 1};
   float A[] = { -0.034, 0.499 };
   int lda = 1;
   float B[] = { 0.888, -0.199 };
   int ldb = 1;
   float B_expected[] = { -0.449878, 0.069109 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1422) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1422) imag");
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
   float alpha[2] = {0, 1};
   float A[] = { 0.927, -0.864 };
   int lda = 1;
   float B[] = { 0.648, 0.97 };
   int ldb = 1;
   float B_expected[] = { -0.97, 0.648 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1423) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1423) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { 0.664, 0.324 };
   int lda = 1;
   float B[] = { 0.706, 0.483 };
   int ldb = 1;
   float B_expected[] = { -0.1967796, 0.0349372 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1424) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1424) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { 0.518, 0.868 };
   int lda = 1;
   float B[] = { 0.609, -0.803 };
   int ldb = 1;
   float B_expected[] = { -0.1024, 0.3018 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1425) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1425) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { 0.724, 0.812 };
   int lda = 1;
   float B[] = { -0.509, 0.329 };
   int ldb = 1;
   float B_expected[] = { -0.03474, -0.205588 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1426) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1426) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { -0.636, -0.471 };
   int lda = 1;
   float B[] = { 0.241, 0.439 };
   int ldb = 1;
   float B_expected[] = { -0.1162, -0.1076 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1427) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1427) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { 0.531, 0.062 };
   int lda = 1;
   float B[] = { 0.795, 0.131 };
   int ldb = 1;
   float B_expected[] = { -0.1311072, 0.0369454 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1428) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1428) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { 0.537, -0.396 };
   int lda = 1;
   float B[] = { 0.873, 0.239 };
   int ldb = 1;
   float B_expected[] = { -0.2858, 0.0156 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1429) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1429) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { 0.165, 0.807 };
   int lda = 1;
   float B[] = { 0.081, 0.489 };
   int ldb = 1;
   float B_expected[] = { -0.1239282, 0.0362034 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1430) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1430) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { -0.954, 0.838 };
   int lda = 1;
   float B[] = { 0.372, 0.804 };
   int ldb = 1;
   float B_expected[] = { -0.192, -0.204 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1431) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1431) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { 0.429, 0.864 };
   int lda = 1;
   float B[] = { 0.15, -0.438 };
   int ldb = 1;
   float B_expected[] = { 0.1259748, 0.0638424 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1432) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1432) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { -0.228, 0.832 };
   int lda = 1;
   float B[] = { -0.074, 0.735 };
   int ldb = 1;
   float B_expected[] = { -0.0513, -0.2279 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1433) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1433) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { 0.332, 0.198 };
   int lda = 1;
   float B[] = { 0.794, -0.319 };
   int ldb = 1;
   float B_expected[] = { -0.0338218, 0.0989806 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1434) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1434) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { 0.746, -0.361 };
   int lda = 1;
   float B[] = { -0.522, 0.578 };
   int ldb = 1;
   float B_expected[] = { 0.0988, -0.2256 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1435) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1435) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { -0.969, 0.961 };
   int lda = 1;
   float B[] = { -0.163, 0.536 };
   int ldb = 1;
   float B_expected[] = { -0.1656388, 0.1761266 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1436) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1436) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { -0.352, 0.187 };
   int lda = 1;
   float B[] = { 0.753, -0.997 };
   int ldb = 1;
   float B_expected[] = { -0.1262, 0.3744 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1437) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1437) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { 0.899, -0.583 };
   int lda = 1;
   float B[] = { -0.781, 0.888 };
   int ldb = 1;
   float B_expected[] = { 0.331648, -0.224879 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1438) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1438) imag");
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
   float alpha[2] = {-0.3, 0.1};
   float A[] = { -0.2, 0.319 };
   int lda = 1;
   float B[] = { 0.29, 0.729 };
   int ldb = 1;
   float B_expected[] = { -0.1599, -0.1897 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1439) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1439) imag");
     };
   };
  };


}
