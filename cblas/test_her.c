#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_her () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int uplo = 121;
   int N = 1;
   int lda = 1;
   float alpha = 1;
   float A[] = { 0.232, 0.478 };
   float X[] = { -0.118, -0.933 };
   int incX = 1;
   float A_expected[] = { 1.116413, 0 };
   cblas_cher(order, uplo, N, alpha, X, incX, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[2*i], A_expected[2*i], flteps, "cher(case 1156) real");
       gsl_test_rel(A[2*i+1], A_expected[2*i+1], flteps, "cher(case 1156) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int N = 1;
   int lda = 1;
   float alpha = 1;
   float A[] = { 0.232, 0.478 };
   float X[] = { -0.118, -0.933 };
   int incX = 1;
   float A_expected[] = { 1.116413, 0 };
   cblas_cher(order, uplo, N, alpha, X, incX, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[2*i], A_expected[2*i], flteps, "cher(case 1157) real");
       gsl_test_rel(A[2*i+1], A_expected[2*i+1], flteps, "cher(case 1157) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int N = 1;
   int lda = 1;
   float alpha = 1;
   float A[] = { 0.232, 0.478 };
   float X[] = { -0.118, -0.933 };
   int incX = 1;
   float A_expected[] = { 1.116413, 0 };
   cblas_cher(order, uplo, N, alpha, X, incX, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[2*i], A_expected[2*i], flteps, "cher(case 1158) real");
       gsl_test_rel(A[2*i+1], A_expected[2*i+1], flteps, "cher(case 1158) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int N = 1;
   int lda = 1;
   float alpha = 1;
   float A[] = { 0.232, 0.478 };
   float X[] = { -0.118, -0.933 };
   int incX = 1;
   float A_expected[] = { 1.116413, 0 };
   cblas_cher(order, uplo, N, alpha, X, incX, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[2*i], A_expected[2*i], flteps, "cher(case 1159) real");
       gsl_test_rel(A[2*i+1], A_expected[2*i+1], flteps, "cher(case 1159) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int N = 1;
   int lda = 1;
   double alpha = 0.1;
   double A[] = { -0.923, 0.819 };
   double X[] = { -0.045, 0.926 };
   int incX = 1;
   double A_expected[] = { -0.8370499, 0 };
   cblas_zher(order, uplo, N, alpha, X, incX, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[2*i], A_expected[2*i], dbleps, "zher(case 1160) real");
       gsl_test_rel(A[2*i+1], A_expected[2*i+1], dbleps, "zher(case 1160) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int N = 1;
   int lda = 1;
   double alpha = 0.1;
   double A[] = { -0.923, 0.819 };
   double X[] = { -0.045, 0.926 };
   int incX = 1;
   double A_expected[] = { -0.8370499, 0 };
   cblas_zher(order, uplo, N, alpha, X, incX, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[2*i], A_expected[2*i], dbleps, "zher(case 1161) real");
       gsl_test_rel(A[2*i+1], A_expected[2*i+1], dbleps, "zher(case 1161) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int N = 1;
   int lda = 1;
   double alpha = 0.1;
   double A[] = { -0.923, 0.819 };
   double X[] = { -0.045, 0.926 };
   int incX = 1;
   double A_expected[] = { -0.8370499, 0 };
   cblas_zher(order, uplo, N, alpha, X, incX, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[2*i], A_expected[2*i], dbleps, "zher(case 1162) real");
       gsl_test_rel(A[2*i+1], A_expected[2*i+1], dbleps, "zher(case 1162) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int N = 1;
   int lda = 1;
   double alpha = 0.1;
   double A[] = { -0.923, 0.819 };
   double X[] = { -0.045, 0.926 };
   int incX = 1;
   double A_expected[] = { -0.8370499, 0 };
   cblas_zher(order, uplo, N, alpha, X, incX, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[2*i], A_expected[2*i], dbleps, "zher(case 1163) real");
       gsl_test_rel(A[2*i+1], A_expected[2*i+1], dbleps, "zher(case 1163) imag");
     };
   };
  };


}
