#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_her2 () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int uplo = 121;
   int N = 1;
   int lda = 1;
   float alpha[2] = {0, 0};
   float A[] = { 0.395, -0.546 };
   float X[] = { 0.14, 0.165 };
   int incX = 1;
   float Y[] = { 0.499, -0.305 };
   int incY = -1;
   float A_expected[] = { 0.395, -0.546 };
   cblas_cher2(order, uplo, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[2*i], A_expected[2*i], flteps, "cher2(case 1196) real");
       gsl_test_rel(A[2*i+1], A_expected[2*i+1], flteps, "cher2(case 1196) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int N = 1;
   int lda = 1;
   float alpha[2] = {0, 0};
   float A[] = { 0.395, -0.546 };
   float X[] = { 0.14, 0.165 };
   int incX = 1;
   float Y[] = { 0.499, -0.305 };
   int incY = -1;
   float A_expected[] = { 0.395, -0.546 };
   cblas_cher2(order, uplo, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[2*i], A_expected[2*i], flteps, "cher2(case 1197) real");
       gsl_test_rel(A[2*i+1], A_expected[2*i+1], flteps, "cher2(case 1197) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int N = 1;
   int lda = 1;
   float alpha[2] = {0, 0};
   float A[] = { 0.395, -0.546 };
   float X[] = { 0.14, 0.165 };
   int incX = 1;
   float Y[] = { 0.499, -0.305 };
   int incY = -1;
   float A_expected[] = { 0.395, -0.546 };
   cblas_cher2(order, uplo, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[2*i], A_expected[2*i], flteps, "cher2(case 1198) real");
       gsl_test_rel(A[2*i+1], A_expected[2*i+1], flteps, "cher2(case 1198) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int N = 1;
   int lda = 1;
   float alpha[2] = {0, 0};
   float A[] = { 0.395, -0.546 };
   float X[] = { 0.14, 0.165 };
   int incX = 1;
   float Y[] = { 0.499, -0.305 };
   int incY = -1;
   float A_expected[] = { 0.395, -0.546 };
   cblas_cher2(order, uplo, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[2*i], A_expected[2*i], flteps, "cher2(case 1199) real");
       gsl_test_rel(A[2*i+1], A_expected[2*i+1], flteps, "cher2(case 1199) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int N = 1;
   int lda = 1;
   double alpha[2] = {1, 0};
   double A[] = { -0.674, -0.552 };
   double X[] = { -0.536, -0.027 };
   int incX = 1;
   double Y[] = { 0.459, -0.968 };
   int incY = -1;
   double A_expected[] = { -1.113776, 0 };
   cblas_zher2(order, uplo, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[2*i], A_expected[2*i], dbleps, "zher2(case 1200) real");
       gsl_test_rel(A[2*i+1], A_expected[2*i+1], dbleps, "zher2(case 1200) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int N = 1;
   int lda = 1;
   double alpha[2] = {1, 0};
   double A[] = { -0.674, -0.552 };
   double X[] = { -0.536, -0.027 };
   int incX = 1;
   double Y[] = { 0.459, -0.968 };
   int incY = -1;
   double A_expected[] = { -1.113776, 0 };
   cblas_zher2(order, uplo, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[2*i], A_expected[2*i], dbleps, "zher2(case 1201) real");
       gsl_test_rel(A[2*i+1], A_expected[2*i+1], dbleps, "zher2(case 1201) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int N = 1;
   int lda = 1;
   double alpha[2] = {1, 0};
   double A[] = { -0.674, -0.552 };
   double X[] = { -0.536, -0.027 };
   int incX = 1;
   double Y[] = { 0.459, -0.968 };
   int incY = -1;
   double A_expected[] = { -1.113776, 0 };
   cblas_zher2(order, uplo, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[2*i], A_expected[2*i], dbleps, "zher2(case 1202) real");
       gsl_test_rel(A[2*i+1], A_expected[2*i+1], dbleps, "zher2(case 1202) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int N = 1;
   int lda = 1;
   double alpha[2] = {1, 0};
   double A[] = { -0.674, -0.552 };
   double X[] = { -0.536, -0.027 };
   int incX = 1;
   double Y[] = { 0.459, -0.968 };
   int incY = -1;
   double A_expected[] = { -1.113776, 0 };
   cblas_zher2(order, uplo, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[2*i], A_expected[2*i], dbleps, "zher2(case 1203) real");
       gsl_test_rel(A[2*i+1], A_expected[2*i+1], dbleps, "zher2(case 1203) imag");
     };
   };
  };


}
