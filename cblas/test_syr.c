#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_syr () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int uplo = 121;
   int N = 1;
   int lda = 1;
   float alpha = 0;
   float A[] = { -0.301 };
   float X[] = { 0.471 };
   int incX = 1;
   float A_expected[] = { -0.301 };
   cblas_ssyr(order, uplo, N, alpha, X, incX, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[i], A_expected[i], flteps, "ssyr(case 1148)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int N = 1;
   int lda = 1;
   float alpha = 0;
   float A[] = { -0.301 };
   float X[] = { 0.471 };
   int incX = 1;
   float A_expected[] = { -0.301 };
   cblas_ssyr(order, uplo, N, alpha, X, incX, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[i], A_expected[i], flteps, "ssyr(case 1149)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int N = 1;
   int lda = 1;
   float alpha = 0;
   float A[] = { -0.301 };
   float X[] = { 0.471 };
   int incX = 1;
   float A_expected[] = { -0.301 };
   cblas_ssyr(order, uplo, N, alpha, X, incX, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[i], A_expected[i], flteps, "ssyr(case 1150)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int N = 1;
   int lda = 1;
   float alpha = 0;
   float A[] = { -0.301 };
   float X[] = { 0.471 };
   int incX = 1;
   float A_expected[] = { -0.301 };
   cblas_ssyr(order, uplo, N, alpha, X, incX, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[i], A_expected[i], flteps, "ssyr(case 1151)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int N = 1;
   int lda = 1;
   double alpha = -0.3;
   double A[] = { -0.952 };
   double X[] = { 0.659 };
   int incX = 1;
   double A_expected[] = { -1.0822843 };
   cblas_dsyr(order, uplo, N, alpha, X, incX, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[i], A_expected[i], dbleps, "dsyr(case 1152)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int N = 1;
   int lda = 1;
   double alpha = -0.3;
   double A[] = { -0.952 };
   double X[] = { 0.659 };
   int incX = 1;
   double A_expected[] = { -1.0822843 };
   cblas_dsyr(order, uplo, N, alpha, X, incX, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[i], A_expected[i], dbleps, "dsyr(case 1153)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int N = 1;
   int lda = 1;
   double alpha = -0.3;
   double A[] = { -0.952 };
   double X[] = { 0.659 };
   int incX = 1;
   double A_expected[] = { -1.0822843 };
   cblas_dsyr(order, uplo, N, alpha, X, incX, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[i], A_expected[i], dbleps, "dsyr(case 1154)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int N = 1;
   int lda = 1;
   double alpha = -0.3;
   double A[] = { -0.952 };
   double X[] = { 0.659 };
   int incX = 1;
   double A_expected[] = { -1.0822843 };
   cblas_dsyr(order, uplo, N, alpha, X, incX, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[i], A_expected[i], dbleps, "dsyr(case 1155)");
     }
   };
  };


}
