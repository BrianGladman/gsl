#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_syr2 () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int uplo = 121;
   int N = 1;
   int lda = 1;
   float alpha = -0.3;
   float A[] = { -0.269 };
   float X[] = { 0.455 };
   int incX = 1;
   float Y[] = { 0.895 };
   int incY = -1;
   float A_expected[] = { -0.513335 };
   cblas_ssyr2(order, uplo, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[i], A_expected[i], flteps, "ssyr2(case 1180)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int N = 1;
   int lda = 1;
   float alpha = -0.3;
   float A[] = { -0.269 };
   float X[] = { 0.455 };
   int incX = 1;
   float Y[] = { 0.895 };
   int incY = -1;
   float A_expected[] = { -0.513335 };
   cblas_ssyr2(order, uplo, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[i], A_expected[i], flteps, "ssyr2(case 1181)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int N = 1;
   int lda = 1;
   float alpha = -0.3;
   float A[] = { -0.269 };
   float X[] = { 0.455 };
   int incX = 1;
   float Y[] = { 0.895 };
   int incY = -1;
   float A_expected[] = { -0.513335 };
   cblas_ssyr2(order, uplo, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[i], A_expected[i], flteps, "ssyr2(case 1182)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int N = 1;
   int lda = 1;
   float alpha = -0.3;
   float A[] = { -0.269 };
   float X[] = { 0.455 };
   int incX = 1;
   float Y[] = { 0.895 };
   int incY = -1;
   float A_expected[] = { -0.513335 };
   cblas_ssyr2(order, uplo, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[i], A_expected[i], flteps, "ssyr2(case 1183)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int N = 1;
   int lda = 1;
   double alpha = 0;
   double A[] = { -0.465 };
   double X[] = { -0.628 };
   int incX = 1;
   double Y[] = { -0.116 };
   int incY = -1;
   double A_expected[] = { -0.465 };
   cblas_dsyr2(order, uplo, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[i], A_expected[i], dbleps, "dsyr2(case 1184)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int N = 1;
   int lda = 1;
   double alpha = 0;
   double A[] = { -0.465 };
   double X[] = { -0.628 };
   int incX = 1;
   double Y[] = { -0.116 };
   int incY = -1;
   double A_expected[] = { -0.465 };
   cblas_dsyr2(order, uplo, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[i], A_expected[i], dbleps, "dsyr2(case 1185)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int N = 1;
   int lda = 1;
   double alpha = 0;
   double A[] = { -0.465 };
   double X[] = { -0.628 };
   int incX = 1;
   double Y[] = { -0.116 };
   int incY = -1;
   double A_expected[] = { -0.465 };
   cblas_dsyr2(order, uplo, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[i], A_expected[i], dbleps, "dsyr2(case 1186)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int N = 1;
   int lda = 1;
   double alpha = 0;
   double A[] = { -0.465 };
   double X[] = { -0.628 };
   int incX = 1;
   double Y[] = { -0.116 };
   int incY = -1;
   double A_expected[] = { -0.465 };
   cblas_dsyr2(order, uplo, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[i], A_expected[i], dbleps, "dsyr2(case 1187)");
     }
   };
  };


}
