#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_sbmv () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int uplo = 121;
   float alpha = -1;
   float beta = -0.3;
   int N = 2;
   int k = 1;
   int lda = 2;
   float A[] = { -0.841, -0.167, -0.847, -0.028 };
   float X[] = { -0.09, 0.589 };
   int incX = 1;
   float Y[] = { -0.904, 0.307 };
   int incY = -1;
   float y_expected[] = { 0.272662, -0.069967 };
   cblas_ssbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "ssbmv(case 848)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   float alpha = -1;
   float beta = -0.3;
   int N = 2;
   int k = 1;
   int lda = 2;
   float A[] = { -0.841, -0.167, -0.847, -0.028 };
   float X[] = { -0.09, 0.589 };
   int incX = 1;
   float Y[] = { -0.904, 0.307 };
   int incY = -1;
   float y_expected[] = { 0.272662, -0.069967 };
   cblas_ssbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "ssbmv(case 849)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   float alpha = -1;
   float beta = -0.3;
   int N = 2;
   int k = 1;
   int lda = 2;
   float A[] = { -0.841, -0.167, -0.847, -0.028 };
   float X[] = { -0.09, 0.589 };
   int incX = 1;
   float Y[] = { -0.904, 0.307 };
   int incY = -1;
   float y_expected[] = { 0.293333, 0.331093 };
   cblas_ssbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "ssbmv(case 850)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   float alpha = -1;
   float beta = -0.3;
   int N = 2;
   int k = 1;
   int lda = 2;
   float A[] = { -0.841, -0.167, -0.847, -0.028 };
   float X[] = { -0.09, 0.589 };
   int incX = 1;
   float Y[] = { -0.904, 0.307 };
   int incY = -1;
   float y_expected[] = { 0.293333, 0.331093 };
   cblas_ssbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "ssbmv(case 851)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   float alpha = -1;
   float beta = -0.3;
   int N = 2;
   int k = 1;
   int lda = 2;
   float A[] = { -0.841, -0.167, -0.847, -0.028 };
   float X[] = { -0.09, 0.589 };
   int incX = 1;
   float Y[] = { -0.904, 0.307 };
   int incY = -1;
   float y_expected[] = { 0.211462, 0.391753 };
   cblas_ssbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "ssbmv(case 852)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   float alpha = -1;
   float beta = -0.3;
   int N = 2;
   int k = 1;
   int lda = 2;
   float A[] = { -0.841, -0.167, -0.847, -0.028 };
   float X[] = { -0.09, 0.589 };
   int incX = 1;
   float Y[] = { -0.904, 0.307 };
   int incY = -1;
   float y_expected[] = { 0.211462, 0.391753 };
   cblas_ssbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "ssbmv(case 853)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   float alpha = -1;
   float beta = -0.3;
   int N = 2;
   int k = 1;
   int lda = 2;
   float A[] = { -0.841, -0.167, -0.847, -0.028 };
   float X[] = { -0.09, 0.589 };
   int incX = 1;
   float Y[] = { -0.904, 0.307 };
   int incY = -1;
   float y_expected[] = { 0.755053, -0.069427 };
   cblas_ssbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "ssbmv(case 854)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   float alpha = -1;
   float beta = -0.3;
   int N = 2;
   int k = 1;
   int lda = 2;
   float A[] = { -0.841, -0.167, -0.847, -0.028 };
   float X[] = { -0.09, 0.589 };
   int incX = 1;
   float Y[] = { -0.904, 0.307 };
   int incY = -1;
   float y_expected[] = { 0.755053, -0.069427 };
   cblas_ssbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "ssbmv(case 855)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   double alpha = -1;
   double beta = -1;
   int N = 2;
   int k = 1;
   int lda = 2;
   double A[] = { -0.655, -0.327, -0.331, 0.299 };
   double X[] = { 0.369, -0.281 };
   int incX = 1;
   double Y[] = { -0.501, -0.218 };
   int incY = -1;
   double y_expected[] = { 0.705682, 0.248252 };
   cblas_dsbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dsbmv(case 856)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   double alpha = -1;
   double beta = -1;
   int N = 2;
   int k = 1;
   int lda = 2;
   double A[] = { -0.655, -0.327, -0.331, 0.299 };
   double X[] = { 0.369, -0.281 };
   int incX = 1;
   double Y[] = { -0.501, -0.218 };
   int incY = -1;
   double y_expected[] = { 0.705682, 0.248252 };
   cblas_dsbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dsbmv(case 857)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   double alpha = -1;
   double beta = -1;
   int N = 2;
   int k = 1;
   int lda = 2;
   double A[] = { -0.655, -0.327, -0.331, 0.299 };
   double X[] = { 0.369, -0.281 };
   int incX = 1;
   double Y[] = { -0.501, -0.218 };
   int incY = -1;
   double y_expected[] = { 0.531252, 0.366684 };
   cblas_dsbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dsbmv(case 858)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   double alpha = -1;
   double beta = -1;
   int N = 2;
   int k = 1;
   int lda = 2;
   double A[] = { -0.655, -0.327, -0.331, 0.299 };
   double X[] = { 0.369, -0.281 };
   int incX = 1;
   double Y[] = { -0.501, -0.218 };
   int incY = -1;
   double y_expected[] = { 0.531252, 0.366684 };
   cblas_dsbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dsbmv(case 859)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   double alpha = -1;
   double beta = -1;
   int N = 2;
   int k = 1;
   int lda = 2;
   double A[] = { -0.655, -0.327, -0.331, 0.299 };
   double X[] = { 0.369, -0.281 };
   int incX = 1;
   double Y[] = { -0.501, -0.218 };
   int incY = -1;
   double y_expected[] = { 0.707158, 0.245652 };
   cblas_dsbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dsbmv(case 860)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   double alpha = -1;
   double beta = -1;
   int N = 2;
   int k = 1;
   int lda = 2;
   double A[] = { -0.655, -0.327, -0.331, 0.299 };
   double X[] = { 0.369, -0.281 };
   int incX = 1;
   double Y[] = { -0.501, -0.218 };
   int incY = -1;
   double y_expected[] = { 0.707158, 0.245652 };
   cblas_dsbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dsbmv(case 861)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   double alpha = -1;
   double beta = -1;
   int N = 2;
   int k = 1;
   int lda = 2;
   double A[] = { -0.655, -0.327, -0.331, 0.299 };
   double X[] = { 0.369, -0.281 };
   int incX = 1;
   double Y[] = { -0.501, -0.218 };
   int incY = -1;
   double y_expected[] = { 0.528652, 0.367808 };
   cblas_dsbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dsbmv(case 862)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   double alpha = -1;
   double beta = -1;
   int N = 2;
   int k = 1;
   int lda = 2;
   double A[] = { -0.655, -0.327, -0.331, 0.299 };
   double X[] = { 0.369, -0.281 };
   int incX = 1;
   double Y[] = { -0.501, -0.218 };
   int incY = -1;
   double y_expected[] = { 0.528652, 0.367808 };
   cblas_dsbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dsbmv(case 863)");
     }
   };
  };


}
