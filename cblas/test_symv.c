#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_symv () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int uplo = 121;
   float alpha = 1;
   float beta = 0.1;
   int N = 1;
   int lda = 1;
   float A[] = { -0.18 };
   float X[] = { 0.519 };
   int incX = 1;
   float Y[] = { -0.991 };
   int incY = -1;
   float y_expected[] = { -0.19252 };
   cblas_ssymv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "ssymv(case 800)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   float alpha = 1;
   float beta = 0.1;
   int N = 1;
   int lda = 1;
   float A[] = { -0.18 };
   float X[] = { 0.519 };
   int incX = 1;
   float Y[] = { -0.991 };
   int incY = -1;
   float y_expected[] = { -0.19252 };
   cblas_ssymv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "ssymv(case 801)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   float alpha = 1;
   float beta = 0.1;
   int N = 1;
   int lda = 1;
   float A[] = { -0.18 };
   float X[] = { 0.519 };
   int incX = 1;
   float Y[] = { -0.991 };
   int incY = -1;
   float y_expected[] = { -0.19252 };
   cblas_ssymv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "ssymv(case 802)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   float alpha = 1;
   float beta = 0.1;
   int N = 1;
   int lda = 1;
   float A[] = { -0.18 };
   float X[] = { 0.519 };
   int incX = 1;
   float Y[] = { -0.991 };
   int incY = -1;
   float y_expected[] = { -0.19252 };
   cblas_ssymv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "ssymv(case 803)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   float alpha = 1;
   float beta = 0.1;
   int N = 1;
   int lda = 1;
   float A[] = { -0.18 };
   float X[] = { 0.519 };
   int incX = 1;
   float Y[] = { -0.991 };
   int incY = -1;
   float y_expected[] = { -0.19252 };
   cblas_ssymv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "ssymv(case 804)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   float alpha = 1;
   float beta = 0.1;
   int N = 1;
   int lda = 1;
   float A[] = { -0.18 };
   float X[] = { 0.519 };
   int incX = 1;
   float Y[] = { -0.991 };
   int incY = -1;
   float y_expected[] = { -0.19252 };
   cblas_ssymv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "ssymv(case 805)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   float alpha = 1;
   float beta = 0.1;
   int N = 1;
   int lda = 1;
   float A[] = { -0.18 };
   float X[] = { 0.519 };
   int incX = 1;
   float Y[] = { -0.991 };
   int incY = -1;
   float y_expected[] = { -0.19252 };
   cblas_ssymv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "ssymv(case 806)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   float alpha = 1;
   float beta = 0.1;
   int N = 1;
   int lda = 1;
   float A[] = { -0.18 };
   float X[] = { 0.519 };
   int incX = 1;
   float Y[] = { -0.991 };
   int incY = -1;
   float y_expected[] = { -0.19252 };
   cblas_ssymv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "ssymv(case 807)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   double alpha = 0.1;
   double beta = 1;
   int N = 1;
   int lda = 1;
   double A[] = { 0.78 };
   double X[] = { -0.372 };
   int incX = 1;
   double Y[] = { -0.399 };
   int incY = -1;
   double y_expected[] = { -0.428016 };
   cblas_dsymv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dsymv(case 808)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   double alpha = 0.1;
   double beta = 1;
   int N = 1;
   int lda = 1;
   double A[] = { 0.78 };
   double X[] = { -0.372 };
   int incX = 1;
   double Y[] = { -0.399 };
   int incY = -1;
   double y_expected[] = { -0.428016 };
   cblas_dsymv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dsymv(case 809)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   double alpha = 0.1;
   double beta = 1;
   int N = 1;
   int lda = 1;
   double A[] = { 0.78 };
   double X[] = { -0.372 };
   int incX = 1;
   double Y[] = { -0.399 };
   int incY = -1;
   double y_expected[] = { -0.428016 };
   cblas_dsymv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dsymv(case 810)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   double alpha = 0.1;
   double beta = 1;
   int N = 1;
   int lda = 1;
   double A[] = { 0.78 };
   double X[] = { -0.372 };
   int incX = 1;
   double Y[] = { -0.399 };
   int incY = -1;
   double y_expected[] = { -0.428016 };
   cblas_dsymv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dsymv(case 811)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   double alpha = 0.1;
   double beta = 1;
   int N = 1;
   int lda = 1;
   double A[] = { 0.78 };
   double X[] = { -0.372 };
   int incX = 1;
   double Y[] = { -0.399 };
   int incY = -1;
   double y_expected[] = { -0.428016 };
   cblas_dsymv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dsymv(case 812)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   double alpha = 0.1;
   double beta = 1;
   int N = 1;
   int lda = 1;
   double A[] = { 0.78 };
   double X[] = { -0.372 };
   int incX = 1;
   double Y[] = { -0.399 };
   int incY = -1;
   double y_expected[] = { -0.428016 };
   cblas_dsymv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dsymv(case 813)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   double alpha = 0.1;
   double beta = 1;
   int N = 1;
   int lda = 1;
   double A[] = { 0.78 };
   double X[] = { -0.372 };
   int incX = 1;
   double Y[] = { -0.399 };
   int incY = -1;
   double y_expected[] = { -0.428016 };
   cblas_dsymv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dsymv(case 814)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   double alpha = 0.1;
   double beta = 1;
   int N = 1;
   int lda = 1;
   double A[] = { 0.78 };
   double X[] = { -0.372 };
   int incX = 1;
   double Y[] = { -0.399 };
   int incY = -1;
   double y_expected[] = { -0.428016 };
   cblas_dsymv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dsymv(case 815)");
     }
   };
  };


}
