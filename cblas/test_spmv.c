#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_spmv () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int uplo = 121;
   float alpha = 0.1;
   float beta = 0;
   int N = 1;
   float A[] = { 0.15 };
   float X[] = { 0.918 };
   int incX = 1;
   float Y[] = { -0.194 };
   int incY = -1;
   float y_expected[] = { 0.01377 };
   cblas_sspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sspmv(case 880)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   float alpha = 0.1;
   float beta = 0;
   int N = 1;
   float A[] = { 0.15 };
   float X[] = { 0.918 };
   int incX = 1;
   float Y[] = { -0.194 };
   int incY = -1;
   float y_expected[] = { 0.01377 };
   cblas_sspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sspmv(case 881)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   float alpha = 0.1;
   float beta = 0;
   int N = 1;
   float A[] = { 0.15 };
   float X[] = { 0.918 };
   int incX = 1;
   float Y[] = { -0.194 };
   int incY = -1;
   float y_expected[] = { 0.01377 };
   cblas_sspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sspmv(case 882)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   float alpha = 0.1;
   float beta = 0;
   int N = 1;
   float A[] = { 0.15 };
   float X[] = { 0.918 };
   int incX = 1;
   float Y[] = { -0.194 };
   int incY = -1;
   float y_expected[] = { 0.01377 };
   cblas_sspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sspmv(case 883)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   float alpha = 0.1;
   float beta = 0;
   int N = 1;
   float A[] = { 0.15 };
   float X[] = { 0.918 };
   int incX = 1;
   float Y[] = { -0.194 };
   int incY = -1;
   float y_expected[] = { 0.01377 };
   cblas_sspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sspmv(case 884)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   float alpha = 0.1;
   float beta = 0;
   int N = 1;
   float A[] = { 0.15 };
   float X[] = { 0.918 };
   int incX = 1;
   float Y[] = { -0.194 };
   int incY = -1;
   float y_expected[] = { 0.01377 };
   cblas_sspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sspmv(case 885)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   float alpha = 0.1;
   float beta = 0;
   int N = 1;
   float A[] = { 0.15 };
   float X[] = { 0.918 };
   int incX = 1;
   float Y[] = { -0.194 };
   int incY = -1;
   float y_expected[] = { 0.01377 };
   cblas_sspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sspmv(case 886)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   float alpha = 0.1;
   float beta = 0;
   int N = 1;
   float A[] = { 0.15 };
   float X[] = { 0.918 };
   int incX = 1;
   float Y[] = { -0.194 };
   int incY = -1;
   float y_expected[] = { 0.01377 };
   cblas_sspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sspmv(case 887)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   double alpha = 1;
   double beta = -0.3;
   int N = 1;
   double A[] = { 0.284 };
   double X[] = { 0.29 };
   int incX = 1;
   double Y[] = { -0.278 };
   int incY = -1;
   double y_expected[] = { 0.16576 };
   cblas_dspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dspmv(case 888)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   double alpha = 1;
   double beta = -0.3;
   int N = 1;
   double A[] = { 0.284 };
   double X[] = { 0.29 };
   int incX = 1;
   double Y[] = { -0.278 };
   int incY = -1;
   double y_expected[] = { 0.16576 };
   cblas_dspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dspmv(case 889)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   double alpha = 1;
   double beta = -0.3;
   int N = 1;
   double A[] = { 0.284 };
   double X[] = { 0.29 };
   int incX = 1;
   double Y[] = { -0.278 };
   int incY = -1;
   double y_expected[] = { 0.16576 };
   cblas_dspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dspmv(case 890)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   double alpha = 1;
   double beta = -0.3;
   int N = 1;
   double A[] = { 0.284 };
   double X[] = { 0.29 };
   int incX = 1;
   double Y[] = { -0.278 };
   int incY = -1;
   double y_expected[] = { 0.16576 };
   cblas_dspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dspmv(case 891)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   double alpha = 1;
   double beta = -0.3;
   int N = 1;
   double A[] = { 0.284 };
   double X[] = { 0.29 };
   int incX = 1;
   double Y[] = { -0.278 };
   int incY = -1;
   double y_expected[] = { 0.16576 };
   cblas_dspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dspmv(case 892)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   double alpha = 1;
   double beta = -0.3;
   int N = 1;
   double A[] = { 0.284 };
   double X[] = { 0.29 };
   int incX = 1;
   double Y[] = { -0.278 };
   int incY = -1;
   double y_expected[] = { 0.16576 };
   cblas_dspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dspmv(case 893)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   double alpha = 1;
   double beta = -0.3;
   int N = 1;
   double A[] = { 0.284 };
   double X[] = { 0.29 };
   int incX = 1;
   double Y[] = { -0.278 };
   int incY = -1;
   double y_expected[] = { 0.16576 };
   cblas_dspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dspmv(case 894)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   double alpha = 1;
   double beta = -0.3;
   int N = 1;
   double A[] = { 0.284 };
   double X[] = { 0.29 };
   int incX = 1;
   double Y[] = { -0.278 };
   int incY = -1;
   double y_expected[] = { 0.16576 };
   cblas_dspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dspmv(case 895)");
     }
   };
  };


}
