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
   float alpha = 1;
   float beta = 0;
   int N = 3;
   int k = 1;
   int lda = 3;
   float A[] = { 0.627, -0.312, 0.031, 0.308, 0.323, -0.578, 0.797, 0.545, -0.476 };
   float X[] = { -0.542, 0.606, 0.727 };
   int incX = -1;
   float Y[] = { 0.755, 0.268, -0.99 };
   int incY = -1;
   float y_expected[] = { -0.236236, -0.215242, 0.266757 };
   cblas_ssbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "ssbmv(case 1102)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   float alpha = 1;
   float beta = 0;
   int N = 3;
   int k = 1;
   int lda = 3;
   float A[] = { 0.627, -0.312, 0.031, 0.308, 0.323, -0.578, 0.797, 0.545, -0.476 };
   float X[] = { -0.542, 0.606, 0.727 };
   int incX = -1;
   float Y[] = { 0.755, 0.268, -0.99 };
   int incY = -1;
   float y_expected[] = { -0.236236, -0.215242, 0.266757 };
   cblas_ssbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "ssbmv(case 1103)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   float alpha = 1;
   float beta = 0;
   int N = 3;
   int k = 1;
   int lda = 3;
   float A[] = { 0.627, -0.312, 0.031, 0.308, 0.323, -0.578, 0.797, 0.545, -0.476 };
   float X[] = { -0.542, 0.606, 0.727 };
   int incX = -1;
   float Y[] = { 0.755, 0.268, -0.99 };
   int incY = -1;
   float y_expected[] = { 0.187592, -0.01232, -0.040176 };
   cblas_ssbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "ssbmv(case 1104)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   float alpha = 1;
   float beta = 0;
   int N = 3;
   int k = 1;
   int lda = 3;
   float A[] = { 0.627, -0.312, 0.031, 0.308, 0.323, -0.578, 0.797, 0.545, -0.476 };
   float X[] = { -0.542, 0.606, 0.727 };
   int incX = -1;
   float Y[] = { 0.755, 0.268, -0.99 };
   int incY = -1;
   float y_expected[] = { 0.187592, -0.01232, -0.040176 };
   cblas_ssbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "ssbmv(case 1105)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   float alpha = 1;
   float beta = 0;
   int N = 3;
   int k = 1;
   int lda = 3;
   float A[] = { 0.627, -0.312, 0.031, 0.308, 0.323, -0.578, 0.797, 0.545, -0.476 };
   float X[] = { -0.542, 0.606, 0.727 };
   int incX = -1;
   float Y[] = { 0.755, 0.268, -0.99 };
   int incY = -1;
   float y_expected[] = { 0.187592, -0.01232, -0.040176 };
   cblas_ssbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "ssbmv(case 1106)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   float alpha = 1;
   float beta = 0;
   int N = 3;
   int k = 1;
   int lda = 3;
   float A[] = { 0.627, -0.312, 0.031, 0.308, 0.323, -0.578, 0.797, 0.545, -0.476 };
   float X[] = { -0.542, 0.606, 0.727 };
   int incX = -1;
   float Y[] = { 0.755, 0.268, -0.99 };
   int incY = -1;
   float y_expected[] = { 0.187592, -0.01232, -0.040176 };
   cblas_ssbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "ssbmv(case 1107)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   float alpha = 1;
   float beta = 0;
   int N = 3;
   int k = 1;
   int lda = 3;
   float A[] = { 0.627, -0.312, 0.031, 0.308, 0.323, -0.578, 0.797, 0.545, -0.476 };
   float X[] = { -0.542, 0.606, 0.727 };
   int incX = -1;
   float Y[] = { 0.755, 0.268, -0.99 };
   int incY = -1;
   float y_expected[] = { -0.236236, -0.215242, 0.266757 };
   cblas_ssbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "ssbmv(case 1108)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   float alpha = 1;
   float beta = 0;
   int N = 3;
   int k = 1;
   int lda = 3;
   float A[] = { 0.627, -0.312, 0.031, 0.308, 0.323, -0.578, 0.797, 0.545, -0.476 };
   float X[] = { -0.542, 0.606, 0.727 };
   int incX = -1;
   float Y[] = { 0.755, 0.268, -0.99 };
   int incY = -1;
   float y_expected[] = { -0.236236, -0.215242, 0.266757 };
   cblas_ssbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "ssbmv(case 1109)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   double alpha = 0;
   double beta = 1;
   int N = 3;
   int k = 1;
   int lda = 3;
   double A[] = { 0.83, -0.568, -0.888, 0.281, -0.779, -0.148, 0.138, 0.053, -0.757 };
   double X[] = { 0.166, 0.808, 0.723 };
   int incX = -1;
   double Y[] = { 0.9, 0.99, -0.578 };
   int incY = -1;
   double y_expected[] = { 0.9, 0.99, -0.578 };
   cblas_dsbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dsbmv(case 1110)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   double alpha = 0;
   double beta = 1;
   int N = 3;
   int k = 1;
   int lda = 3;
   double A[] = { 0.83, -0.568, -0.888, 0.281, -0.779, -0.148, 0.138, 0.053, -0.757 };
   double X[] = { 0.166, 0.808, 0.723 };
   int incX = -1;
   double Y[] = { 0.9, 0.99, -0.578 };
   int incY = -1;
   double y_expected[] = { 0.9, 0.99, -0.578 };
   cblas_dsbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dsbmv(case 1111)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   double alpha = 0;
   double beta = 1;
   int N = 3;
   int k = 1;
   int lda = 3;
   double A[] = { 0.83, -0.568, -0.888, 0.281, -0.779, -0.148, 0.138, 0.053, -0.757 };
   double X[] = { 0.166, 0.808, 0.723 };
   int incX = -1;
   double Y[] = { 0.9, 0.99, -0.578 };
   int incY = -1;
   double y_expected[] = { 0.9, 0.99, -0.578 };
   cblas_dsbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dsbmv(case 1112)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   double alpha = 0;
   double beta = 1;
   int N = 3;
   int k = 1;
   int lda = 3;
   double A[] = { 0.83, -0.568, -0.888, 0.281, -0.779, -0.148, 0.138, 0.053, -0.757 };
   double X[] = { 0.166, 0.808, 0.723 };
   int incX = -1;
   double Y[] = { 0.9, 0.99, -0.578 };
   int incY = -1;
   double y_expected[] = { 0.9, 0.99, -0.578 };
   cblas_dsbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dsbmv(case 1113)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   double alpha = 0;
   double beta = 1;
   int N = 3;
   int k = 1;
   int lda = 3;
   double A[] = { 0.83, -0.568, -0.888, 0.281, -0.779, -0.148, 0.138, 0.053, -0.757 };
   double X[] = { 0.166, 0.808, 0.723 };
   int incX = -1;
   double Y[] = { 0.9, 0.99, -0.578 };
   int incY = -1;
   double y_expected[] = { 0.9, 0.99, -0.578 };
   cblas_dsbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dsbmv(case 1114)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   double alpha = 0;
   double beta = 1;
   int N = 3;
   int k = 1;
   int lda = 3;
   double A[] = { 0.83, -0.568, -0.888, 0.281, -0.779, -0.148, 0.138, 0.053, -0.757 };
   double X[] = { 0.166, 0.808, 0.723 };
   int incX = -1;
   double Y[] = { 0.9, 0.99, -0.578 };
   int incY = -1;
   double y_expected[] = { 0.9, 0.99, -0.578 };
   cblas_dsbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dsbmv(case 1115)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   double alpha = 0;
   double beta = 1;
   int N = 3;
   int k = 1;
   int lda = 3;
   double A[] = { 0.83, -0.568, -0.888, 0.281, -0.779, -0.148, 0.138, 0.053, -0.757 };
   double X[] = { 0.166, 0.808, 0.723 };
   int incX = -1;
   double Y[] = { 0.9, 0.99, -0.578 };
   int incY = -1;
   double y_expected[] = { 0.9, 0.99, -0.578 };
   cblas_dsbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dsbmv(case 1116)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   double alpha = 0;
   double beta = 1;
   int N = 3;
   int k = 1;
   int lda = 3;
   double A[] = { 0.83, -0.568, -0.888, 0.281, -0.779, -0.148, 0.138, 0.053, -0.757 };
   double X[] = { 0.166, 0.808, 0.723 };
   int incX = -1;
   double Y[] = { 0.9, 0.99, -0.578 };
   int incY = -1;
   double y_expected[] = { 0.9, 0.99, -0.578 };
   cblas_dsbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dsbmv(case 1117)");
     }
   };
  };


}
