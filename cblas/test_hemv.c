#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_hemv () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int uplo = 121;
   float alpha[2] = {0, 1};
   float beta[2] = {-0.3, 0.1};
   int N = 1;
   int lda = 1;
   float A[] = { -0.975, 0.108 };
   float X[] = { -0.459, 0.905 };
   int incX = 1;
   float Y[] = { -0.787, 0.329 };
   int incY = -1;
   float y_expected[] = { 1.085575, 0.270125 };
   cblas_chemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chemv(case 816) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chemv(case 816) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   float alpha[2] = {0, 1};
   float beta[2] = {-0.3, 0.1};
   int N = 1;
   int lda = 1;
   float A[] = { -0.975, 0.108 };
   float X[] = { -0.459, 0.905 };
   int incX = 1;
   float Y[] = { -0.787, 0.329 };
   int incY = -1;
   float y_expected[] = { 1.085575, 0.270125 };
   cblas_chemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chemv(case 817) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chemv(case 817) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   float alpha[2] = {0, 1};
   float beta[2] = {-0.3, 0.1};
   int N = 1;
   int lda = 1;
   float A[] = { -0.975, 0.108 };
   float X[] = { -0.459, 0.905 };
   int incX = 1;
   float Y[] = { -0.787, 0.329 };
   int incY = -1;
   float y_expected[] = { 1.085575, 0.270125 };
   cblas_chemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chemv(case 818) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chemv(case 818) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   float alpha[2] = {0, 1};
   float beta[2] = {-0.3, 0.1};
   int N = 1;
   int lda = 1;
   float A[] = { -0.975, 0.108 };
   float X[] = { -0.459, 0.905 };
   int incX = 1;
   float Y[] = { -0.787, 0.329 };
   int incY = -1;
   float y_expected[] = { 1.085575, 0.270125 };
   cblas_chemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chemv(case 819) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chemv(case 819) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   float alpha[2] = {0, 1};
   float beta[2] = {-0.3, 0.1};
   int N = 1;
   int lda = 1;
   float A[] = { -0.975, 0.108 };
   float X[] = { -0.459, 0.905 };
   int incX = 1;
   float Y[] = { -0.787, 0.329 };
   int incY = -1;
   float y_expected[] = { 1.085575, 0.270125 };
   cblas_chemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chemv(case 820) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chemv(case 820) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   float alpha[2] = {0, 1};
   float beta[2] = {-0.3, 0.1};
   int N = 1;
   int lda = 1;
   float A[] = { -0.975, 0.108 };
   float X[] = { -0.459, 0.905 };
   int incX = 1;
   float Y[] = { -0.787, 0.329 };
   int incY = -1;
   float y_expected[] = { 1.085575, 0.270125 };
   cblas_chemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chemv(case 821) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chemv(case 821) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   float alpha[2] = {0, 1};
   float beta[2] = {-0.3, 0.1};
   int N = 1;
   int lda = 1;
   float A[] = { -0.975, 0.108 };
   float X[] = { -0.459, 0.905 };
   int incX = 1;
   float Y[] = { -0.787, 0.329 };
   int incY = -1;
   float y_expected[] = { 1.085575, 0.270125 };
   cblas_chemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chemv(case 822) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chemv(case 822) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   float alpha[2] = {0, 1};
   float beta[2] = {-0.3, 0.1};
   int N = 1;
   int lda = 1;
   float A[] = { -0.975, 0.108 };
   float X[] = { -0.459, 0.905 };
   int incX = 1;
   float Y[] = { -0.787, 0.329 };
   int incY = -1;
   float y_expected[] = { 1.085575, 0.270125 };
   cblas_chemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chemv(case 823) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chemv(case 823) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   double alpha[2] = {0, 0};
   double beta[2] = {-0.3, 0.1};
   int N = 1;
   int lda = 1;
   double A[] = { -0.967, -0.253 };
   double X[] = { 0.446, 0.982 };
   int incX = 1;
   double Y[] = { 0.637, 0.517 };
   int incY = -1;
   double y_expected[] = { -0.2428, -0.0914 };
   cblas_zhemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhemv(case 824) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhemv(case 824) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   double alpha[2] = {0, 0};
   double beta[2] = {-0.3, 0.1};
   int N = 1;
   int lda = 1;
   double A[] = { -0.967, -0.253 };
   double X[] = { 0.446, 0.982 };
   int incX = 1;
   double Y[] = { 0.637, 0.517 };
   int incY = -1;
   double y_expected[] = { -0.2428, -0.0914 };
   cblas_zhemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhemv(case 825) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhemv(case 825) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   double alpha[2] = {0, 0};
   double beta[2] = {-0.3, 0.1};
   int N = 1;
   int lda = 1;
   double A[] = { -0.967, -0.253 };
   double X[] = { 0.446, 0.982 };
   int incX = 1;
   double Y[] = { 0.637, 0.517 };
   int incY = -1;
   double y_expected[] = { -0.2428, -0.0914 };
   cblas_zhemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhemv(case 826) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhemv(case 826) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   double alpha[2] = {0, 0};
   double beta[2] = {-0.3, 0.1};
   int N = 1;
   int lda = 1;
   double A[] = { -0.967, -0.253 };
   double X[] = { 0.446, 0.982 };
   int incX = 1;
   double Y[] = { 0.637, 0.517 };
   int incY = -1;
   double y_expected[] = { -0.2428, -0.0914 };
   cblas_zhemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhemv(case 827) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhemv(case 827) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   double alpha[2] = {0, 0};
   double beta[2] = {-0.3, 0.1};
   int N = 1;
   int lda = 1;
   double A[] = { -0.967, -0.253 };
   double X[] = { 0.446, 0.982 };
   int incX = 1;
   double Y[] = { 0.637, 0.517 };
   int incY = -1;
   double y_expected[] = { -0.2428, -0.0914 };
   cblas_zhemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhemv(case 828) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhemv(case 828) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   double alpha[2] = {0, 0};
   double beta[2] = {-0.3, 0.1};
   int N = 1;
   int lda = 1;
   double A[] = { -0.967, -0.253 };
   double X[] = { 0.446, 0.982 };
   int incX = 1;
   double Y[] = { 0.637, 0.517 };
   int incY = -1;
   double y_expected[] = { -0.2428, -0.0914 };
   cblas_zhemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhemv(case 829) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhemv(case 829) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   double alpha[2] = {0, 0};
   double beta[2] = {-0.3, 0.1};
   int N = 1;
   int lda = 1;
   double A[] = { -0.967, -0.253 };
   double X[] = { 0.446, 0.982 };
   int incX = 1;
   double Y[] = { 0.637, 0.517 };
   int incY = -1;
   double y_expected[] = { -0.2428, -0.0914 };
   cblas_zhemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhemv(case 830) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhemv(case 830) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   double alpha[2] = {0, 0};
   double beta[2] = {-0.3, 0.1};
   int N = 1;
   int lda = 1;
   double A[] = { -0.967, -0.253 };
   double X[] = { 0.446, 0.982 };
   int incX = 1;
   double Y[] = { 0.637, 0.517 };
   int incY = -1;
   double y_expected[] = { -0.2428, -0.0914 };
   cblas_zhemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhemv(case 831) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhemv(case 831) imag");
     };
   };
  };


}
