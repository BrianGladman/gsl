#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

#include "tests.h"

void
test_hemv (void) {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int uplo = 121;
   float alpha[2] = {1, 0};
   float beta[2] = {-0.3, 0.1};
   int N = 1;
   int lda = 1;
   float A[] = { -0.434, 0.837 };
   float X[] = { 0.209, -0.935 };
   int incX = -1;
   float Y[] = { 0.346, -0.412 };
   int incY = -1;
   float y_expected[] = { -0.153306, 0.56399 };
   cblas_chemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chemv(case 1070) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chemv(case 1070) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   float alpha[2] = {1, 0};
   float beta[2] = {-0.3, 0.1};
   int N = 1;
   int lda = 1;
   float A[] = { -0.434, 0.837 };
   float X[] = { 0.209, -0.935 };
   int incX = -1;
   float Y[] = { 0.346, -0.412 };
   int incY = -1;
   float y_expected[] = { -0.153306, 0.56399 };
   cblas_chemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chemv(case 1071) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chemv(case 1071) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   float alpha[2] = {1, 0};
   float beta[2] = {-0.3, 0.1};
   int N = 1;
   int lda = 1;
   float A[] = { -0.434, 0.837 };
   float X[] = { 0.209, -0.935 };
   int incX = -1;
   float Y[] = { 0.346, -0.412 };
   int incY = -1;
   float y_expected[] = { -0.153306, 0.56399 };
   cblas_chemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chemv(case 1072) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chemv(case 1072) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   float alpha[2] = {1, 0};
   float beta[2] = {-0.3, 0.1};
   int N = 1;
   int lda = 1;
   float A[] = { -0.434, 0.837 };
   float X[] = { 0.209, -0.935 };
   int incX = -1;
   float Y[] = { 0.346, -0.412 };
   int incY = -1;
   float y_expected[] = { -0.153306, 0.56399 };
   cblas_chemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chemv(case 1073) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chemv(case 1073) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   float alpha[2] = {1, 0};
   float beta[2] = {-0.3, 0.1};
   int N = 1;
   int lda = 1;
   float A[] = { -0.434, 0.837 };
   float X[] = { 0.209, -0.935 };
   int incX = -1;
   float Y[] = { 0.346, -0.412 };
   int incY = -1;
   float y_expected[] = { -0.153306, 0.56399 };
   cblas_chemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chemv(case 1074) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chemv(case 1074) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   float alpha[2] = {1, 0};
   float beta[2] = {-0.3, 0.1};
   int N = 1;
   int lda = 1;
   float A[] = { -0.434, 0.837 };
   float X[] = { 0.209, -0.935 };
   int incX = -1;
   float Y[] = { 0.346, -0.412 };
   int incY = -1;
   float y_expected[] = { -0.153306, 0.56399 };
   cblas_chemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chemv(case 1075) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chemv(case 1075) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   float alpha[2] = {1, 0};
   float beta[2] = {-0.3, 0.1};
   int N = 1;
   int lda = 1;
   float A[] = { -0.434, 0.837 };
   float X[] = { 0.209, -0.935 };
   int incX = -1;
   float Y[] = { 0.346, -0.412 };
   int incY = -1;
   float y_expected[] = { -0.153306, 0.56399 };
   cblas_chemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chemv(case 1076) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chemv(case 1076) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   float alpha[2] = {1, 0};
   float beta[2] = {-0.3, 0.1};
   int N = 1;
   int lda = 1;
   float A[] = { -0.434, 0.837 };
   float X[] = { 0.209, -0.935 };
   int incX = -1;
   float Y[] = { 0.346, -0.412 };
   int incY = -1;
   float y_expected[] = { -0.153306, 0.56399 };
   cblas_chemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chemv(case 1077) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chemv(case 1077) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   double alpha[2] = {0, 0};
   double beta[2] = {1, 0};
   int N = 1;
   int lda = 1;
   double A[] = { 0.036, -0.966 };
   double X[] = { -0.695, 0.886 };
   int incX = -1;
   double Y[] = { 0.486, 0.629 };
   int incY = -1;
   double y_expected[] = { 0.486, 0.629 };
   cblas_zhemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhemv(case 1078) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhemv(case 1078) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   double alpha[2] = {0, 0};
   double beta[2] = {1, 0};
   int N = 1;
   int lda = 1;
   double A[] = { 0.036, -0.966 };
   double X[] = { -0.695, 0.886 };
   int incX = -1;
   double Y[] = { 0.486, 0.629 };
   int incY = -1;
   double y_expected[] = { 0.486, 0.629 };
   cblas_zhemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhemv(case 1079) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhemv(case 1079) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   double alpha[2] = {0, 0};
   double beta[2] = {1, 0};
   int N = 1;
   int lda = 1;
   double A[] = { 0.036, -0.966 };
   double X[] = { -0.695, 0.886 };
   int incX = -1;
   double Y[] = { 0.486, 0.629 };
   int incY = -1;
   double y_expected[] = { 0.486, 0.629 };
   cblas_zhemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhemv(case 1080) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhemv(case 1080) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   double alpha[2] = {0, 0};
   double beta[2] = {1, 0};
   int N = 1;
   int lda = 1;
   double A[] = { 0.036, -0.966 };
   double X[] = { -0.695, 0.886 };
   int incX = -1;
   double Y[] = { 0.486, 0.629 };
   int incY = -1;
   double y_expected[] = { 0.486, 0.629 };
   cblas_zhemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhemv(case 1081) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhemv(case 1081) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   double alpha[2] = {0, 0};
   double beta[2] = {1, 0};
   int N = 1;
   int lda = 1;
   double A[] = { 0.036, -0.966 };
   double X[] = { -0.695, 0.886 };
   int incX = -1;
   double Y[] = { 0.486, 0.629 };
   int incY = -1;
   double y_expected[] = { 0.486, 0.629 };
   cblas_zhemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhemv(case 1082) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhemv(case 1082) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   double alpha[2] = {0, 0};
   double beta[2] = {1, 0};
   int N = 1;
   int lda = 1;
   double A[] = { 0.036, -0.966 };
   double X[] = { -0.695, 0.886 };
   int incX = -1;
   double Y[] = { 0.486, 0.629 };
   int incY = -1;
   double y_expected[] = { 0.486, 0.629 };
   cblas_zhemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhemv(case 1083) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhemv(case 1083) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   double alpha[2] = {0, 0};
   double beta[2] = {1, 0};
   int N = 1;
   int lda = 1;
   double A[] = { 0.036, -0.966 };
   double X[] = { -0.695, 0.886 };
   int incX = -1;
   double Y[] = { 0.486, 0.629 };
   int incY = -1;
   double y_expected[] = { 0.486, 0.629 };
   cblas_zhemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhemv(case 1084) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhemv(case 1084) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   double alpha[2] = {0, 0};
   double beta[2] = {1, 0};
   int N = 1;
   int lda = 1;
   double A[] = { 0.036, -0.966 };
   double X[] = { -0.695, 0.886 };
   int incX = -1;
   double Y[] = { 0.486, 0.629 };
   int incY = -1;
   double y_expected[] = { 0.486, 0.629 };
   cblas_zhemv(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhemv(case 1085) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhemv(case 1085) imag");
     };
   };
  };


}
