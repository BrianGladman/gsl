#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_hpmv () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int uplo = 121;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {0, 1};
   int N = 1;
   float A[] = { 0.838, 0.503 };
   float X[] = { 0.943, -0.41 };
   int incX = 1;
   float Y[] = { 0.937, -0.344 };
   int incY = -1;
   float y_expected[] = { 0.1412878, 1.1190974 };
   cblas_chpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chpmv(case 864) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chpmv(case 864) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {0, 1};
   int N = 1;
   float A[] = { 0.838, 0.503 };
   float X[] = { 0.943, -0.41 };
   int incX = 1;
   float Y[] = { 0.937, -0.344 };
   int incY = -1;
   float y_expected[] = { 0.1412878, 1.1190974 };
   cblas_chpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chpmv(case 865) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chpmv(case 865) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {0, 1};
   int N = 1;
   float A[] = { 0.838, 0.503 };
   float X[] = { 0.943, -0.41 };
   int incX = 1;
   float Y[] = { 0.937, -0.344 };
   int incY = -1;
   float y_expected[] = { 0.1412878, 1.1190974 };
   cblas_chpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chpmv(case 866) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chpmv(case 866) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {0, 1};
   int N = 1;
   float A[] = { 0.838, 0.503 };
   float X[] = { 0.943, -0.41 };
   int incX = 1;
   float Y[] = { 0.937, -0.344 };
   int incY = -1;
   float y_expected[] = { 0.1412878, 1.1190974 };
   cblas_chpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chpmv(case 867) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chpmv(case 867) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {0, 1};
   int N = 1;
   float A[] = { 0.838, 0.503 };
   float X[] = { 0.943, -0.41 };
   int incX = 1;
   float Y[] = { 0.937, -0.344 };
   int incY = -1;
   float y_expected[] = { 0.1412878, 1.1190974 };
   cblas_chpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chpmv(case 868) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chpmv(case 868) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {0, 1};
   int N = 1;
   float A[] = { 0.838, 0.503 };
   float X[] = { 0.943, -0.41 };
   int incX = 1;
   float Y[] = { 0.937, -0.344 };
   int incY = -1;
   float y_expected[] = { 0.1412878, 1.1190974 };
   cblas_chpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chpmv(case 869) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chpmv(case 869) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {0, 1};
   int N = 1;
   float A[] = { 0.838, 0.503 };
   float X[] = { 0.943, -0.41 };
   int incX = 1;
   float Y[] = { 0.937, -0.344 };
   int incY = -1;
   float y_expected[] = { 0.1412878, 1.1190974 };
   cblas_chpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chpmv(case 870) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chpmv(case 870) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {0, 1};
   int N = 1;
   float A[] = { 0.838, 0.503 };
   float X[] = { 0.943, -0.41 };
   int incX = 1;
   float Y[] = { 0.937, -0.344 };
   int incY = -1;
   float y_expected[] = { 0.1412878, 1.1190974 };
   cblas_chpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chpmv(case 871) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chpmv(case 871) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   double alpha[2] = {1, 0};
   double beta[2] = {0, 0};
   int N = 1;
   double A[] = { 0.103, -0.985 };
   double X[] = { -0.33, -0.51 };
   int incX = 1;
   double Y[] = { -0.399, -0.017 };
   int incY = -1;
   double y_expected[] = { -0.03399, -0.05253 };
   cblas_zhpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhpmv(case 872) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhpmv(case 872) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   double alpha[2] = {1, 0};
   double beta[2] = {0, 0};
   int N = 1;
   double A[] = { 0.103, -0.985 };
   double X[] = { -0.33, -0.51 };
   int incX = 1;
   double Y[] = { -0.399, -0.017 };
   int incY = -1;
   double y_expected[] = { -0.03399, -0.05253 };
   cblas_zhpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhpmv(case 873) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhpmv(case 873) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   double alpha[2] = {1, 0};
   double beta[2] = {0, 0};
   int N = 1;
   double A[] = { 0.103, -0.985 };
   double X[] = { -0.33, -0.51 };
   int incX = 1;
   double Y[] = { -0.399, -0.017 };
   int incY = -1;
   double y_expected[] = { -0.03399, -0.05253 };
   cblas_zhpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhpmv(case 874) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhpmv(case 874) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   double alpha[2] = {1, 0};
   double beta[2] = {0, 0};
   int N = 1;
   double A[] = { 0.103, -0.985 };
   double X[] = { -0.33, -0.51 };
   int incX = 1;
   double Y[] = { -0.399, -0.017 };
   int incY = -1;
   double y_expected[] = { -0.03399, -0.05253 };
   cblas_zhpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhpmv(case 875) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhpmv(case 875) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   double alpha[2] = {1, 0};
   double beta[2] = {0, 0};
   int N = 1;
   double A[] = { 0.103, -0.985 };
   double X[] = { -0.33, -0.51 };
   int incX = 1;
   double Y[] = { -0.399, -0.017 };
   int incY = -1;
   double y_expected[] = { -0.03399, -0.05253 };
   cblas_zhpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhpmv(case 876) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhpmv(case 876) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   double alpha[2] = {1, 0};
   double beta[2] = {0, 0};
   int N = 1;
   double A[] = { 0.103, -0.985 };
   double X[] = { -0.33, -0.51 };
   int incX = 1;
   double Y[] = { -0.399, -0.017 };
   int incY = -1;
   double y_expected[] = { -0.03399, -0.05253 };
   cblas_zhpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhpmv(case 877) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhpmv(case 877) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   double alpha[2] = {1, 0};
   double beta[2] = {0, 0};
   int N = 1;
   double A[] = { 0.103, -0.985 };
   double X[] = { -0.33, -0.51 };
   int incX = 1;
   double Y[] = { -0.399, -0.017 };
   int incY = -1;
   double y_expected[] = { -0.03399, -0.05253 };
   cblas_zhpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhpmv(case 878) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhpmv(case 878) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   double alpha[2] = {1, 0};
   double beta[2] = {0, 0};
   int N = 1;
   double A[] = { 0.103, -0.985 };
   double X[] = { -0.33, -0.51 };
   int incX = 1;
   double Y[] = { -0.399, -0.017 };
   int incY = -1;
   double y_expected[] = { -0.03399, -0.05253 };
   cblas_zhpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhpmv(case 879) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhpmv(case 879) imag");
     };
   };
  };


}
