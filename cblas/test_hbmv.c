#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_hbmv () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int uplo = 121;
   float alpha[2] = {0, 0.1};
   float beta[2] = {1, 0};
   int N = 2;
   int k = 1;
   int lda = 2;
   float A[] = { 0.278, 0.503, 0.255, 0.301, -0.255, -0.211, -0.837, -0.849 };
   float X[] = { 0.735, 0.979, -0.69, 0.714 };
   int incX = 1;
   float Y[] = { -0.387, -0.543, -0.028, 0.547 };
   int incY = -1;
   float y_expected[] = { -0.3300792, -0.4370366, -0.0004735, 0.4891711 };
   cblas_chbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chbmv(case 832) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chbmv(case 832) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   float alpha[2] = {0, 0.1};
   float beta[2] = {1, 0};
   int N = 2;
   int k = 1;
   int lda = 2;
   float A[] = { 0.278, 0.503, 0.255, 0.301, -0.255, -0.211, -0.837, -0.849 };
   float X[] = { 0.735, 0.979, -0.69, 0.714 };
   int incX = 1;
   float Y[] = { -0.387, -0.543, -0.028, 0.547 };
   int incY = -1;
   float y_expected[] = { -0.3300792, -0.4370366, -0.0004735, 0.4891711 };
   cblas_chbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chbmv(case 833) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chbmv(case 833) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   float alpha[2] = {0, 0.1};
   float beta[2] = {1, 0};
   int N = 2;
   int k = 1;
   int lda = 2;
   float A[] = { 0.278, 0.503, 0.255, 0.301, -0.255, -0.211, -0.837, -0.849 };
   float X[] = { 0.735, 0.979, -0.69, 0.714 };
   int incX = 1;
   float Y[] = { -0.387, -0.543, -0.028, 0.547 };
   int incY = -1;
   float y_expected[] = { -0.364734, -0.5586806, -0.0224502, 0.5699626 };
   cblas_chbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chbmv(case 834) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chbmv(case 834) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   float alpha[2] = {0, 0.1};
   float beta[2] = {1, 0};
   int N = 2;
   int k = 1;
   int lda = 2;
   float A[] = { 0.278, 0.503, 0.255, 0.301, -0.255, -0.211, -0.837, -0.849 };
   float X[] = { 0.735, 0.979, -0.69, 0.714 };
   int incX = 1;
   float Y[] = { -0.387, -0.543, -0.028, 0.547 };
   int incY = -1;
   float y_expected[] = { -0.364734, -0.5586806, -0.0224502, 0.5699626 };
   cblas_chbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chbmv(case 835) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chbmv(case 835) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   float alpha[2] = {0, 0.1};
   float beta[2] = {1, 0};
   int N = 2;
   int k = 1;
   int lda = 2;
   float A[] = { 0.278, 0.503, 0.255, 0.301, -0.255, -0.211, -0.837, -0.849 };
   float X[] = { 0.735, 0.979, -0.69, 0.714 };
   int incX = 1;
   float Y[] = { -0.387, -0.543, -0.028, 0.547 };
   int incY = -1;
   float y_expected[] = { -0.3177822, -0.5246464, -0.0493165, 0.5984029 };
   cblas_chbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chbmv(case 836) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chbmv(case 836) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   float alpha[2] = {0, 0.1};
   float beta[2] = {1, 0};
   int N = 2;
   int k = 1;
   int lda = 2;
   float A[] = { 0.278, 0.503, 0.255, 0.301, -0.255, -0.211, -0.837, -0.849 };
   float X[] = { 0.735, 0.979, -0.69, 0.714 };
   int incX = 1;
   float Y[] = { -0.387, -0.543, -0.028, 0.547 };
   int incY = -1;
   float y_expected[] = { -0.3177822, -0.5246464, -0.0493165, 0.5984029 };
   cblas_chbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chbmv(case 837) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chbmv(case 837) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   float alpha[2] = {0, 0.1};
   float beta[2] = {1, 0};
   int N = 2;
   int k = 1;
   int lda = 2;
   float A[] = { 0.278, 0.503, 0.255, 0.301, -0.255, -0.211, -0.837, -0.849 };
   float X[] = { 0.735, 0.979, -0.69, 0.714 };
   int incX = 1;
   float Y[] = { -0.387, -0.543, -0.028, 0.547 };
   int incY = -1;
   float y_expected[] = { -0.415881, -0.5361304, -0.0941922, 0.5713294 };
   cblas_chbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chbmv(case 838) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chbmv(case 838) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   float alpha[2] = {0, 0.1};
   float beta[2] = {1, 0};
   int N = 2;
   int k = 1;
   int lda = 2;
   float A[] = { 0.278, 0.503, 0.255, 0.301, -0.255, -0.211, -0.837, -0.849 };
   float X[] = { 0.735, 0.979, -0.69, 0.714 };
   int incX = 1;
   float Y[] = { -0.387, -0.543, -0.028, 0.547 };
   int incY = -1;
   float y_expected[] = { -0.415881, -0.5361304, -0.0941922, 0.5713294 };
   cblas_chbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chbmv(case 839) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chbmv(case 839) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   double alpha[2] = {1, 0};
   double beta[2] = {0, 0.1};
   int N = 2;
   int k = 1;
   int lda = 2;
   double A[] = { -0.997, 0.203, -0.799, -0.03, -0.101, -0.319, 0.808, -0.728 };
   double X[] = { 0.171, -0.195, -0.019, -0.447 };
   int incX = 1;
   double Y[] = { 0.946, 0.785, -0.921, 0.532 };
   int incY = -1;
   double y_expected[] = { -0.224631, -0.105641, -0.0687, 0.285318 };
   cblas_zhbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhbmv(case 840) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhbmv(case 840) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   double alpha[2] = {1, 0};
   double beta[2] = {0, 0.1};
   int N = 2;
   int k = 1;
   int lda = 2;
   double A[] = { -0.997, 0.203, -0.799, -0.03, -0.101, -0.319, 0.808, -0.728 };
   double X[] = { 0.171, -0.195, -0.019, -0.447 };
   int incX = 1;
   double Y[] = { 0.946, 0.785, -0.921, 0.532 };
   int incY = -1;
   double y_expected[] = { -0.224631, -0.105641, -0.0687, 0.285318 };
   cblas_zhbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhbmv(case 841) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhbmv(case 841) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   double alpha[2] = {1, 0};
   double beta[2] = {0, 0.1};
   int N = 2;
   int k = 1;
   int lda = 2;
   double A[] = { -0.997, 0.203, -0.799, -0.03, -0.101, -0.319, 0.808, -0.728 };
   double X[] = { 0.171, -0.195, -0.019, -0.447 };
   int incX = 1;
   double Y[] = { 0.946, 0.785, -0.921, 0.532 };
   int incY = -1;
   double y_expected[] = { -0.142795, 0.416899, -0.079175, 0.141401 };
   cblas_zhbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhbmv(case 842) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhbmv(case 842) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   double alpha[2] = {1, 0};
   double beta[2] = {0, 0.1};
   int N = 2;
   int k = 1;
   int lda = 2;
   double A[] = { -0.997, 0.203, -0.799, -0.03, -0.101, -0.319, 0.808, -0.728 };
   double X[] = { 0.171, -0.195, -0.019, -0.447 };
   int incX = 1;
   double Y[] = { 0.946, 0.785, -0.921, 0.532 };
   int incY = -1;
   double y_expected[] = { -0.142795, 0.416899, -0.079175, 0.141401 };
   cblas_zhbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhbmv(case 843) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhbmv(case 843) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   double alpha[2] = {1, 0};
   double beta[2] = {0, 0.1};
   int N = 2;
   int k = 1;
   int lda = 2;
   double A[] = { -0.997, 0.203, -0.799, -0.03, -0.101, -0.319, 0.808, -0.728 };
   double X[] = { 0.171, -0.195, -0.019, -0.447 };
   int incX = 1;
   double Y[] = { 0.946, 0.785, -0.921, 0.532 };
   int incY = -1;
   double y_expected[] = { -0.048918, -0.192332, -0.330503, 0.114913 };
   cblas_zhbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhbmv(case 844) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhbmv(case 844) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   double alpha[2] = {1, 0};
   double beta[2] = {0, 0.1};
   int N = 2;
   int k = 1;
   int lda = 2;
   double A[] = { -0.997, 0.203, -0.799, -0.03, -0.101, -0.319, 0.808, -0.728 };
   double X[] = { 0.171, -0.195, -0.019, -0.447 };
   int incX = 1;
   double Y[] = { 0.946, 0.785, -0.921, 0.532 };
   int incY = -1;
   double y_expected[] = { -0.048918, -0.192332, -0.330503, 0.114913 };
   cblas_zhbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhbmv(case 845) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhbmv(case 845) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   double alpha[2] = {1, 0};
   double beta[2] = {0, 0.1};
   int N = 2;
   int k = 1;
   int lda = 2;
   double A[] = { -0.997, 0.203, -0.799, -0.03, -0.101, -0.319, 0.808, -0.728 };
   double X[] = { 0.171, -0.195, -0.019, -0.447 };
   int incX = 1;
   double Y[] = { 0.946, 0.785, -0.921, 0.532 };
   int incY = -1;
   double y_expected[] = { -0.21906, 0.290422, -0.195096, 0.458898 };
   cblas_zhbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhbmv(case 846) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhbmv(case 846) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   double alpha[2] = {1, 0};
   double beta[2] = {0, 0.1};
   int N = 2;
   int k = 1;
   int lda = 2;
   double A[] = { -0.997, 0.203, -0.799, -0.03, -0.101, -0.319, 0.808, -0.728 };
   double X[] = { 0.171, -0.195, -0.019, -0.447 };
   int incX = 1;
   double Y[] = { 0.946, 0.785, -0.921, 0.532 };
   int incY = -1;
   double y_expected[] = { -0.21906, 0.290422, -0.195096, 0.458898 };
   cblas_zhbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhbmv(case 847) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhbmv(case 847) imag");
     };
   };
  };


}
