#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_gemv () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int trans = 111;
   int M = 1;
   int N = 1;
   int lda = 1;
   float alpha = 0.1;
   float beta = 0;
   float A[] = { -0.603 };
   float X[] = { 0.339 };
   int incX = 1;
   float Y[] = { 0.239 };
   int incY = -1;
   float y_expected[] = { -0.0204417 };
   cblas_sgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sgemv(case 520)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int M = 1;
   int N = 1;
   int lda = 1;
   float alpha = 0.1;
   float beta = 0;
   float A[] = { -0.603 };
   float X[] = { 0.339 };
   int incX = 1;
   float Y[] = { 0.239 };
   int incY = -1;
   float y_expected[] = { -0.0204417 };
   cblas_sgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sgemv(case 521)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int M = 1;
   int N = 1;
   int lda = 1;
   float alpha = -1;
   float beta = -0.3;
   float A[] = { -0.603 };
   float X[] = { 0.339 };
   int incX = 1;
   float Y[] = { 0.239 };
   int incY = -1;
   float y_expected[] = { 0.132717 };
   cblas_sgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sgemv(case 522)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int M = 1;
   int N = 1;
   int lda = 1;
   float alpha = -1;
   float beta = -0.3;
   float A[] = { -0.603 };
   float X[] = { 0.339 };
   int incX = 1;
   float Y[] = { 0.239 };
   int incY = -1;
   float y_expected[] = { 0.132717 };
   cblas_sgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sgemv(case 523)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int M = 1;
   int N = 1;
   int lda = 1;
   double alpha = -1;
   double beta = 0;
   double A[] = { 0.503 };
   double X[] = { 0.629 };
   int incX = 1;
   double Y[] = { -0.419 };
   int incY = -1;
   double y_expected[] = { -0.316387 };
   cblas_dgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dgemv(case 524)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int M = 1;
   int N = 1;
   int lda = 1;
   double alpha = -1;
   double beta = 0;
   double A[] = { 0.503 };
   double X[] = { 0.629 };
   int incX = 1;
   double Y[] = { -0.419 };
   int incY = -1;
   double y_expected[] = { -0.316387 };
   cblas_dgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dgemv(case 525)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int M = 1;
   int N = 1;
   int lda = 1;
   double alpha = -0.3;
   double beta = 0.1;
   double A[] = { 0.503 };
   double X[] = { 0.629 };
   int incX = 1;
   double Y[] = { -0.419 };
   int incY = -1;
   double y_expected[] = { -0.1368161 };
   cblas_dgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dgemv(case 526)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int M = 1;
   int N = 1;
   int lda = 1;
   double alpha = -0.3;
   double beta = 0.1;
   double A[] = { 0.503 };
   double X[] = { 0.629 };
   int incX = 1;
   double Y[] = { -0.419 };
   int incY = -1;
   double y_expected[] = { -0.1368161 };
   cblas_dgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dgemv(case 527)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int M = 1;
   int N = 1;
   int lda = 1;
   float alpha[2] = {1, 0};
   float beta[2] = {0, 0.1};
   float A[] = { 0.331, 0.622 };
   float X[] = { 0.521, 0.556 };
   int incX = 1;
   float Y[] = { -0.811, -0.147 };
   int incY = -1;
   float y_expected[] = { -0.158681, 0.426998 };
   cblas_cgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "cgemv(case 528) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "cgemv(case 528) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int M = 1;
   int N = 1;
   int lda = 1;
   float alpha[2] = {1, 0};
   float beta[2] = {0, 0.1};
   float A[] = { 0.331, 0.622 };
   float X[] = { 0.521, 0.556 };
   int incX = 1;
   float Y[] = { -0.811, -0.147 };
   int incY = -1;
   float y_expected[] = { -0.158681, 0.426998 };
   cblas_cgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "cgemv(case 529) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "cgemv(case 529) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int M = 1;
   int N = 1;
   int lda = 1;
   float alpha[2] = {0, 0.1};
   float beta[2] = {1, 0};
   float A[] = { 0.331, 0.622 };
   float X[] = { 0.521, 0.556 };
   int incX = 1;
   float Y[] = { -0.811, -0.147 };
   int incY = -1;
   float y_expected[] = { -0.8618098, -0.1643381 };
   cblas_cgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "cgemv(case 530) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "cgemv(case 530) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int M = 1;
   int N = 1;
   int lda = 1;
   float alpha[2] = {0, 0.1};
   float beta[2] = {1, 0};
   float A[] = { 0.331, 0.622 };
   float X[] = { 0.521, 0.556 };
   int incX = 1;
   float Y[] = { -0.811, -0.147 };
   int incY = -1;
   float y_expected[] = { -0.8618098, -0.1643381 };
   cblas_cgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "cgemv(case 531) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "cgemv(case 531) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int M = 1;
   int N = 1;
   int lda = 1;
   float alpha[2] = {0, 0};
   float beta[2] = {0, 0.1};
   float A[] = { 0.331, 0.622 };
   float X[] = { 0.521, 0.556 };
   int incX = 1;
   float Y[] = { -0.811, -0.147 };
   int incY = -1;
   float y_expected[] = { 0.0147, -0.0811 };
   cblas_cgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "cgemv(case 532) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "cgemv(case 532) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int M = 1;
   int N = 1;
   int lda = 1;
   float alpha[2] = {0, 0};
   float beta[2] = {0, 0.1};
   float A[] = { 0.331, 0.622 };
   float X[] = { 0.521, 0.556 };
   int incX = 1;
   float Y[] = { -0.811, -0.147 };
   int incY = -1;
   float y_expected[] = { 0.0147, -0.0811 };
   cblas_cgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "cgemv(case 533) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "cgemv(case 533) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int M = 1;
   int N = 1;
   int lda = 1;
   double alpha[2] = {-1, 0};
   double beta[2] = {0, 0};
   double A[] = { -0.176, -0.165 };
   double X[] = { -0.201, 0.087 };
   int incX = 1;
   double Y[] = { -0.464, 0.7 };
   int incY = -1;
   double y_expected[] = { -0.049731, -0.017853 };
   cblas_zgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zgemv(case 534) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zgemv(case 534) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int M = 1;
   int N = 1;
   int lda = 1;
   double alpha[2] = {-1, 0};
   double beta[2] = {0, 0};
   double A[] = { -0.176, -0.165 };
   double X[] = { -0.201, 0.087 };
   int incX = 1;
   double Y[] = { -0.464, 0.7 };
   int incY = -1;
   double y_expected[] = { -0.049731, -0.017853 };
   cblas_zgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zgemv(case 535) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zgemv(case 535) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int M = 1;
   int N = 1;
   int lda = 1;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {-1, 0};
   double A[] = { -0.176, -0.165 };
   double X[] = { -0.201, 0.087 };
   int incX = 1;
   double Y[] = { -0.464, 0.7 };
   int incY = -1;
   double y_expected[] = { 0.4472954, -0.7003828 };
   cblas_zgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zgemv(case 536) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zgemv(case 536) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int M = 1;
   int N = 1;
   int lda = 1;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {-1, 0};
   double A[] = { -0.176, -0.165 };
   double X[] = { -0.201, 0.087 };
   int incX = 1;
   double Y[] = { -0.464, 0.7 };
   int incY = -1;
   double y_expected[] = { 0.4472954, -0.7003828 };
   cblas_zgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zgemv(case 537) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zgemv(case 537) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int M = 1;
   int N = 1;
   int lda = 1;
   double alpha[2] = {-1, 0};
   double beta[2] = {-0.3, 0.1};
   double A[] = { -0.176, -0.165 };
   double X[] = { -0.201, 0.087 };
   int incX = 1;
   double Y[] = { -0.464, 0.7 };
   int incY = -1;
   double y_expected[] = { 0.048179, -0.207923 };
   cblas_zgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zgemv(case 538) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zgemv(case 538) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int M = 1;
   int N = 1;
   int lda = 1;
   double alpha[2] = {-1, 0};
   double beta[2] = {-0.3, 0.1};
   double A[] = { -0.176, -0.165 };
   double X[] = { -0.201, 0.087 };
   int incX = 1;
   double Y[] = { -0.464, 0.7 };
   int incY = -1;
   double y_expected[] = { 0.048179, -0.207923 };
   cblas_zgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zgemv(case 539) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zgemv(case 539) imag");
     };
   };
  };


}
