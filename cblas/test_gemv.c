#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

#include "tests.h"

void
test_gemv (void) {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int trans = 111;
   int M = 1;
   int N = 1;
   int lda = 1;
   float alpha = 1;
   float beta = -0.3;
   float A[] = { -0.805 };
   float X[] = { -0.965 };
   int incX = -1;
   float Y[] = { 0.537 };
   int incY = -1;
   float y_expected[] = { 0.615725 };
   cblas_sgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sgemv(case 774)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int M = 1;
   int N = 1;
   int lda = 1;
   float alpha = 1;
   float beta = -0.3;
   float A[] = { -0.805 };
   float X[] = { -0.965 };
   int incX = -1;
   float Y[] = { 0.537 };
   int incY = -1;
   float y_expected[] = { 0.615725 };
   cblas_sgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sgemv(case 775)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int M = 1;
   int N = 1;
   int lda = 1;
   float alpha = 1;
   float beta = 0;
   float A[] = { -0.805 };
   float X[] = { -0.965 };
   int incX = -1;
   float Y[] = { 0.537 };
   int incY = -1;
   float y_expected[] = { 0.776825 };
   cblas_sgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sgemv(case 776)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int M = 1;
   int N = 1;
   int lda = 1;
   float alpha = 1;
   float beta = 0;
   float A[] = { -0.805 };
   float X[] = { -0.965 };
   int incX = -1;
   float Y[] = { 0.537 };
   int incY = -1;
   float y_expected[] = { 0.776825 };
   cblas_sgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sgemv(case 777)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int M = 1;
   int N = 1;
   int lda = 1;
   double alpha = -0.3;
   double beta = -1;
   double A[] = { -0.047 };
   double X[] = { 0.672 };
   int incX = -1;
   double Y[] = { 0.554 };
   int incY = -1;
   double y_expected[] = { -0.5445248 };
   cblas_dgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dgemv(case 778)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int M = 1;
   int N = 1;
   int lda = 1;
   double alpha = -0.3;
   double beta = -1;
   double A[] = { -0.047 };
   double X[] = { 0.672 };
   int incX = -1;
   double Y[] = { 0.554 };
   int incY = -1;
   double y_expected[] = { -0.5445248 };
   cblas_dgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dgemv(case 779)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int M = 1;
   int N = 1;
   int lda = 1;
   double alpha = -1;
   double beta = 1;
   double A[] = { -0.047 };
   double X[] = { 0.672 };
   int incX = -1;
   double Y[] = { 0.554 };
   int incY = -1;
   double y_expected[] = { 0.585584 };
   cblas_dgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dgemv(case 780)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int M = 1;
   int N = 1;
   int lda = 1;
   double alpha = -1;
   double beta = 1;
   double A[] = { -0.047 };
   double X[] = { 0.672 };
   int incX = -1;
   double Y[] = { 0.554 };
   int incY = -1;
   double y_expected[] = { 0.585584 };
   cblas_dgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dgemv(case 781)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int M = 1;
   int N = 1;
   int lda = 1;
   float alpha[2] = {0, 0.1};
   float beta[2] = {0, 1};
   float A[] = { 0.629, 0.801 };
   float X[] = { 0.778, -0.073 };
   int incX = -1;
   float Y[] = { -0.976, -0.682 };
   int incY = -1;
   float y_expected[] = { 0.624274, -0.921216 };
   cblas_cgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "cgemv(case 782) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "cgemv(case 782) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int M = 1;
   int N = 1;
   int lda = 1;
   float alpha[2] = {0, 0.1};
   float beta[2] = {0, 1};
   float A[] = { 0.629, 0.801 };
   float X[] = { 0.778, -0.073 };
   int incX = -1;
   float Y[] = { -0.976, -0.682 };
   int incY = -1;
   float y_expected[] = { 0.624274, -0.921216 };
   cblas_cgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "cgemv(case 783) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "cgemv(case 783) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int M = 1;
   int N = 1;
   int lda = 1;
   float alpha[2] = {0, 1};
   float beta[2] = {-0.3, 0.1};
   float A[] = { 0.629, 0.801 };
   float X[] = { 0.778, -0.073 };
   int incX = -1;
   float Y[] = { -0.976, -0.682 };
   int incY = -1;
   float y_expected[] = { -0.216261, 0.654835 };
   cblas_cgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "cgemv(case 784) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "cgemv(case 784) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int M = 1;
   int N = 1;
   int lda = 1;
   float alpha[2] = {0, 1};
   float beta[2] = {-0.3, 0.1};
   float A[] = { 0.629, 0.801 };
   float X[] = { 0.778, -0.073 };
   int incX = -1;
   float Y[] = { -0.976, -0.682 };
   int incY = -1;
   float y_expected[] = { -0.216261, 0.654835 };
   cblas_cgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "cgemv(case 785) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "cgemv(case 785) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int M = 1;
   int N = 1;
   int lda = 1;
   float alpha[2] = {0, 0.1};
   float beta[2] = {-0.3, 0.1};
   float A[] = { 0.629, 0.801 };
   float X[] = { 0.778, -0.073 };
   int incX = -1;
   float Y[] = { -0.976, -0.682 };
   int incY = -1;
   float y_expected[] = { 0.427909, 0.150089 };
   cblas_cgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "cgemv(case 786) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "cgemv(case 786) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int M = 1;
   int N = 1;
   int lda = 1;
   float alpha[2] = {0, 0.1};
   float beta[2] = {-0.3, 0.1};
   float A[] = { 0.629, 0.801 };
   float X[] = { 0.778, -0.073 };
   int incX = -1;
   float Y[] = { -0.976, -0.682 };
   int incY = -1;
   float y_expected[] = { 0.427909, 0.150089 };
   cblas_cgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "cgemv(case 787) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "cgemv(case 787) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int M = 1;
   int N = 1;
   int lda = 1;
   double alpha[2] = {0, 0.1};
   double beta[2] = {1, 0};
   double A[] = { 0.932, -0.724 };
   double X[] = { 0.334, -0.317 };
   int incX = -1;
   double Y[] = { 0.348, 0.07 };
   int incY = -1;
   double y_expected[] = { 0.401726, 0.078178 };
   cblas_zgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zgemv(case 788) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zgemv(case 788) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int M = 1;
   int N = 1;
   int lda = 1;
   double alpha[2] = {0, 0.1};
   double beta[2] = {1, 0};
   double A[] = { 0.932, -0.724 };
   double X[] = { 0.334, -0.317 };
   int incX = -1;
   double Y[] = { 0.348, 0.07 };
   int incY = -1;
   double y_expected[] = { 0.401726, 0.078178 };
   cblas_zgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zgemv(case 789) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zgemv(case 789) imag");
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
   double beta[2] = {0, 1};
   double A[] = { 0.932, -0.724 };
   double X[] = { 0.334, -0.317 };
   int incX = -1;
   double Y[] = { 0.348, 0.07 };
   int incY = -1;
   double y_expected[] = { -0.040808, 0.517356 };
   cblas_zgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zgemv(case 790) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zgemv(case 790) imag");
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
   double beta[2] = {0, 1};
   double A[] = { 0.932, -0.724 };
   double X[] = { 0.334, -0.317 };
   int incX = -1;
   double Y[] = { 0.348, 0.07 };
   int incY = -1;
   double y_expected[] = { -0.040808, 0.517356 };
   cblas_zgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zgemv(case 791) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zgemv(case 791) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int M = 1;
   int N = 1;
   int lda = 1;
   double alpha[2] = {1, 0};
   double beta[2] = {0, 0};
   double A[] = { 0.932, -0.724 };
   double X[] = { 0.334, -0.317 };
   int incX = -1;
   double Y[] = { 0.348, 0.07 };
   int incY = -1;
   double y_expected[] = { 0.540796, -0.053628 };
   cblas_zgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zgemv(case 792) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zgemv(case 792) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int M = 1;
   int N = 1;
   int lda = 1;
   double alpha[2] = {1, 0};
   double beta[2] = {0, 0};
   double A[] = { 0.932, -0.724 };
   double X[] = { 0.334, -0.317 };
   int incX = -1;
   double Y[] = { 0.348, 0.07 };
   int incY = -1;
   double y_expected[] = { 0.540796, -0.053628 };
   cblas_zgemv(order, trans, M, N, alpha, A, lda, X,                                  incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zgemv(case 793) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zgemv(case 793) imag");
     };
   };
  };


}
