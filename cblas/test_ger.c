#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_ger () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int M = 1;
   int N = 1;
   int lda = 1;
   float alpha = 1;
   float A[] = { 0.108 };
   float X[] = { -0.496 };
   int incX = 1;
   float Y[] = { 0.927 };
   int incY = -1;
   float A_expected[] = { -0.351792 };
   cblas_sger(order, M, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[i], A_expected[i], flteps, "sger(case 1136)");
     }
   };
  };


  {
   int order = 102;
   int M = 1;
   int N = 1;
   int lda = 1;
   float alpha = 1;
   float A[] = { 0.108 };
   float X[] = { -0.496 };
   int incX = 1;
   float Y[] = { 0.927 };
   int incY = -1;
   float A_expected[] = { -0.351792 };
   cblas_sger(order, M, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[i], A_expected[i], flteps, "sger(case 1137)");
     }
   };
  };


  {
   int order = 101;
   int M = 1;
   int N = 1;
   int lda = 1;
   double alpha = -0.3;
   double A[] = { 0.729 };
   double X[] = { -0.92 };
   int incX = 1;
   double Y[] = { -0.469 };
   int incY = -1;
   double A_expected[] = { 0.599556 };
   cblas_dger(order, M, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[i], A_expected[i], dbleps, "dger(case 1138)");
     }
   };
  };


  {
   int order = 102;
   int M = 1;
   int N = 1;
   int lda = 1;
   double alpha = -0.3;
   double A[] = { 0.729 };
   double X[] = { -0.92 };
   int incX = 1;
   double Y[] = { -0.469 };
   int incY = -1;
   double A_expected[] = { 0.599556 };
   cblas_dger(order, M, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[i], A_expected[i], dbleps, "dger(case 1139)");
     }
   };
  };


  {
   int order = 101;
   int M = 1;
   int N = 1;
   int lda = 1;
   float alpha[2] = {1, 0};
   float A[] = { -0.877, -0.308 };
   float X[] = { -0.806, 0.771 };
   int incX = 1;
   float Y[] = { 0.637, -0.177 };
   int incY = -1;
   float A_expected[] = { -1.253955, 0.325789 };
   cblas_cgeru(order, M, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[2*i], A_expected[2*i], flteps, "cgeru(case 1140) real");
       gsl_test_rel(A[2*i+1], A_expected[2*i+1], flteps, "cgeru(case 1140) imag");
     };
   };
  };


  {
   int order = 101;
   int M = 1;
   int N = 1;
   int lda = 1;
   float alpha[2] = {1, 0};
   float A[] = { -0.877, -0.308 };
   float X[] = { -0.806, 0.771 };
   int incX = 1;
   float Y[] = { 0.637, -0.177 };
   int incY = -1;
   float A_expected[] = { -1.526889, 0.040465 };
   cblas_cgerc(order, M, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[2*i], A_expected[2*i], flteps, "cgerc(case 1141) real");
       gsl_test_rel(A[2*i+1], A_expected[2*i+1], flteps, "cgerc(case 1141) imag");
     };
   };
  };


  {
   int order = 102;
   int M = 1;
   int N = 1;
   int lda = 1;
   float alpha[2] = {1, 0};
   float A[] = { -0.877, -0.308 };
   float X[] = { -0.806, 0.771 };
   int incX = 1;
   float Y[] = { 0.637, -0.177 };
   int incY = -1;
   float A_expected[] = { -1.253955, 0.325789 };
   cblas_cgeru(order, M, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[2*i], A_expected[2*i], flteps, "cgeru(case 1142) real");
       gsl_test_rel(A[2*i+1], A_expected[2*i+1], flteps, "cgeru(case 1142) imag");
     };
   };
  };


  {
   int order = 102;
   int M = 1;
   int N = 1;
   int lda = 1;
   float alpha[2] = {1, 0};
   float A[] = { -0.877, -0.308 };
   float X[] = { -0.806, 0.771 };
   int incX = 1;
   float Y[] = { 0.637, -0.177 };
   int incY = -1;
   float A_expected[] = { -1.526889, 0.040465 };
   cblas_cgerc(order, M, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[2*i], A_expected[2*i], flteps, "cgerc(case 1143) real");
       gsl_test_rel(A[2*i+1], A_expected[2*i+1], flteps, "cgerc(case 1143) imag");
     };
   };
  };


  {
   int order = 101;
   int M = 1;
   int N = 1;
   int lda = 1;
   double alpha[2] = {1, 0};
   double A[] = { 0.662, -0.887 };
   double X[] = { 0.261, 0.791 };
   int incX = 1;
   double Y[] = { 0.66, -0.107 };
   int incY = -1;
   double A_expected[] = { 0.918897, -0.392867 };
   cblas_zgeru(order, M, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[2*i], A_expected[2*i], dbleps, "zgeru(case 1144) real");
       gsl_test_rel(A[2*i+1], A_expected[2*i+1], dbleps, "zgeru(case 1144) imag");
     };
   };
  };


  {
   int order = 101;
   int M = 1;
   int N = 1;
   int lda = 1;
   double alpha[2] = {1, 0};
   double A[] = { 0.662, -0.887 };
   double X[] = { 0.261, 0.791 };
   int incX = 1;
   double Y[] = { 0.66, -0.107 };
   int incY = -1;
   double A_expected[] = { 0.749623, -0.337013 };
   cblas_zgerc(order, M, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[2*i], A_expected[2*i], dbleps, "zgerc(case 1145) real");
       gsl_test_rel(A[2*i+1], A_expected[2*i+1], dbleps, "zgerc(case 1145) imag");
     };
   };
  };


  {
   int order = 102;
   int M = 1;
   int N = 1;
   int lda = 1;
   double alpha[2] = {1, 0};
   double A[] = { 0.662, -0.887 };
   double X[] = { 0.261, 0.791 };
   int incX = 1;
   double Y[] = { 0.66, -0.107 };
   int incY = -1;
   double A_expected[] = { 0.918897, -0.392867 };
   cblas_zgeru(order, M, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[2*i], A_expected[2*i], dbleps, "zgeru(case 1146) real");
       gsl_test_rel(A[2*i+1], A_expected[2*i+1], dbleps, "zgeru(case 1146) imag");
     };
   };
  };


  {
   int order = 102;
   int M = 1;
   int N = 1;
   int lda = 1;
   double alpha[2] = {1, 0};
   double A[] = { 0.662, -0.887 };
   double X[] = { 0.261, 0.791 };
   int incX = 1;
   double Y[] = { 0.66, -0.107 };
   int incY = -1;
   double A_expected[] = { 0.749623, -0.337013 };
   cblas_zgerc(order, M, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[2*i], A_expected[2*i], dbleps, "zgerc(case 1147) real");
       gsl_test_rel(A[2*i+1], A_expected[2*i+1], dbleps, "zgerc(case 1147) imag");
     };
   };
  };


}
