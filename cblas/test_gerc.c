#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_gerc () {
const double flteps = 1e-4, dbleps = 1e-6;
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
