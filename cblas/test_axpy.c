#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_axpy () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int N = 1;
   float alpha = 0.1;
   float X[] = { -0.434 };
   int incX = 1;
   float Y[] = { -0.402 };
   int incY = -1;
   float expected[] = { -0.4454 };
   cblas_saxpy(N, alpha, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], expected[i], flteps, "saxpy(case 22)");
     }
   };
  };


  {
   int N = 1;
   double alpha = 0;
   double X[] = { 0.899 };
   int incX = 1;
   double Y[] = { -0.113 };
   int incY = -1;
   double expected[] = { -0.113 };
   cblas_daxpy(N, alpha, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], expected[i], dbleps, "daxpy(case 23)");
     }
   };
  };


  {
   int N = 1;
   float alpha[2] = {0, 1};
   float X[] = { -0.897, -0.204 };
   int incX = 1;
   float Y[] = { -0.759, 0.557 };
   int incY = -1;
   float expected[] = { -0.555, -0.34 };
   cblas_caxpy(N, alpha, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], expected[2*i], flteps, "caxpy(case 24) real");
       gsl_test_rel(Y[2*i+1], expected[2*i+1], flteps, "caxpy(case 24) imag");
     };
   };
  };


  {
   int N = 1;
   double alpha[2] = {0, 0.1};
   double X[] = { 0.071, 0.776 };
   int incX = 1;
   double Y[] = { 0.983, 0.549 };
   int incY = -1;
   double expected[] = { 0.9054, 0.5561 };
   cblas_zaxpy(N, alpha, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], expected[2*i], dbleps, "zaxpy(case 25) real");
       gsl_test_rel(Y[2*i+1], expected[2*i+1], dbleps, "zaxpy(case 25) imag");
     };
   };
  };


}
