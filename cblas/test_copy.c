#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_copy () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int N = 1;
   float X[] = { 0.14 };
   int incX = 1;
   float Y[] = { -0.632 };
   int incY = -1;
   float expected[] = { 0.14 };
   cblas_scopy(N, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], expected[i], flteps, "scopy(case 26)");
     }
   };
  };


  {
   int N = 1;
   double X[] = { 0.696 };
   int incX = 1;
   double Y[] = { -0.804 };
   int incY = -1;
   double expected[] = { 0.696 };
   cblas_dcopy(N, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], expected[i], dbleps, "dcopy(case 27)");
     }
   };
  };


  {
   int N = 1;
   float X[] = { 0.281, 0.367 };
   int incX = 1;
   float Y[] = { -0.063, 0.232 };
   int incY = -1;
   float expected[] = { 0.281, 0.367 };
   cblas_ccopy(N, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], expected[2*i], flteps, "ccopy(case 28) real");
       gsl_test_rel(Y[2*i+1], expected[2*i+1], flteps, "ccopy(case 28) imag");
     };
   };
  };


  {
   int N = 1;
   double X[] = { -0.359, -0.906 };
   int incX = 1;
   double Y[] = { -0.76, -0.108 };
   int incY = -1;
   double expected[] = { -0.359, -0.906 };
   cblas_zcopy(N, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], expected[2*i], dbleps, "zcopy(case 29) real");
       gsl_test_rel(Y[2*i+1], expected[2*i+1], dbleps, "zcopy(case 29) imag");
     };
   };
  };


}
