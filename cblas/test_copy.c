#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

#include "tests.h"

void
test_copy (void) {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int N = 1;
   float X[] = { 0.898 };
   int incX = 1;
   float Y[] = { 0.699 };
   int incY = -1;
   float expected[] = { 0.898 };
   cblas_scopy(N, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], expected[i], flteps, "scopy(case 76)");
     }
   };
  };


  {
   int N = 1;
   double X[] = { 0.002 };
   int incX = 1;
   double Y[] = { -0.921 };
   int incY = -1;
   double expected[] = { 0.002 };
   cblas_dcopy(N, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], expected[i], dbleps, "dcopy(case 77)");
     }
   };
  };


  {
   int N = 1;
   float X[] = { -0.166, 0.639 };
   int incX = 1;
   float Y[] = { 0.863, 0.613 };
   int incY = -1;
   float expected[] = { -0.166, 0.639 };
   cblas_ccopy(N, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], expected[2*i], flteps, "ccopy(case 78) real");
       gsl_test_rel(Y[2*i+1], expected[2*i+1], flteps, "ccopy(case 78) imag");
     };
   };
  };


  {
   int N = 1;
   double X[] = { 0.315, -0.324 };
   int incX = 1;
   double Y[] = { -0.312, -0.748 };
   int incY = -1;
   double expected[] = { 0.315, -0.324 };
   cblas_zcopy(N, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], expected[2*i], dbleps, "zcopy(case 79) real");
       gsl_test_rel(Y[2*i+1], expected[2*i+1], dbleps, "zcopy(case 79) imag");
     };
   };
  };


  {
   int N = 1;
   float X[] = { 0.222 };
   int incX = -1;
   float Y[] = { 0.522 };
   int incY = 1;
   float expected[] = { 0.222 };
   cblas_scopy(N, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], expected[i], flteps, "scopy(case 80)");
     }
   };
  };


  {
   int N = 1;
   double X[] = { 0.021 };
   int incX = -1;
   double Y[] = { 0.898 };
   int incY = 1;
   double expected[] = { 0.021 };
   cblas_dcopy(N, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], expected[i], dbleps, "dcopy(case 81)");
     }
   };
  };


  {
   int N = 1;
   float X[] = { 0.376, 0.229 };
   int incX = -1;
   float Y[] = { 0.143, -0.955 };
   int incY = 1;
   float expected[] = { 0.376, 0.229 };
   cblas_ccopy(N, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], expected[2*i], flteps, "ccopy(case 82) real");
       gsl_test_rel(Y[2*i+1], expected[2*i+1], flteps, "ccopy(case 82) imag");
     };
   };
  };


  {
   int N = 1;
   double X[] = { -0.265, -0.84 };
   int incX = -1;
   double Y[] = { -0.156, 0.939 };
   int incY = 1;
   double expected[] = { -0.265, -0.84 };
   cblas_zcopy(N, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], expected[2*i], dbleps, "zcopy(case 83) real");
       gsl_test_rel(Y[2*i+1], expected[2*i+1], dbleps, "zcopy(case 83) imag");
     };
   };
  };


  {
   int N = 1;
   float X[] = { 0.074 };
   int incX = -1;
   float Y[] = { -0.802 };
   int incY = -1;
   float expected[] = { 0.074 };
   cblas_scopy(N, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], expected[i], flteps, "scopy(case 84)");
     }
   };
  };


  {
   int N = 1;
   double X[] = { -0.374 };
   int incX = -1;
   double Y[] = { -0.161 };
   int incY = -1;
   double expected[] = { -0.374 };
   cblas_dcopy(N, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], expected[i], dbleps, "dcopy(case 85)");
     }
   };
  };


  {
   int N = 1;
   float X[] = { 0.084, 0.778 };
   int incX = -1;
   float Y[] = { 0.31, -0.797 };
   int incY = -1;
   float expected[] = { 0.084, 0.778 };
   cblas_ccopy(N, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], expected[2*i], flteps, "ccopy(case 86) real");
       gsl_test_rel(Y[2*i+1], expected[2*i+1], flteps, "ccopy(case 86) imag");
     };
   };
  };


  {
   int N = 1;
   double X[] = { 0.831, -0.282 };
   int incX = -1;
   double Y[] = { -0.62, 0.32 };
   int incY = -1;
   double expected[] = { 0.831, -0.282 };
   cblas_zcopy(N, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], expected[2*i], dbleps, "zcopy(case 87) real");
       gsl_test_rel(Y[2*i+1], expected[2*i+1], dbleps, "zcopy(case 87) imag");
     };
   };
  };


}
