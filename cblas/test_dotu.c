#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_dotu () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int N = 1;
   float X[] = { 0.785, -0.7 };
   float Y[] = { 0.79, -0.679 };
   int incX = 1;
   int incY = -1;
   float expected[2] = {0.14485, -1.086015};
   float f[2];
   cblas_cdotu_sub(N, X, incX, Y, incY, &f);
   gsl_test_rel(f[0], expected[0], flteps, "cdotu(case 6) real");
   gsl_test_rel(f[1], expected[1], flteps, "cdotu(case 6) imag");
  };


  {
   int N = 1;
   double X[] = { 0.474, -0.27 };
   double Y[] = { -0.144, -0.392 };
   int incX = 1;
   int incY = -1;
   double expected[2] = {-0.174096, -0.146928};
   double f[2];
   cblas_zdotu_sub(N, X, incX, Y, incY, &f);
   gsl_test_rel(f[0], expected[0], dbleps, "zdotu(case 8) real");
   gsl_test_rel(f[1], expected[1], dbleps, "zdotu(case 8) imag");
  };


}
