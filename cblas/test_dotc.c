#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_dotc () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int N = 1;
   float X[] = { 0.785, -0.7 };
   float Y[] = { 0.79, -0.679 };
   int incX = 1;
   int incY = -1;
   float expected[2] = {1.09545, 0.019985};
   float f[2];
   cblas_cdotc_sub(N, X, incX, Y, incY, &f);
   gsl_test_rel(f[0], expected[0], flteps, "cdotc(case 7) real");
   gsl_test_rel(f[1], expected[1], flteps, "cdotc(case 7) imag");
  };


  {
   int N = 1;
   double X[] = { 0.474, -0.27 };
   double Y[] = { -0.144, -0.392 };
   int incX = 1;
   int incY = -1;
   double expected[2] = {0.037584, -0.224688};
   double f[2];
   cblas_zdotc_sub(N, X, incX, Y, incY, &f);
   gsl_test_rel(f[0], expected[0], dbleps, "zdotc(case 9) real");
   gsl_test_rel(f[1], expected[1], dbleps, "zdotc(case 9) imag");
  };


}
