#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_asum () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int N = 1;
   float X[] = { 0.949 };
   int incX = -1;
   float expected = 0;
   float f;
   f = cblas_sasum(N, X, incX);
   gsl_test_rel(f, expected, flteps, "sasum(case 14)");
  };


  {
   int N = 1;
   double X[] = { -0.873 };
   int incX = -1;
   double expected = 0;
   double f;
   f = cblas_dasum(N, X, incX);
   gsl_test_rel(f, expected, dbleps, "dasum(case 15)");
  };


  {
   int N = 1;
   float X[] = { 0.852, -0.045 };
   int incX = -1;
   float expected = 0;
   float f;
   f = cblas_scasum(N, X, incX);
   gsl_test_rel(f, expected, flteps, "scasum(case 16)");
  };


  {
   int N = 1;
   double X[] = { 0.626, -0.164 };
   int incX = -1;
   double expected = 0;
   double f;
   f = cblas_dzasum(N, X, incX);
   gsl_test_rel(f, expected, dbleps, "dzasum(case 17)");
  };


}
