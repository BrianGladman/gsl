#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_nrm2 () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int N = 1;
   float X[] = { -0.87 };
   int incX = -1;
   float expected = 0;
   float f;
   f = cblas_snrm2(N, X, incX);
   gsl_test_rel(f, expected, flteps, "snrm2(case 10)");
  };


  {
   int N = 1;
   double X[] = { -0.631 };
   int incX = -1;
   double expected = 0;
   double f;
   f = cblas_dnrm2(N, X, incX);
   gsl_test_rel(f, expected, dbleps, "dnrm2(case 11)");
  };


  {
   int N = 1;
   float X[] = { -0.7, -0.224 };
   int incX = -1;
   float expected = 0;
   float f;
   f = cblas_scnrm2(N, X, incX);
   gsl_test_rel(f, expected, flteps, "scnrm2(case 12)");
  };


  {
   int N = 1;
   double X[] = { -0.457, 0.839 };
   int incX = -1;
   double expected = 0;
   double f;
   f = cblas_dznrm2(N, X, incX);
   gsl_test_rel(f, expected, dbleps, "dznrm2(case 13)");
  };


}
