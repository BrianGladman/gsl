#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

#include "tests.h"

void
test_nrm2 (void) {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int N = 1;
   float X[] = { 0.317 };
   int incX = -1;
   float expected = 0;
   float f;
   f = cblas_snrm2(N, X, incX);
   gsl_test_rel(f, expected, flteps, "snrm2(case 28)");
  };


  {
   int N = 1;
   double X[] = { 0.071 };
   int incX = -1;
   double expected = 0;
   double f;
   f = cblas_dnrm2(N, X, incX);
   gsl_test_rel(f, expected, dbleps, "dnrm2(case 29)");
  };


  {
   int N = 1;
   float X[] = { 0.776, 0.983 };
   int incX = -1;
   float expected = 0;
   float f;
   f = cblas_scnrm2(N, X, incX);
   gsl_test_rel(f, expected, flteps, "scnrm2(case 30)");
  };


  {
   int N = 1;
   double X[] = { 0.549, -0.354 };
   int incX = -1;
   double expected = 0;
   double f;
   f = cblas_dznrm2(N, X, incX);
   gsl_test_rel(f, expected, dbleps, "dznrm2(case 31)");
  };


  {
   int N = 2;
   float X[] = { 0.14, -0.632 };
   int incX = 1;
   float expected = 0.647320631527;
   float f;
   f = cblas_snrm2(N, X, incX);
   gsl_test_rel(f, expected, flteps, "snrm2(case 32)");
  };


  {
   int N = 2;
   double X[] = { 0.696, -0.804 };
   int incX = 1;
   double expected = 1.06340584915;
   double f;
   f = cblas_dnrm2(N, X, incX);
   gsl_test_rel(f, expected, dbleps, "dnrm2(case 33)");
  };


  {
   int N = 2;
   float X[] = { 0.281, -0.063, 0.367, 0.232 };
   int incX = 1;
   float expected = 0.521001919382;
   float f;
   f = cblas_scnrm2(N, X, incX);
   gsl_test_rel(f, expected, flteps, "scnrm2(case 34)");
  };


  {
   int N = 2;
   double X[] = { -0.359, -0.76, -0.906, -0.108 };
   int incX = 1;
   double expected = 1.24055672986;
   double f;
   f = cblas_dznrm2(N, X, incX);
   gsl_test_rel(f, expected, dbleps, "dznrm2(case 35)");
  };


  {
   int N = 2;
   float X[] = { 0.918, -0.126 };
   int incX = -1;
   float expected = 0;
   float f;
   f = cblas_snrm2(N, X, incX);
   gsl_test_rel(f, expected, flteps, "snrm2(case 36)");
  };


  {
   int N = 2;
   double X[] = { 0.217, -0.588 };
   int incX = -1;
   double expected = 0;
   double f;
   f = cblas_dnrm2(N, X, incX);
   gsl_test_rel(f, expected, dbleps, "dnrm2(case 37)");
  };


  {
   int N = 2;
   float X[] = { 0.31, 0.059, -0.442, 0.987 };
   int incX = -1;
   float expected = 0;
   float f;
   f = cblas_scnrm2(N, X, incX);
   gsl_test_rel(f, expected, flteps, "scnrm2(case 38)");
  };


  {
   int N = 2;
   double X[] = { 0.609, 0.615, -0.143, -0.957 };
   int incX = -1;
   double expected = 0;
   double f;
   f = cblas_dznrm2(N, X, incX);
   gsl_test_rel(f, expected, dbleps, "dznrm2(case 39)");
  };


}
