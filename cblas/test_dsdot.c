#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_dsdot () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int N = 1;
   float alpha = 0;
   float X[] = { 0.733 };
   float Y[] = { 0.825 };
   int incX = 1;
   int incY = -1;
   float expected = 0.604725;
   float f;
   f = cblas_sdsdot (N, alpha, X, incX, Y, incY);
   gsl_test_rel(f, expected, flteps, "sdsdot(case 1)");
  };


  {
   int N = 1;
   float alpha = 0.1;
   float X[] = { 0.733 };
   float Y[] = { 0.825 };
   int incX = 1;
   int incY = -1;
   float expected = 0.704725;
   float f;
   f = cblas_sdsdot (N, alpha, X, incX, Y, incY);
   gsl_test_rel(f, expected, flteps, "sdsdot(case 2)");
  };


  {
   int N = 1;
   float alpha = 1;
   float X[] = { 0.733 };
   float Y[] = { 0.825 };
   int incX = 1;
   int incY = -1;
   float expected = 1.604725;
   float f;
   f = cblas_sdsdot (N, alpha, X, incX, Y, incY);
   gsl_test_rel(f, expected, flteps, "sdsdot(case 3)");
  };


}
