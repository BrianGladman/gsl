#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_dot () {
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


  {
   int N = 1;
   float X[] = { -0.812 };
   float Y[] = { -0.667 };
   int incX = 1;
   int incY = -1;
   float expected = 0.541604;
   float f;
   f = cblas_sdot(N, X, incX, Y, incY);
   gsl_test_rel(f, expected, flteps, "sdot(case 4)");
  };


  {
   int N = 1;
   double X[] = { 0.481 };
   double Y[] = { 0.523 };
   int incX = 1;
   int incY = -1;
   double expected = 0.251563;
   double f;
   f = cblas_ddot(N, X, incX, Y, incY);
   gsl_test_rel(f, expected, dbleps, "ddot(case 5)");
  };


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
   double expected[2] = {-0.174096, -0.146928};
   double f[2];
   cblas_zdotu_sub(N, X, incX, Y, incY, &f);
   gsl_test_rel(f[0], expected[0], dbleps, "zdotu(case 8) real");
   gsl_test_rel(f[1], expected[1], dbleps, "zdotu(case 8) imag");
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
