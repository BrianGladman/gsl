#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_swap () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int N = 1;
   float X[] = { 0.918 };
   int incX = 1;
   float Y[] = { -0.126 };
   int incY = -1;
   float expected1[] = { -0.126 };
   float expected2[] = { 0.918 };
   cblas_sswap(N, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], expected1[i], flteps, "sswap(case 30)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], expected2[i], flteps, "sswap(case 31)");
     }
   };
  };


  {
   int N = 1;
   double X[] = { 0.217 };
   int incX = 1;
   double Y[] = { -0.588 };
   int incY = -1;
   double expected1[] = { -0.588 };
   double expected2[] = { 0.217 };
   cblas_dswap(N, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], expected1[i], dbleps, "dswap(case 32)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], expected2[i], dbleps, "dswap(case 33)");
     }
   };
  };


  {
   int N = 1;
   float X[] = { 0.31, -0.442 };
   int incX = 1;
   float Y[] = { 0.059, 0.987 };
   int incY = -1;
   float expected1[] = { 0.059, 0.987 };
   float expected2[] = { 0.31, -0.442 };
   cblas_cswap(N, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], expected1[2*i], flteps, "cswap(case 34) real");
       gsl_test_rel(X[2*i+1], expected1[2*i+1], flteps, "cswap(case 34) imag");
     };
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], expected2[2*i], flteps, "cswap(case 35) real");
       gsl_test_rel(Y[2*i+1], expected2[2*i+1], flteps, "cswap(case 35) imag");
     };
   };
  };


  {
   int N = 1;
   double X[] = { 0.609, -0.143 };
   int incX = 1;
   double Y[] = { 0.615, -0.957 };
   int incY = -1;
   double expected1[] = { 0.615, -0.957 };
   double expected2[] = { 0.609, -0.143 };
   cblas_zswap(N, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], expected1[2*i], dbleps, "zswap(case 36) real");
       gsl_test_rel(X[2*i+1], expected1[2*i+1], dbleps, "zswap(case 36) imag");
     };
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], expected2[2*i], dbleps, "zswap(case 37) real");
       gsl_test_rel(Y[2*i+1], expected2[2*i+1], dbleps, "zswap(case 37) imag");
     };
   };
  };


}
