#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_rot () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int N = 1;
   float c = 0;
   float s = 0;
   float X[] = { 0.899 };
   int incX = 1;
   float Y[] = { -0.72 };
   int incY = -1;
   float x_expected[] = { 0.000000000000e+00 };
   float y_expected[] = { 0.000000000000e+00 };
   cblas_srot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "srot(case 448)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "srot(case 449)");
     }
   };
  };


  {
   int N = 1;
   float c = 0.866025403784;
   float s = 0.5;
   float X[] = { 0.899 };
   int incX = 1;
   float Y[] = { -0.72 };
   int incY = -1;
   float x_expected[] = { 0.418556838002 };
   float y_expected[] = { -1.07303829072 };
   cblas_srot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "srot(case 450)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "srot(case 451)");
     }
   };
  };


  {
   int N = 1;
   float c = 0;
   float s = -1;
   float X[] = { 0.899 };
   int incX = 1;
   float Y[] = { -0.72 };
   int incY = -1;
   float x_expected[] = { 0.72 };
   float y_expected[] = { 0.899 };
   cblas_srot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "srot(case 452)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "srot(case 453)");
     }
   };
  };


  {
   int N = 1;
   float c = -1;
   float s = 0;
   float X[] = { 0.899 };
   int incX = 1;
   float Y[] = { -0.72 };
   int incY = -1;
   float x_expected[] = { -0.899 };
   float y_expected[] = { 0.72 };
   cblas_srot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "srot(case 454)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "srot(case 455)");
     }
   };
  };


  {
   int N = 1;
   double c = 0;
   double s = 0;
   double X[] = { 0.271 };
   int incX = 1;
   double Y[] = { -0.012 };
   int incY = -1;
   double x_expected[] = { 0.000000000000e+00 };
   double y_expected[] = { 0.000000000000e+00 };
   cblas_drot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "drot(case 456)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "drot(case 457)");
     }
   };
  };


  {
   int N = 1;
   double c = 0.866025403784;
   double s = 0.5;
   double X[] = { 0.271 };
   int incX = 1;
   double Y[] = { -0.012 };
   int incY = -1;
   double x_expected[] = { 0.228692884426 };
   double y_expected[] = { -0.145892304845 };
   cblas_drot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "drot(case 458)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "drot(case 459)");
     }
   };
  };


  {
   int N = 1;
   double c = 0;
   double s = -1;
   double X[] = { 0.271 };
   int incX = 1;
   double Y[] = { -0.012 };
   int incY = -1;
   double x_expected[] = { 0.012 };
   double y_expected[] = { 0.271 };
   cblas_drot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "drot(case 460)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "drot(case 461)");
     }
   };
  };


  {
   int N = 1;
   double c = -1;
   double s = 0;
   double X[] = { 0.271 };
   int incX = 1;
   double Y[] = { -0.012 };
   int incY = -1;
   double x_expected[] = { -0.271 };
   double y_expected[] = { 0.012 };
   cblas_drot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "drot(case 462)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "drot(case 463)");
     }
   };
  };


}
