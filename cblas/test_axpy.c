#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_axpy () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int N = 1;
   float alpha = 0;
   float X[] = { 0.018 };
   int incX = 1;
   float Y[] = { -0.417 };
   int incY = -1;
   float expected[] = { -0.417 };
   cblas_saxpy(N, alpha, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], expected[i], flteps, "saxpy(case 64)");
     }
   };
  };


  {
   int N = 1;
   double alpha = 0;
   double X[] = { 0.071 };
   int incX = 1;
   double Y[] = { -0.888 };
   int incY = -1;
   double expected[] = { -0.888 };
   cblas_daxpy(N, alpha, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], expected[i], dbleps, "daxpy(case 65)");
     }
   };
  };


  {
   int N = 1;
   float alpha[2] = {1, 0};
   float X[] = { -0.542, -0.362 };
   int incX = 1;
   float Y[] = { -0.459, -0.433 };
   int incY = -1;
   float expected[] = { -1.001, -0.795 };
   cblas_caxpy(N, alpha, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], expected[2*i], flteps, "caxpy(case 66) real");
       gsl_test_rel(Y[2*i+1], expected[2*i+1], flteps, "caxpy(case 66) imag");
     };
   };
  };


  {
   int N = 1;
   double alpha[2] = {-1, 0};
   double X[] = { 0.003, -0.514 };
   int incX = 1;
   double Y[] = { -0.529, 0.743 };
   int incY = -1;
   double expected[] = { -0.532, 1.257 };
   cblas_zaxpy(N, alpha, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], expected[2*i], dbleps, "zaxpy(case 67) real");
       gsl_test_rel(Y[2*i+1], expected[2*i+1], dbleps, "zaxpy(case 67) imag");
     };
   };
  };


  {
   int N = 1;
   float alpha = 0.1;
   float X[] = { 0.771 };
   int incX = -1;
   float Y[] = { 0.507 };
   int incY = 1;
   float expected[] = { 0.5841 };
   cblas_saxpy(N, alpha, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], expected[i], flteps, "saxpy(case 68)");
     }
   };
  };


  {
   int N = 1;
   double alpha = -0.3;
   double X[] = { 0.029 };
   int incX = -1;
   double Y[] = { -0.992 };
   int incY = 1;
   double expected[] = { -1.0007 };
   cblas_daxpy(N, alpha, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], expected[i], dbleps, "daxpy(case 69)");
     }
   };
  };


  {
   int N = 1;
   float alpha[2] = {-0.3, 0.1};
   float X[] = { 0.194, -0.959 };
   int incX = -1;
   float Y[] = { 0.096, 0.032 };
   int incY = 1;
   float expected[] = { 0.1337, 0.3391 };
   cblas_caxpy(N, alpha, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], expected[2*i], flteps, "caxpy(case 70) real");
       gsl_test_rel(Y[2*i+1], expected[2*i+1], flteps, "caxpy(case 70) imag");
     };
   };
  };


  {
   int N = 1;
   double alpha[2] = {0, 1};
   double X[] = { 0.776, -0.671 };
   int incX = -1;
   double Y[] = { 0.39, 0.404 };
   int incY = 1;
   double expected[] = { 1.061, 1.18 };
   cblas_zaxpy(N, alpha, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], expected[2*i], dbleps, "zaxpy(case 71) real");
       gsl_test_rel(Y[2*i+1], expected[2*i+1], dbleps, "zaxpy(case 71) imag");
     };
   };
  };


  {
   int N = 1;
   float alpha = 1;
   float X[] = { 0.647 };
   int incX = -1;
   float Y[] = { 0.016 };
   int incY = -1;
   float expected[] = { 0.663 };
   cblas_saxpy(N, alpha, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], expected[i], flteps, "saxpy(case 72)");
     }
   };
  };


  {
   int N = 1;
   double alpha = -1;
   double X[] = { -0.558 };
   int incX = -1;
   double Y[] = { 0.308 };
   int incY = -1;
   double expected[] = { 0.866 };
   cblas_daxpy(N, alpha, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], expected[i], dbleps, "daxpy(case 73)");
     }
   };
  };


  {
   int N = 1;
   float alpha[2] = {-0.3, 0.1};
   float X[] = { 0.899, -0.624 };
   int incX = -1;
   float Y[] = { 0.155, -0.33 };
   int incY = -1;
   float expected[] = { -0.0523, -0.0529 };
   cblas_caxpy(N, alpha, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], expected[2*i], flteps, "caxpy(case 74) real");
       gsl_test_rel(Y[2*i+1], expected[2*i+1], flteps, "caxpy(case 74) imag");
     };
   };
  };


  {
   int N = 1;
   double alpha[2] = {0, 1};
   double X[] = { -0.451, 0.768 };
   int incX = -1;
   double Y[] = { 0.007, 0.732 };
   int incY = -1;
   double expected[] = { -0.761, 0.281 };
   cblas_zaxpy(N, alpha, X, incX, Y, incY);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[2*i], expected[2*i], dbleps, "zaxpy(case 75) real");
       gsl_test_rel(Y[2*i+1], expected[2*i+1], dbleps, "zaxpy(case 75) imag");
     };
   };
  };


}
