#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_scal () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int N = 1;
   float alpha = 0;
   float X[] = { 0.239 };
   int incX = -1;
   float expected[] = { 0.239 };
   cblas_sscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], expected[i], flteps, "sscal(case 38)");
     }
   };
  };


  {
   int N = 1;
   float alpha = 0.1;
   float X[] = { 0.239 };
   int incX = -1;
   float expected[] = { 0.239 };
   cblas_sscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], expected[i], flteps, "sscal(case 39)");
     }
   };
  };


  {
   int N = 1;
   float alpha = 1;
   float X[] = { 0.239 };
   int incX = -1;
   float expected[] = { 0.239 };
   cblas_sscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], expected[i], flteps, "sscal(case 40)");
     }
   };
  };


  {
   int N = 1;
   double alpha = 0;
   double X[] = { -0.413 };
   int incX = -1;
   double expected[] = { -0.413 };
   cblas_dscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], expected[i], dbleps, "dscal(case 41)");
     }
   };
  };


  {
   int N = 1;
   double alpha = 0.1;
   double X[] = { -0.413 };
   int incX = -1;
   double expected[] = { -0.413 };
   cblas_dscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], expected[i], dbleps, "dscal(case 42)");
     }
   };
  };


  {
   int N = 1;
   double alpha = 1;
   double X[] = { -0.413 };
   int incX = -1;
   double expected[] = { -0.413 };
   cblas_dscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], expected[i], dbleps, "dscal(case 43)");
     }
   };
  };


  {
   int N = 1;
   float alpha[2] = {0, 0};
   float X[] = { 0.1, 0.017 };
   int incX = -1;
   float expected[] = { 0.1, 0.017 };
   cblas_cscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], expected[2*i], flteps, "cscal(case 44) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], flteps, "cscal(case 44) imag");
     };
   };
  };


  {
   int N = 1;
   float alpha[2] = {0.1, 0};
   float X[] = { 0.1, 0.017 };
   int incX = -1;
   float expected[] = { 0.1, 0.017 };
   cblas_cscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], expected[2*i], flteps, "cscal(case 45) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], flteps, "cscal(case 45) imag");
     };
   };
  };


  {
   int N = 1;
   float alpha[2] = {1, 0};
   float X[] = { 0.1, 0.017 };
   int incX = -1;
   float expected[] = { 0.1, 0.017 };
   cblas_cscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], expected[2*i], flteps, "cscal(case 46) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], flteps, "cscal(case 46) imag");
     };
   };
  };


  {
   int N = 1;
   float alpha[2] = {0, 0.1};
   float X[] = { 0.1, 0.017 };
   int incX = -1;
   float expected[] = { 0.1, 0.017 };
   cblas_cscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], expected[2*i], flteps, "cscal(case 47) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], flteps, "cscal(case 47) imag");
     };
   };
  };


  {
   int N = 1;
   float alpha[2] = {0.1, 0.2};
   float X[] = { 0.1, 0.017 };
   int incX = -1;
   float expected[] = { 0.1, 0.017 };
   cblas_cscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], expected[2*i], flteps, "cscal(case 48) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], flteps, "cscal(case 48) imag");
     };
   };
  };


  {
   int N = 1;
   float alpha[2] = {1, 0.3};
   float X[] = { 0.1, 0.017 };
   int incX = -1;
   float expected[] = { 0.1, 0.017 };
   cblas_cscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], expected[2*i], flteps, "cscal(case 49) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], flteps, "cscal(case 49) imag");
     };
   };
  };


  {
   int N = 1;
   double alpha[2] = {0, 0};
   double X[] = { -0.651, 0.079 };
   int incX = -1;
   double expected[] = { -0.651, 0.079 };
   cblas_zscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], expected[2*i], dbleps, "zscal(case 50) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], dbleps, "zscal(case 50) imag");
     };
   };
  };


  {
   int N = 1;
   double alpha[2] = {0.1, 0};
   double X[] = { -0.651, 0.079 };
   int incX = -1;
   double expected[] = { -0.651, 0.079 };
   cblas_zscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], expected[2*i], dbleps, "zscal(case 51) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], dbleps, "zscal(case 51) imag");
     };
   };
  };


  {
   int N = 1;
   double alpha[2] = {1, 0};
   double X[] = { -0.651, 0.079 };
   int incX = -1;
   double expected[] = { -0.651, 0.079 };
   cblas_zscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], expected[2*i], dbleps, "zscal(case 52) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], dbleps, "zscal(case 52) imag");
     };
   };
  };


  {
   int N = 1;
   double alpha[2] = {0, 0.1};
   double X[] = { -0.651, 0.079 };
   int incX = -1;
   double expected[] = { -0.651, 0.079 };
   cblas_zscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], expected[2*i], dbleps, "zscal(case 53) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], dbleps, "zscal(case 53) imag");
     };
   };
  };


  {
   int N = 1;
   double alpha[2] = {0.1, 0.2};
   double X[] = { -0.651, 0.079 };
   int incX = -1;
   double expected[] = { -0.651, 0.079 };
   cblas_zscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], expected[2*i], dbleps, "zscal(case 54) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], dbleps, "zscal(case 54) imag");
     };
   };
  };


  {
   int N = 1;
   double alpha[2] = {1, 0.3};
   double X[] = { -0.651, 0.079 };
   int incX = -1;
   double expected[] = { -0.651, 0.079 };
   cblas_zscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], expected[2*i], dbleps, "zscal(case 55) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], dbleps, "zscal(case 55) imag");
     };
   };
  };


}
