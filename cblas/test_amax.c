#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_amax () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int N = 1;
   float X[] = { -0.786 };
   int incX = -1;
   float expected = 0;
   int k;
   k = cblas_isamax(N, X, incX);
   gsl_test_int(k, expected, "samax(case 18)");
  };


  {
   int N = 1;
   double X[] = { -0.341 };
   int incX = -1;
   double expected = 0;
   int k;
   k = cblas_idamax(N, X, incX);
   gsl_test_int(k, expected, "damax(case 19)");
  };


  {
   int N = 1;
   float X[] = { -0.271, -0.896 };
   int incX = -1;
   float expected = 0;
   int k;
   k = cblas_icamax(N, X, incX);
   gsl_test_int(k, expected, "camax(case 20)");
  };


  {
   int N = 1;
   double X[] = { -0.088, -0.165 };
   int incX = -1;
   double expected = 0;
   int k;
   k = cblas_izamax(N, X, incX);
   gsl_test_int(k, expected, "zamax(case 21)");
  };


}
