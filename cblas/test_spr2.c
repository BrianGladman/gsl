#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_spr2 () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int uplo = 121;
   int N = 1;
   float alpha = 0;
   float Ap[] = { -0.232 };
   float X[] = { -0.744 };
   int incX = 1;
   float Y[] = { -0.064 };
   int incY = -1;
   float Ap_expected[] = { -0.232 };
   cblas_sspr2(order, uplo, N, alpha, X, incX, Y, incY, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[i], Ap_expected[i], flteps, "sspr2(case 1188)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int N = 1;
   float alpha = 0;
   float Ap[] = { -0.232 };
   float X[] = { -0.744 };
   int incX = 1;
   float Y[] = { -0.064 };
   int incY = -1;
   float Ap_expected[] = { -0.232 };
   cblas_sspr2(order, uplo, N, alpha, X, incX, Y, incY, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[i], Ap_expected[i], flteps, "sspr2(case 1189)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int N = 1;
   float alpha = 0;
   float Ap[] = { -0.232 };
   float X[] = { -0.744 };
   int incX = 1;
   float Y[] = { -0.064 };
   int incY = -1;
   float Ap_expected[] = { -0.232 };
   cblas_sspr2(order, uplo, N, alpha, X, incX, Y, incY, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[i], Ap_expected[i], flteps, "sspr2(case 1190)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int N = 1;
   float alpha = 0;
   float Ap[] = { -0.232 };
   float X[] = { -0.744 };
   int incX = 1;
   float Y[] = { -0.064 };
   int incY = -1;
   float Ap_expected[] = { -0.232 };
   cblas_sspr2(order, uplo, N, alpha, X, incX, Y, incY, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[i], Ap_expected[i], flteps, "sspr2(case 1191)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int N = 1;
   double alpha = 1;
   double Ap[] = { 0.111 };
   double X[] = { -0.415 };
   int incX = 1;
   double Y[] = { -0.258 };
   int incY = -1;
   double Ap_expected[] = { 0.32514 };
   cblas_dspr2(order, uplo, N, alpha, X, incX, Y, incY, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[i], Ap_expected[i], dbleps, "dspr2(case 1192)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int N = 1;
   double alpha = 1;
   double Ap[] = { 0.111 };
   double X[] = { -0.415 };
   int incX = 1;
   double Y[] = { -0.258 };
   int incY = -1;
   double Ap_expected[] = { 0.32514 };
   cblas_dspr2(order, uplo, N, alpha, X, incX, Y, incY, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[i], Ap_expected[i], dbleps, "dspr2(case 1193)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int N = 1;
   double alpha = 1;
   double Ap[] = { 0.111 };
   double X[] = { -0.415 };
   int incX = 1;
   double Y[] = { -0.258 };
   int incY = -1;
   double Ap_expected[] = { 0.32514 };
   cblas_dspr2(order, uplo, N, alpha, X, incX, Y, incY, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[i], Ap_expected[i], dbleps, "dspr2(case 1194)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int N = 1;
   double alpha = 1;
   double Ap[] = { 0.111 };
   double X[] = { -0.415 };
   int incX = 1;
   double Y[] = { -0.258 };
   int incY = -1;
   double Ap_expected[] = { 0.32514 };
   cblas_dspr2(order, uplo, N, alpha, X, incX, Y, incY, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[i], Ap_expected[i], dbleps, "dspr2(case 1195)");
     }
   };
  };


}
