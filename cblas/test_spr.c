#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_spr () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int uplo = 121;
   int N = 1;
   float alpha = 0.1;
   float Ap[] = { -0.263 };
   float X[] = { -0.88 };
   int incX = 1;
   float Ap_expected[] = { -0.18556 };
   cblas_sspr(order, uplo, N, alpha, X, incX, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[i], Ap_expected[i], flteps, "sspr(case 1172)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int N = 1;
   float alpha = 0.1;
   float Ap[] = { -0.263 };
   float X[] = { -0.88 };
   int incX = 1;
   float Ap_expected[] = { -0.18556 };
   cblas_sspr(order, uplo, N, alpha, X, incX, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[i], Ap_expected[i], flteps, "sspr(case 1173)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int N = 1;
   float alpha = 0.1;
   float Ap[] = { -0.263 };
   float X[] = { -0.88 };
   int incX = 1;
   float Ap_expected[] = { -0.18556 };
   cblas_sspr(order, uplo, N, alpha, X, incX, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[i], Ap_expected[i], flteps, "sspr(case 1174)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int N = 1;
   float alpha = 0.1;
   float Ap[] = { -0.263 };
   float X[] = { -0.88 };
   int incX = 1;
   float Ap_expected[] = { -0.18556 };
   cblas_sspr(order, uplo, N, alpha, X, incX, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[i], Ap_expected[i], flteps, "sspr(case 1175)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int N = 1;
   double alpha = -1;
   double Ap[] = { 0.034 };
   double X[] = { 0.407 };
   int incX = 1;
   double Ap_expected[] = { -0.131649 };
   cblas_dspr(order, uplo, N, alpha, X, incX, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[i], Ap_expected[i], dbleps, "dspr(case 1176)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int N = 1;
   double alpha = -1;
   double Ap[] = { 0.034 };
   double X[] = { 0.407 };
   int incX = 1;
   double Ap_expected[] = { -0.131649 };
   cblas_dspr(order, uplo, N, alpha, X, incX, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[i], Ap_expected[i], dbleps, "dspr(case 1177)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int N = 1;
   double alpha = -1;
   double Ap[] = { 0.034 };
   double X[] = { 0.407 };
   int incX = 1;
   double Ap_expected[] = { -0.131649 };
   cblas_dspr(order, uplo, N, alpha, X, incX, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[i], Ap_expected[i], dbleps, "dspr(case 1178)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int N = 1;
   double alpha = -1;
   double Ap[] = { 0.034 };
   double X[] = { 0.407 };
   int incX = 1;
   double Ap_expected[] = { -0.131649 };
   cblas_dspr(order, uplo, N, alpha, X, incX, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[i], Ap_expected[i], dbleps, "dspr(case 1179)");
     }
   };
  };


}
