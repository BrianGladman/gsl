#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_hpr () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int uplo = 121;
   int N = 1;
   float alpha = 1;
   float Ap[] = { 0.276, -0.084 };
   float X[] = { 0.082, -0.394 };
   int incX = 1;
   float Ap_expected[] = { 0.43796, 0 };
   cblas_chpr(order, uplo, N, alpha, X, incX, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[2*i], Ap_expected[2*i], flteps, "chpr(case 1164) real");
       gsl_test_rel(Ap[2*i+1], Ap_expected[2*i+1], flteps, "chpr(case 1164) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int N = 1;
   float alpha = 1;
   float Ap[] = { 0.276, -0.084 };
   float X[] = { 0.082, -0.394 };
   int incX = 1;
   float Ap_expected[] = { 0.43796, 0 };
   cblas_chpr(order, uplo, N, alpha, X, incX, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[2*i], Ap_expected[2*i], flteps, "chpr(case 1165) real");
       gsl_test_rel(Ap[2*i+1], Ap_expected[2*i+1], flteps, "chpr(case 1165) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int N = 1;
   float alpha = 1;
   float Ap[] = { 0.276, -0.084 };
   float X[] = { 0.082, -0.394 };
   int incX = 1;
   float Ap_expected[] = { 0.43796, 0 };
   cblas_chpr(order, uplo, N, alpha, X, incX, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[2*i], Ap_expected[2*i], flteps, "chpr(case 1166) real");
       gsl_test_rel(Ap[2*i+1], Ap_expected[2*i+1], flteps, "chpr(case 1166) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int N = 1;
   float alpha = 1;
   float Ap[] = { 0.276, -0.084 };
   float X[] = { 0.082, -0.394 };
   int incX = 1;
   float Ap_expected[] = { 0.43796, 0 };
   cblas_chpr(order, uplo, N, alpha, X, incX, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[2*i], Ap_expected[2*i], flteps, "chpr(case 1167) real");
       gsl_test_rel(Ap[2*i+1], Ap_expected[2*i+1], flteps, "chpr(case 1167) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int N = 1;
   double alpha = -0.3;
   double Ap[] = { 0.923, -0.788 };
   double X[] = { 0.775, 0.691 };
   int incX = 1;
   double Ap_expected[] = { 0.5995682, 0 };
   cblas_zhpr(order, uplo, N, alpha, X, incX, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[2*i], Ap_expected[2*i], dbleps, "zhpr(case 1168) real");
       gsl_test_rel(Ap[2*i+1], Ap_expected[2*i+1], dbleps, "zhpr(case 1168) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int N = 1;
   double alpha = -0.3;
   double Ap[] = { 0.923, -0.788 };
   double X[] = { 0.775, 0.691 };
   int incX = 1;
   double Ap_expected[] = { 0.5995682, 0 };
   cblas_zhpr(order, uplo, N, alpha, X, incX, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[2*i], Ap_expected[2*i], dbleps, "zhpr(case 1169) real");
       gsl_test_rel(Ap[2*i+1], Ap_expected[2*i+1], dbleps, "zhpr(case 1169) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int N = 1;
   double alpha = -0.3;
   double Ap[] = { 0.923, -0.788 };
   double X[] = { 0.775, 0.691 };
   int incX = 1;
   double Ap_expected[] = { 0.5995682, 0 };
   cblas_zhpr(order, uplo, N, alpha, X, incX, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[2*i], Ap_expected[2*i], dbleps, "zhpr(case 1170) real");
       gsl_test_rel(Ap[2*i+1], Ap_expected[2*i+1], dbleps, "zhpr(case 1170) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int N = 1;
   double alpha = -0.3;
   double Ap[] = { 0.923, -0.788 };
   double X[] = { 0.775, 0.691 };
   int incX = 1;
   double Ap_expected[] = { 0.5995682, 0 };
   cblas_zhpr(order, uplo, N, alpha, X, incX, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[2*i], Ap_expected[2*i], dbleps, "zhpr(case 1171) real");
       gsl_test_rel(Ap[2*i+1], Ap_expected[2*i+1], dbleps, "zhpr(case 1171) imag");
     };
   };
  };


}
