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
   int N = 2;
   float alpha = 0.1;
   float Ap[] = { -0.273, -0.499, -0.305, -0.277, 0.238, -0.369 };
   float X[] = { 0.638, -0.905, 0.224, 0.182 };
   int incX = -1;
   float Ap_expected[] = { -0.26467, 0, -0.30718, -0.245116, 0.360607, 0 };
   cblas_chpr(order, uplo, N, alpha, X, incX, Ap);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Ap[2*i], Ap_expected[2*i], flteps, "chpr(case 1418) real");
       gsl_test_rel(Ap[2*i+1], Ap_expected[2*i+1], flteps, "chpr(case 1418) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int N = 2;
   float alpha = 0.1;
   float Ap[] = { -0.273, -0.499, -0.305, -0.277, 0.238, -0.369 };
   float X[] = { 0.638, -0.905, 0.224, 0.182 };
   int incX = -1;
   float Ap_expected[] = { -0.26467, 0, -0.30718, -0.308884, 0.360607, 0 };
   cblas_chpr(order, uplo, N, alpha, X, incX, Ap);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Ap[2*i], Ap_expected[2*i], flteps, "chpr(case 1419) real");
       gsl_test_rel(Ap[2*i+1], Ap_expected[2*i+1], flteps, "chpr(case 1419) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int N = 2;
   float alpha = 0.1;
   float Ap[] = { -0.273, -0.499, -0.305, -0.277, 0.238, -0.369 };
   float X[] = { 0.638, -0.905, 0.224, 0.182 };
   int incX = -1;
   float Ap_expected[] = { -0.26467, 0, -0.30718, -0.245116, 0.360607, 0 };
   cblas_chpr(order, uplo, N, alpha, X, incX, Ap);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Ap[2*i], Ap_expected[2*i], flteps, "chpr(case 1420) real");
       gsl_test_rel(Ap[2*i+1], Ap_expected[2*i+1], flteps, "chpr(case 1420) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int N = 2;
   float alpha = 0.1;
   float Ap[] = { -0.273, -0.499, -0.305, -0.277, 0.238, -0.369 };
   float X[] = { 0.638, -0.905, 0.224, 0.182 };
   int incX = -1;
   float Ap_expected[] = { -0.26467, 0, -0.30718, -0.308884, 0.360607, 0 };
   cblas_chpr(order, uplo, N, alpha, X, incX, Ap);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Ap[2*i], Ap_expected[2*i], flteps, "chpr(case 1421) real");
       gsl_test_rel(Ap[2*i+1], Ap_expected[2*i+1], flteps, "chpr(case 1421) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int N = 2;
   double alpha = 1;
   double Ap[] = { 0.265, 0.362, -0.855, 0.035, 0.136, 0.133 };
   double X[] = { -0.278, -0.686, -0.736, -0.918 };
   int incX = -1;
   double Ap_expected[] = { 1.64942, 0, -0.020644, -0.214692, 0.68388, 0 };
   cblas_zhpr(order, uplo, N, alpha, X, incX, Ap);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Ap[2*i], Ap_expected[2*i], dbleps, "zhpr(case 1422) real");
       gsl_test_rel(Ap[2*i+1], Ap_expected[2*i+1], dbleps, "zhpr(case 1422) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int N = 2;
   double alpha = 1;
   double Ap[] = { 0.265, 0.362, -0.855, 0.035, 0.136, 0.133 };
   double X[] = { -0.278, -0.686, -0.736, -0.918 };
   int incX = -1;
   double Ap_expected[] = { 1.64942, 0, -0.020644, 0.284692, 0.68388, 0 };
   cblas_zhpr(order, uplo, N, alpha, X, incX, Ap);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Ap[2*i], Ap_expected[2*i], dbleps, "zhpr(case 1423) real");
       gsl_test_rel(Ap[2*i+1], Ap_expected[2*i+1], dbleps, "zhpr(case 1423) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int N = 2;
   double alpha = 1;
   double Ap[] = { 0.265, 0.362, -0.855, 0.035, 0.136, 0.133 };
   double X[] = { -0.278, -0.686, -0.736, -0.918 };
   int incX = -1;
   double Ap_expected[] = { 1.64942, 0, -0.020644, -0.214692, 0.68388, 0 };
   cblas_zhpr(order, uplo, N, alpha, X, incX, Ap);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Ap[2*i], Ap_expected[2*i], dbleps, "zhpr(case 1424) real");
       gsl_test_rel(Ap[2*i+1], Ap_expected[2*i+1], dbleps, "zhpr(case 1424) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int N = 2;
   double alpha = 1;
   double Ap[] = { 0.265, 0.362, -0.855, 0.035, 0.136, 0.133 };
   double X[] = { -0.278, -0.686, -0.736, -0.918 };
   int incX = -1;
   double Ap_expected[] = { 1.64942, 0, -0.020644, 0.284692, 0.68388, 0 };
   cblas_zhpr(order, uplo, N, alpha, X, incX, Ap);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Ap[2*i], Ap_expected[2*i], dbleps, "zhpr(case 1425) real");
       gsl_test_rel(Ap[2*i+1], Ap_expected[2*i+1], dbleps, "zhpr(case 1425) imag");
     };
   };
  };


}
