#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_hpr2 () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int uplo = 121;
   int N = 1;
   float alpha[2] = {1, 0};
   float Ap[] = { 0.838, -0.066 };
   float X[] = { 0.982, 0.565 };
   int incX = 1;
   float Y[] = { -0.18, -0.186 };
   int incY = -1;
   float Ap_expected[] = { 0.2743, 0 };
   cblas_chpr2(order, uplo, N, alpha, X, incX, Y, incY, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[2*i], Ap_expected[2*i], flteps, "chpr2(case 1204) real");
       gsl_test_rel(Ap[2*i+1], Ap_expected[2*i+1], flteps, "chpr2(case 1204) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int N = 1;
   float alpha[2] = {1, 0};
   float Ap[] = { 0.838, -0.066 };
   float X[] = { 0.982, 0.565 };
   int incX = 1;
   float Y[] = { -0.18, -0.186 };
   int incY = -1;
   float Ap_expected[] = { 0.2743, 0 };
   cblas_chpr2(order, uplo, N, alpha, X, incX, Y, incY, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[2*i], Ap_expected[2*i], flteps, "chpr2(case 1205) real");
       gsl_test_rel(Ap[2*i+1], Ap_expected[2*i+1], flteps, "chpr2(case 1205) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int N = 1;
   float alpha[2] = {1, 0};
   float Ap[] = { 0.838, -0.066 };
   float X[] = { 0.982, 0.565 };
   int incX = 1;
   float Y[] = { -0.18, -0.186 };
   int incY = -1;
   float Ap_expected[] = { 0.2743, 0 };
   cblas_chpr2(order, uplo, N, alpha, X, incX, Y, incY, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[2*i], Ap_expected[2*i], flteps, "chpr2(case 1206) real");
       gsl_test_rel(Ap[2*i+1], Ap_expected[2*i+1], flteps, "chpr2(case 1206) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int N = 1;
   float alpha[2] = {1, 0};
   float Ap[] = { 0.838, -0.066 };
   float X[] = { 0.982, 0.565 };
   int incX = 1;
   float Y[] = { -0.18, -0.186 };
   int incY = -1;
   float Ap_expected[] = { 0.2743, 0 };
   cblas_chpr2(order, uplo, N, alpha, X, incX, Y, incY, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[2*i], Ap_expected[2*i], flteps, "chpr2(case 1207) real");
       gsl_test_rel(Ap[2*i+1], Ap_expected[2*i+1], flteps, "chpr2(case 1207) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int N = 1;
   double alpha[2] = {-1, 0};
   double Ap[] = { -0.653, -0.466 };
   double X[] = { -0.013, -0.807 };
   int incX = 1;
   double Y[] = { 0.864, -0.606 };
   int incY = -1;
   double Ap_expected[] = { -1.60862, 0 };
   cblas_zhpr2(order, uplo, N, alpha, X, incX, Y, incY, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[2*i], Ap_expected[2*i], dbleps, "zhpr2(case 1208) real");
       gsl_test_rel(Ap[2*i+1], Ap_expected[2*i+1], dbleps, "zhpr2(case 1208) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int N = 1;
   double alpha[2] = {-1, 0};
   double Ap[] = { -0.653, -0.466 };
   double X[] = { -0.013, -0.807 };
   int incX = 1;
   double Y[] = { 0.864, -0.606 };
   int incY = -1;
   double Ap_expected[] = { -1.60862, 0 };
   cblas_zhpr2(order, uplo, N, alpha, X, incX, Y, incY, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[2*i], Ap_expected[2*i], dbleps, "zhpr2(case 1209) real");
       gsl_test_rel(Ap[2*i+1], Ap_expected[2*i+1], dbleps, "zhpr2(case 1209) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int N = 1;
   double alpha[2] = {-1, 0};
   double Ap[] = { -0.653, -0.466 };
   double X[] = { -0.013, -0.807 };
   int incX = 1;
   double Y[] = { 0.864, -0.606 };
   int incY = -1;
   double Ap_expected[] = { -1.60862, 0 };
   cblas_zhpr2(order, uplo, N, alpha, X, incX, Y, incY, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[2*i], Ap_expected[2*i], dbleps, "zhpr2(case 1210) real");
       gsl_test_rel(Ap[2*i+1], Ap_expected[2*i+1], dbleps, "zhpr2(case 1210) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int N = 1;
   double alpha[2] = {-1, 0};
   double Ap[] = { -0.653, -0.466 };
   double X[] = { -0.013, -0.807 };
   int incX = 1;
   double Y[] = { 0.864, -0.606 };
   int incY = -1;
   double Ap_expected[] = { -1.60862, 0 };
   cblas_zhpr2(order, uplo, N, alpha, X, incX, Y, incY, Ap);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Ap[2*i], Ap_expected[2*i], dbleps, "zhpr2(case 1211) real");
       gsl_test_rel(Ap[2*i+1], Ap_expected[2*i+1], dbleps, "zhpr2(case 1211) imag");
     };
   };
  };


}
