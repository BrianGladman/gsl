#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

#include "tests.h"

void
test_hpmv (void) {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int uplo = 121;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {0, 1};
   int N = 2;
   float A[] = { 0.339, -0.102, 0.908, 0.097, -0.808, 0.236 };
   float X[] = { 0.993, -0.502, -0.653, 0.796 };
   int incX = -1;
   float Y[] = { -0.35, 0.339, -0.269, -0.122 };
   int incY = -1;
   float y_expected[] = { -0.0627557, -0.839323, -0.0877262, -0.169208 };
   cblas_chpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chpmv(case 1118) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chpmv(case 1118) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {0, 1};
   int N = 2;
   float A[] = { 0.339, -0.102, 0.908, 0.097, -0.808, 0.236 };
   float X[] = { 0.993, -0.502, -0.653, 0.796 };
   int incX = -1;
   float Y[] = { -0.35, 0.339, -0.269, -0.122 };
   int incY = -1;
   float y_expected[] = { -0.0627557, -0.839323, -0.0877262, -0.169208 };
   cblas_chpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chpmv(case 1119) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chpmv(case 1119) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {0, 1};
   int N = 2;
   float A[] = { 0.339, -0.102, 0.908, 0.097, -0.808, 0.236 };
   float X[] = { 0.993, -0.502, -0.653, 0.796 };
   int incX = -1;
   float Y[] = { -0.35, 0.339, -0.269, -0.122 };
   int incY = -1;
   float y_expected[] = { -0.0037603, -0.816761, -0.0392456, -0.121154 };
   cblas_chpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chpmv(case 1120) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chpmv(case 1120) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {0, 1};
   int N = 2;
   float A[] = { 0.339, -0.102, 0.908, 0.097, -0.808, 0.236 };
   float X[] = { 0.993, -0.502, -0.653, 0.796 };
   int incX = -1;
   float Y[] = { -0.35, 0.339, -0.269, -0.122 };
   int incY = -1;
   float y_expected[] = { -0.0037603, -0.816761, -0.0392456, -0.121154 };
   cblas_chpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chpmv(case 1121) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chpmv(case 1121) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {0, 1};
   int N = 2;
   float A[] = { 0.339, -0.102, 0.908, 0.097, -0.808, 0.236 };
   float X[] = { 0.993, -0.502, -0.653, 0.796 };
   int incX = -1;
   float Y[] = { -0.35, 0.339, -0.269, -0.122 };
   int incY = -1;
   float y_expected[] = { -0.0627557, -0.839323, -0.0877262, -0.169208 };
   cblas_chpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chpmv(case 1122) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chpmv(case 1122) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {0, 1};
   int N = 2;
   float A[] = { 0.339, -0.102, 0.908, 0.097, -0.808, 0.236 };
   float X[] = { 0.993, -0.502, -0.653, 0.796 };
   int incX = -1;
   float Y[] = { -0.35, 0.339, -0.269, -0.122 };
   int incY = -1;
   float y_expected[] = { -0.0627557, -0.839323, -0.0877262, -0.169208 };
   cblas_chpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chpmv(case 1123) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chpmv(case 1123) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {0, 1};
   int N = 2;
   float A[] = { 0.339, -0.102, 0.908, 0.097, -0.808, 0.236 };
   float X[] = { 0.993, -0.502, -0.653, 0.796 };
   int incX = -1;
   float Y[] = { -0.35, 0.339, -0.269, -0.122 };
   int incY = -1;
   float y_expected[] = { -0.0037603, -0.816761, -0.0392456, -0.121154 };
   cblas_chpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chpmv(case 1124) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chpmv(case 1124) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {0, 1};
   int N = 2;
   float A[] = { 0.339, -0.102, 0.908, 0.097, -0.808, 0.236 };
   float X[] = { 0.993, -0.502, -0.653, 0.796 };
   int incX = -1;
   float Y[] = { -0.35, 0.339, -0.269, -0.122 };
   int incY = -1;
   float y_expected[] = { -0.0037603, -0.816761, -0.0392456, -0.121154 };
   cblas_chpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chpmv(case 1125) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chpmv(case 1125) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   double alpha[2] = {1, 0};
   double beta[2] = {0, 1};
   int N = 2;
   double A[] = { 0.543, -0.737, 0.281, -0.053, -0.098, -0.826 };
   double X[] = { 0.67, -0.857, -0.613, -0.927 };
   int incX = -1;
   double Y[] = { -0.398, -0.934, -0.204, 0.183 };
   int incY = -1;
   double y_expected[] = { 0.745218, -0.60699, -0.37301, -0.983688 };
   cblas_zhpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhpmv(case 1126) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhpmv(case 1126) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   double alpha[2] = {1, 0};
   double beta[2] = {0, 1};
   int N = 2;
   double A[] = { 0.543, -0.737, 0.281, -0.053, -0.098, -0.826 };
   double X[] = { 0.67, -0.857, -0.613, -0.927 };
   int incX = -1;
   double Y[] = { -0.398, -0.934, -0.204, 0.183 };
   int incY = -1;
   double y_expected[] = { 0.745218, -0.60699, -0.37301, -0.983688 };
   cblas_zhpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhpmv(case 1127) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhpmv(case 1127) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   double alpha[2] = {1, 0};
   double beta[2] = {0, 1};
   int N = 2;
   double A[] = { 0.543, -0.737, 0.281, -0.053, -0.098, -0.826 };
   double X[] = { 0.67, -0.857, -0.613, -0.927 };
   int incX = -1;
   double Y[] = { -0.398, -0.934, -0.204, 0.183 };
   int incY = -1;
   double y_expected[] = { 0.646956, -0.542012, -0.282168, -0.912668 };
   cblas_zhpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhpmv(case 1128) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhpmv(case 1128) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   double alpha[2] = {1, 0};
   double beta[2] = {0, 1};
   int N = 2;
   double A[] = { 0.543, -0.737, 0.281, -0.053, -0.098, -0.826 };
   double X[] = { 0.67, -0.857, -0.613, -0.927 };
   int incX = -1;
   double Y[] = { -0.398, -0.934, -0.204, 0.183 };
   int incY = -1;
   double y_expected[] = { 0.646956, -0.542012, -0.282168, -0.912668 };
   cblas_zhpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhpmv(case 1129) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhpmv(case 1129) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   double alpha[2] = {1, 0};
   double beta[2] = {0, 1};
   int N = 2;
   double A[] = { 0.543, -0.737, 0.281, -0.053, -0.098, -0.826 };
   double X[] = { 0.67, -0.857, -0.613, -0.927 };
   int incX = -1;
   double Y[] = { -0.398, -0.934, -0.204, 0.183 };
   int incY = -1;
   double y_expected[] = { 0.745218, -0.60699, -0.37301, -0.983688 };
   cblas_zhpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhpmv(case 1130) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhpmv(case 1130) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   double alpha[2] = {1, 0};
   double beta[2] = {0, 1};
   int N = 2;
   double A[] = { 0.543, -0.737, 0.281, -0.053, -0.098, -0.826 };
   double X[] = { 0.67, -0.857, -0.613, -0.927 };
   int incX = -1;
   double Y[] = { -0.398, -0.934, -0.204, 0.183 };
   int incY = -1;
   double y_expected[] = { 0.745218, -0.60699, -0.37301, -0.983688 };
   cblas_zhpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhpmv(case 1131) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhpmv(case 1131) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   double alpha[2] = {1, 0};
   double beta[2] = {0, 1};
   int N = 2;
   double A[] = { 0.543, -0.737, 0.281, -0.053, -0.098, -0.826 };
   double X[] = { 0.67, -0.857, -0.613, -0.927 };
   int incX = -1;
   double Y[] = { -0.398, -0.934, -0.204, 0.183 };
   int incY = -1;
   double y_expected[] = { 0.646956, -0.542012, -0.282168, -0.912668 };
   cblas_zhpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhpmv(case 1132) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhpmv(case 1132) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   double alpha[2] = {1, 0};
   double beta[2] = {0, 1};
   int N = 2;
   double A[] = { 0.543, -0.737, 0.281, -0.053, -0.098, -0.826 };
   double X[] = { 0.67, -0.857, -0.613, -0.927 };
   int incX = -1;
   double Y[] = { -0.398, -0.934, -0.204, 0.183 };
   int incY = -1;
   double y_expected[] = { 0.646956, -0.542012, -0.282168, -0.912668 };
   cblas_zhpmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhpmv(case 1133) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhpmv(case 1133) imag");
     };
   };
  };


}
