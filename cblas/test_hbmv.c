#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_hbmv () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int uplo = 121;
   float alpha[2] = {0, 1};
   float beta[2] = {-0.3, 0.1};
   int N = 3;
   int k = 1;
   int lda = 3;
   float A[] = { 0.937, -0.035, 0.339, 0.847, 0.022, 0.153, -0.785, 0.193, -0.731, -0.166, -0.243, -0.319, 0.173, -0.24, 0.079, -0.058, 0.124, 0.445 };
   float X[] = { -0.093, -0.103, -0.537, -0.151, 0.094, 0.954 };
   int incX = -1;
   float Y[] = { 0.029, -0.391, -0.256, 0.031, -0.478, 0.098 };
   int incY = -1;
   float y_expected[] = { 0.02698, 0.521724, -0.379354, 1.27743, -0.25427, -0.043268 };
   cblas_chbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chbmv(case 1086) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chbmv(case 1086) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   float alpha[2] = {0, 1};
   float beta[2] = {-0.3, 0.1};
   int N = 3;
   int k = 1;
   int lda = 3;
   float A[] = { 0.937, -0.035, 0.339, 0.847, 0.022, 0.153, -0.785, 0.193, -0.731, -0.166, -0.243, -0.319, 0.173, -0.24, 0.079, -0.058, 0.124, 0.445 };
   float X[] = { -0.093, -0.103, -0.537, -0.151, 0.094, 0.954 };
   int incX = -1;
   float Y[] = { 0.029, -0.391, -0.256, 0.031, -0.478, 0.098 };
   int incY = -1;
   float y_expected[] = { 0.02698, 0.521724, -0.379354, 1.27743, -0.25427, -0.043268 };
   cblas_chbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chbmv(case 1087) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chbmv(case 1087) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   float alpha[2] = {0, 1};
   float beta[2] = {-0.3, 0.1};
   int N = 3;
   int k = 1;
   int lda = 3;
   float A[] = { 0.937, -0.035, 0.339, 0.847, 0.022, 0.153, -0.785, 0.193, -0.731, -0.166, -0.243, -0.319, 0.173, -0.24, 0.079, -0.058, 0.124, 0.445 };
   float X[] = { -0.093, -0.103, -0.537, -0.151, 0.094, 0.954 };
   int incX = -1;
   float Y[] = { 0.029, -0.391, -0.256, 0.031, -0.478, 0.098 };
   int incY = -1;
   float y_expected[] = { -0.06422, -0.016288, 0.734206, 0.108366, -0.411982, 0.347068 };
   cblas_chbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chbmv(case 1088) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chbmv(case 1088) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   float alpha[2] = {0, 1};
   float beta[2] = {-0.3, 0.1};
   int N = 3;
   int k = 1;
   int lda = 3;
   float A[] = { 0.937, -0.035, 0.339, 0.847, 0.022, 0.153, -0.785, 0.193, -0.731, -0.166, -0.243, -0.319, 0.173, -0.24, 0.079, -0.058, 0.124, 0.445 };
   float X[] = { -0.093, -0.103, -0.537, -0.151, 0.094, 0.954 };
   int incX = -1;
   float Y[] = { 0.029, -0.391, -0.256, 0.031, -0.478, 0.098 };
   int incY = -1;
   float y_expected[] = { -0.06422, -0.016288, 0.734206, 0.108366, -0.411982, 0.347068 };
   cblas_chbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chbmv(case 1089) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chbmv(case 1089) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   float alpha[2] = {0, 1};
   float beta[2] = {-0.3, 0.1};
   int N = 3;
   int k = 1;
   int lda = 3;
   float A[] = { 0.937, -0.035, 0.339, 0.847, 0.022, 0.153, -0.785, 0.193, -0.731, -0.166, -0.243, -0.319, 0.173, -0.24, 0.079, -0.058, 0.124, 0.445 };
   float X[] = { -0.093, -0.103, -0.537, -0.151, 0.094, 0.954 };
   int incX = -1;
   float Y[] = { 0.029, -0.391, -0.256, 0.031, -0.478, 0.098 };
   int incY = -1;
   float y_expected[] = { 0.19354, 0.056192, 0.72585, 0.42717, -0.2047, 0.405354 };
   cblas_chbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chbmv(case 1090) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chbmv(case 1090) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   float alpha[2] = {0, 1};
   float beta[2] = {-0.3, 0.1};
   int N = 3;
   int k = 1;
   int lda = 3;
   float A[] = { 0.937, -0.035, 0.339, 0.847, 0.022, 0.153, -0.785, 0.193, -0.731, -0.166, -0.243, -0.319, 0.173, -0.24, 0.079, -0.058, 0.124, 0.445 };
   float X[] = { -0.093, -0.103, -0.537, -0.151, 0.094, 0.954 };
   int incX = -1;
   float Y[] = { 0.029, -0.391, -0.256, 0.031, -0.478, 0.098 };
   int incY = -1;
   float y_expected[] = { 0.19354, 0.056192, 0.72585, 0.42717, -0.2047, 0.405354 };
   cblas_chbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chbmv(case 1091) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chbmv(case 1091) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   float alpha[2] = {0, 1};
   float beta[2] = {-0.3, 0.1};
   int N = 3;
   int k = 1;
   int lda = 3;
   float A[] = { 0.937, -0.035, 0.339, 0.847, 0.022, 0.153, -0.785, 0.193, -0.731, -0.166, -0.243, -0.319, 0.173, -0.24, 0.079, -0.058, 0.124, 0.445 };
   float X[] = { -0.093, -0.103, -0.537, -0.151, 0.094, 0.954 };
   int incX = -1;
   float Y[] = { 0.029, -0.391, -0.256, 0.031, -0.478, 0.098 };
   int incY = -1;
   float y_expected[] = { -0.151304, 0.471592, -0.507714, -0.304446, -1.16395, -0.299062 };
   cblas_chbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chbmv(case 1092) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chbmv(case 1092) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   float alpha[2] = {0, 1};
   float beta[2] = {-0.3, 0.1};
   int N = 3;
   int k = 1;
   int lda = 3;
   float A[] = { 0.937, -0.035, 0.339, 0.847, 0.022, 0.153, -0.785, 0.193, -0.731, -0.166, -0.243, -0.319, 0.173, -0.24, 0.079, -0.058, 0.124, 0.445 };
   float X[] = { -0.093, -0.103, -0.537, -0.151, 0.094, 0.954 };
   int incX = -1;
   float Y[] = { 0.029, -0.391, -0.256, 0.031, -0.478, 0.098 };
   int incY = -1;
   float y_expected[] = { -0.151304, 0.471592, -0.507714, -0.304446, -1.16395, -0.299062 };
   cblas_chbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], flteps, "chbmv(case 1093) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], flteps, "chbmv(case 1093) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   double alpha[2] = {1, 0};
   double beta[2] = {1, 0};
   int N = 3;
   int k = 1;
   int lda = 3;
   double A[] = { 0.662, 0.24, -0.311, -0.345, -0.782, 0.904, -0.842, 0.065, -0.168, -0.855, -0.692, 0.113, 0.009, -0.707, -0.981, 0.019, -0.687, 0.861 };
   double X[] = { 0.873, -0.509, 0.398, 0.471, 0.214, 0.878 };
   int incX = -1;
   double Y[] = { -0.441, -0.781, 0.979, -0.911, 0.879, 0.807 };
   int incY = -1;
   double y_expected[] = { -0.902712, -0.524419, -0.307439, -2.167713, 1.059385, 1.104445 };
   cblas_zhbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhbmv(case 1094) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhbmv(case 1094) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   double alpha[2] = {1, 0};
   double beta[2] = {1, 0};
   int N = 3;
   int k = 1;
   int lda = 3;
   double A[] = { 0.662, 0.24, -0.311, -0.345, -0.782, 0.904, -0.842, 0.065, -0.168, -0.855, -0.692, 0.113, 0.009, -0.707, -0.981, 0.019, -0.687, 0.861 };
   double X[] = { 0.873, -0.509, 0.398, 0.471, 0.214, 0.878 };
   int incX = -1;
   double Y[] = { -0.441, -0.781, 0.979, -0.911, 0.879, 0.807 };
   int incY = -1;
   double y_expected[] = { -0.902712, -0.524419, -0.307439, -2.167713, 1.059385, 1.104445 };
   cblas_zhbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhbmv(case 1095) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhbmv(case 1095) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   double alpha[2] = {1, 0};
   double beta[2] = {1, 0};
   int N = 3;
   int k = 1;
   int lda = 3;
   double A[] = { 0.662, 0.24, -0.311, -0.345, -0.782, 0.904, -0.842, 0.065, -0.168, -0.855, -0.692, 0.113, 0.009, -0.707, -0.981, 0.019, -0.687, 0.861 };
   double X[] = { 0.873, -0.509, 0.398, 0.471, 0.214, 0.878 };
   int incX = -1;
   double Y[] = { -0.441, -0.781, 0.979, -0.911, 0.879, 0.807 };
   int incY = -1;
   double y_expected[] = { -0.960834, -0.558818, 1.042598, -1.102864, 0.507945, 0.11149 };
   cblas_zhbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhbmv(case 1096) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhbmv(case 1096) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   double alpha[2] = {1, 0};
   double beta[2] = {1, 0};
   int N = 3;
   int k = 1;
   int lda = 3;
   double A[] = { 0.662, 0.24, -0.311, -0.345, -0.782, 0.904, -0.842, 0.065, -0.168, -0.855, -0.692, 0.113, 0.009, -0.707, -0.981, 0.019, -0.687, 0.861 };
   double X[] = { 0.873, -0.509, 0.398, 0.471, 0.214, 0.878 };
   int incX = -1;
   double Y[] = { -0.441, -0.781, 0.979, -0.911, 0.879, 0.807 };
   int incY = -1;
   double y_expected[] = { -0.960834, -0.558818, 1.042598, -1.102864, 0.507945, 0.11149 };
   cblas_zhbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhbmv(case 1097) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhbmv(case 1097) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   double alpha[2] = {1, 0};
   double beta[2] = {1, 0};
   int N = 3;
   int k = 1;
   int lda = 3;
   double A[] = { 0.662, 0.24, -0.311, -0.345, -0.782, 0.904, -0.842, 0.065, -0.168, -0.855, -0.692, 0.113, 0.009, -0.707, -0.981, 0.019, -0.687, 0.861 };
   double X[] = { 0.873, -0.509, 0.398, 0.471, 0.214, 0.878 };
   int incX = -1;
   double Y[] = { -0.441, -0.781, 0.979, -0.911, 0.879, 0.807 };
   int incY = -1;
   double y_expected[] = { -1.626828, 0.003954, 0.437012, -2.365106, 0.446715, 0.16323 };
   cblas_zhbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhbmv(case 1098) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhbmv(case 1098) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   double alpha[2] = {1, 0};
   double beta[2] = {1, 0};
   int N = 3;
   int k = 1;
   int lda = 3;
   double A[] = { 0.662, 0.24, -0.311, -0.345, -0.782, 0.904, -0.842, 0.065, -0.168, -0.855, -0.692, 0.113, 0.009, -0.707, -0.981, 0.019, -0.687, 0.861 };
   double X[] = { 0.873, -0.509, 0.398, 0.471, 0.214, 0.878 };
   int incX = -1;
   double Y[] = { -0.441, -0.781, 0.979, -0.911, 0.879, 0.807 };
   int incY = -1;
   double y_expected[] = { -1.626828, 0.003954, 0.437012, -2.365106, 0.446715, 0.16323 };
   cblas_zhbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhbmv(case 1099) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhbmv(case 1099) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   double alpha[2] = {1, 0};
   double beta[2] = {1, 0};
   int N = 3;
   int k = 1;
   int lda = 3;
   double A[] = { 0.662, 0.24, -0.311, -0.345, -0.782, 0.904, -0.842, 0.065, -0.168, -0.855, -0.692, 0.113, 0.009, -0.707, -0.981, 0.019, -0.687, 0.861 };
   double X[] = { 0.873, -0.509, 0.398, 0.471, 0.214, 0.878 };
   int incX = -1;
   double Y[] = { -0.441, -0.781, 0.979, -0.911, 0.879, 0.807 };
   int incY = -1;
   double y_expected[] = { -0.097302, -1.204999, 1.168771, -0.822543, 0.734395, 1.379065 };
   cblas_zhbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhbmv(case 1100) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhbmv(case 1100) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   double alpha[2] = {1, 0};
   double beta[2] = {1, 0};
   int N = 3;
   int k = 1;
   int lda = 3;
   double A[] = { 0.662, 0.24, -0.311, -0.345, -0.782, 0.904, -0.842, 0.065, -0.168, -0.855, -0.692, 0.113, 0.009, -0.707, -0.981, 0.019, -0.687, 0.861 };
   double X[] = { 0.873, -0.509, 0.398, 0.471, 0.214, 0.878 };
   int incX = -1;
   double Y[] = { -0.441, -0.781, 0.979, -0.911, 0.879, 0.807 };
   int incY = -1;
   double y_expected[] = { -0.097302, -1.204999, 1.168771, -0.822543, 0.734395, 1.379065 };
   cblas_zhbmv(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(Y[2*i], y_expected[2*i], dbleps, "zhbmv(case 1101) real");
       gsl_test_rel(Y[2*i+1], y_expected[2*i+1], dbleps, "zhbmv(case 1101) imag");
     };
   };
  };


}
