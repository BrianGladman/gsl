#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_hemm () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int M = 1;
   int N = 1;
   float alpha[2] = {1, 0};
   float beta[2] = {0, 1};
   float A[] = { 0.478, 0.503 };
   int lda = 1;
   float B[] = { 0.313, -0.565 };
   int ldb = 1;
   float C[] = { -0.174, 0.878 };
   int ldc = 1;
   float C_expected[] = { -0.728386, -0.44407 };
   cblas_chemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "chemm(case 1296) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "chemm(case 1296) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int M = 1;
   int N = 1;
   float alpha[2] = {1, 0};
   float beta[2] = {0, 1};
   float A[] = { 0.653, -0.181 };
   int lda = 1;
   float B[] = { -0.071, -0.038 };
   int ldb = 1;
   float C[] = { -0.109, -0.359 };
   int ldc = 1;
   float C_expected[] = { 0.312637, -0.133814 };
   cblas_chemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "chemm(case 1297) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "chemm(case 1297) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int M = 1;
   int N = 1;
   float alpha[2] = {1, 0};
   float beta[2] = {0, 1};
   float A[] = { 0.388, -0.44 };
   int lda = 1;
   float B[] = { 0.995, 0.348 };
   int ldb = 1;
   float C[] = { -0.449, -0.219 };
   int ldc = 1;
   float C_expected[] = { 0.60506, -0.313976 };
   cblas_chemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "chemm(case 1298) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "chemm(case 1298) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int M = 1;
   int N = 1;
   float alpha[2] = {1, 0};
   float beta[2] = {0, 1};
   float A[] = { -0.21, 0.473 };
   int lda = 1;
   float B[] = { 0.748, 0.979 };
   int ldb = 1;
   float C[] = { 0.793, 0.338 };
   int ldc = 1;
   float C_expected[] = { -0.49508, 0.58741 };
   cblas_chemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "chemm(case 1299) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "chemm(case 1299) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int M = 1;
   int N = 1;
   float alpha[2] = {-1, 0};
   float beta[2] = {0, 0.1};
   float A[] = { 0.6, 0.041 };
   int lda = 1;
   float B[] = { 0.896, -0.447 };
   int ldb = 1;
   float C[] = { -0.95, 0.343 };
   int ldc = 1;
   float C_expected[] = { -0.5719, 0.1732 };
   cblas_chemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "chemm(case 1300) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "chemm(case 1300) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int M = 1;
   int N = 1;
   float alpha[2] = {-1, 0};
   float beta[2] = {0, 0.1};
   float A[] = { 0.765, 0.18 };
   int lda = 1;
   float B[] = { 0.977, -0.955 };
   int ldb = 1;
   float C[] = { 0.397, 0.683 };
   int ldc = 1;
   float C_expected[] = { -0.815705, 0.770275 };
   cblas_chemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "chemm(case 1301) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "chemm(case 1301) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float beta[2] = {0, 1};
   float A[] = { 0.3, -0.874 };
   int lda = 1;
   float B[] = { 0.372, -0.745 };
   int ldb = 1;
   float C[] = { 0.076, -0.16 };
   int ldc = 1;
   float C_expected[] = { 0.18235, 0.08716 };
   cblas_chemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "chemm(case 1302) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "chemm(case 1302) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 0.1};
   float beta[2] = {0, 1};
   float A[] = { -0.681, 0.209 };
   int lda = 1;
   float B[] = { 0.436, -0.369 };
   int ldb = 1;
   float C[] = { -0.085, -0.303 };
   int ldc = 1;
   float C_expected[] = { 0.2778711, -0.1146916 };
   cblas_chemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "chemm(case 1303) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "chemm(case 1303) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int M = 1;
   int N = 1;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {0, 0};
   double A[] = { 0.734, -0.305 };
   int lda = 1;
   double B[] = { 0.61, -0.831 };
   int ldb = 1;
   double C[] = { 0.86, -0.233 };
   int ldc = 1;
   double C_expected[] = { -0.0733266, 0.2277602 };
   cblas_zhemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zhemm(case 1304) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zhemm(case 1304) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int M = 1;
   int N = 1;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {0, 0};
   double A[] = { -0.389, -0.205 };
   int lda = 1;
   double B[] = { -0.121, 0.323 };
   int ldb = 1;
   double C[] = { 0.022, 0.795 };
   int ldc = 1;
   double C_expected[] = { -0.001556, 0.042401 };
   cblas_zhemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zhemm(case 1305) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zhemm(case 1305) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int M = 1;
   int N = 1;
   double alpha[2] = {0, 1};
   double beta[2] = {-0.3, 0.1};
   double A[] = { 0.635, 0.116 };
   int lda = 1;
   double B[] = { 0.619, -0.443 };
   int ldb = 1;
   double C[] = { 0.742, 0.144 };
   int ldc = 1;
   double C_expected[] = { 0.044305, 0.424065 };
   cblas_zhemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zhemm(case 1306) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zhemm(case 1306) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int M = 1;
   int N = 1;
   double alpha[2] = {0, 1};
   double beta[2] = {-0.3, 0.1};
   double A[] = { -0.698, 0.783 };
   int lda = 1;
   double B[] = { -0.343, -0.603 };
   int ldb = 1;
   double C[] = { 0.957, -0.633 };
   int ldc = 1;
   double C_expected[] = { -0.644694, 0.525014 };
   cblas_zhemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zhemm(case 1307) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zhemm(case 1307) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int M = 1;
   int N = 1;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {0, 1};
   double A[] = { -0.199, 0.303 };
   int lda = 1;
   double B[] = { -0.705, -0.013 };
   int ldb = 1;
   double C[] = { 0.588, 0.252 };
   int ldc = 1;
   double C_expected[] = { -0.2943472, 0.6012534 };
   cblas_zhemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zhemm(case 1308) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zhemm(case 1308) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int M = 1;
   int N = 1;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {0, 1};
   double A[] = { 0.756, -0.177 };
   int lda = 1;
   double B[] = { -0.079, 0.58 };
   int ldb = 1;
   double C[] = { -0.678, 0.547 };
   int ldc = 1;
   double C_expected[] = { -0.5729308, -0.8155164 };
   cblas_zhemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zhemm(case 1309) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zhemm(case 1309) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int M = 1;
   int N = 1;
   double alpha[2] = {-1, 0};
   double beta[2] = {0, 0.1};
   double A[] = { -0.961, -0.503 };
   int lda = 1;
   double B[] = { 0.81, -0.027 };
   int ldb = 1;
   double C[] = { -0.975, 0.813 };
   int ldc = 1;
   double C_expected[] = { 0.69711, -0.123447 };
   cblas_zhemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zhemm(case 1310) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zhemm(case 1310) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int M = 1;
   int N = 1;
   double alpha[2] = {-1, 0};
   double beta[2] = {0, 0.1};
   double A[] = { 0.798, -0.667 };
   int lda = 1;
   double B[] = { -0.962, 0.226 };
   int ldb = 1;
   double C[] = { 0.212, 0.446 };
   int ldc = 1;
   double C_expected[] = { 0.723076, -0.159148 };
   cblas_zhemm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zhemm(case 1311) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zhemm(case 1311) imag");
     };
   };
  };


}
