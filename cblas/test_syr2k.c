#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_syr2k () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 1;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {0, 0};
   float A[] = { -0.831, -0.163 };
   int lda = 1;
   float B[] = { 0.489, 0.154 };
   int ldb = 1;
   float C[] = { 0.493, -0.175 };
   int ldc = 1;
   float C_expected[] = { 0.2702904, 0.0483572 };
   cblas_csyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyr2k(case 1360) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyr2k(case 1360) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 1;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {0, 0};
   float A[] = { 0.938, 0.342 };
   int lda = 1;
   float B[] = { 0.74, 0.216 };
   int ldb = 1;
   float C[] = { 0.769, 0.572 };
   int ldc = 1;
   float C_expected[] = { -0.4632864, -0.1493632 };
   cblas_csyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyr2k(case 1361) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyr2k(case 1361) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 1;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {0, 0};
   float A[] = { 0.282, -0.77 };
   int lda = 1;
   float B[] = { -0.821, 0.954 };
   int ldb = 1;
   float C[] = { -0.566, -0.845 };
   int ldc = 1;
   float C_expected[] = { -0.4820744, -0.4401072 };
   cblas_csyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyr2k(case 1362) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyr2k(case 1362) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 1;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {0, 0};
   float A[] = { 0.016, -0.334 };
   int lda = 1;
   float B[] = { 0.424, -0.334 };
   int ldb = 1;
   float C[] = { 0.532, 0.802 };
   int ldc = 1;
   float C_expected[] = { 0.0922552, 0.0672216 };
   cblas_csyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyr2k(case 1363) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyr2k(case 1363) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 112;
   int N = 1;
   int K = 1;
   float alpha[2] = {0, 1};
   float beta[2] = {1, 0};
   float A[] = { 0.358, 0.783 };
   int lda = 1;
   float B[] = { 0.159, -0.13 };
   int ldb = 1;
   float C[] = { -0.135, 0.455 };
   int ldc = 1;
   float C_expected[] = { -0.290914, 0.772424 };
   cblas_csyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyr2k(case 1364) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyr2k(case 1364) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 112;
   int N = 1;
   int K = 1;
   float alpha[2] = {0, 1};
   float beta[2] = {1, 0};
   float A[] = { 0.526, -0.267 };
   int lda = 1;
   float B[] = { 0.397, 0.772 };
   int ldb = 1;
   float C[] = { 0.854, 0.851 };
   int ldc = 1;
   float C_expected[] = { 0.253854, 1.680892 };
   cblas_csyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyr2k(case 1365) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyr2k(case 1365) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 112;
   int N = 1;
   int K = 1;
   float alpha[2] = {0, 1};
   float beta[2] = {1, 0};
   float A[] = { -0.839, 0.941 };
   int lda = 1;
   float B[] = { -0.422, 0.891 };
   int ldb = 1;
   float C[] = { 0.997, -0.173 };
   int ldc = 1;
   float C_expected[] = { 3.286302, -1.141746 };
   cblas_csyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyr2k(case 1366) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyr2k(case 1366) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 112;
   int N = 1;
   int K = 1;
   float alpha[2] = {0, 1};
   float beta[2] = {1, 0};
   float A[] = { 0.498, -0.925 };
   int lda = 1;
   float B[] = { 0.199, 0.237 };
   int ldb = 1;
   float C[] = { 0.161, -0.703 };
   int ldc = 1;
   float C_expected[] = { 0.293098, -0.066346 };
   cblas_csyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyr2k(case 1367) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyr2k(case 1367) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 1;
   double alpha[2] = {-1, 0};
   double beta[2] = {-1, 0};
   double A[] = { 0.323, -0.957 };
   int lda = 1;
   double B[] = { -0.303, -0.873 };
   int ldb = 1;
   double C[] = { 0.842, -0.734 };
   int ldc = 1;
   double C_expected[] = { 1.02466, 0.718016 };
   cblas_zsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyr2k(case 1368) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyr2k(case 1368) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 1;
   double alpha[2] = {-1, 0};
   double beta[2] = {-1, 0};
   double A[] = { 0.583, 0.522 };
   int lda = 1;
   double B[] = { -0.83, 0.922 };
   int ldb = 1;
   double C[] = { -0.871, -0.819 };
   int ldc = 1;
   double C_expected[] = { 2.801348, 0.610468 };
   cblas_zsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyr2k(case 1369) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyr2k(case 1369) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 1;
   double alpha[2] = {-1, 0};
   double beta[2] = {-1, 0};
   double A[] = { 0.072, 0.345 };
   int lda = 1;
   double B[] = { 0.944, -0.39 };
   int ldb = 1;
   double C[] = { -0.228, -0.003 };
   int ldc = 1;
   double C_expected[] = { -0.177036, -0.5922 };
   cblas_zsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyr2k(case 1370) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyr2k(case 1370) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 1;
   double alpha[2] = {-1, 0};
   double beta[2] = {-1, 0};
   double A[] = { -0.693, -0.453 };
   int lda = 1;
   double B[] = { -0.434, -0.293 };
   int ldb = 1;
   double C[] = { -0.577, 0.656 };
   int ldc = 1;
   double C_expected[] = { 0.240934, -1.455302 };
   cblas_zsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyr2k(case 1371) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyr2k(case 1371) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 112;
   int N = 1;
   int K = 1;
   double alpha[2] = {0, 0.1};
   double beta[2] = {0, 0.1};
   double A[] = { -0.494, 0.304 };
   int lda = 1;
   double B[] = { 0.147, 0.134 };
   int ldb = 1;
   double C[] = { -0.838, 0.622 };
   int ldc = 1;
   double C_expected[] = { -0.0578984, -0.1064708 };
   cblas_zsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyr2k(case 1372) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyr2k(case 1372) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 112;
   int N = 1;
   int K = 1;
   double alpha[2] = {0, 0.1};
   double beta[2] = {0, 0.1};
   double A[] = { -0.7, 0.541 };
   int lda = 1;
   double B[] = { -0.794, -0.256 };
   int ldb = 1;
   double C[] = { 0.169, 0.734 };
   int ldc = 1;
   double C_expected[] = { -0.0233292, 0.1557592 };
   cblas_zsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyr2k(case 1373) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyr2k(case 1373) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 112;
   int N = 1;
   int K = 1;
   double alpha[2] = {0, 0.1};
   double beta[2] = {0, 0.1};
   double A[] = { -0.937, 0.635 };
   int lda = 1;
   double B[] = { 0.596, -0.51 };
   int ldb = 1;
   double C[] = { 0.844, 0.999 };
   int ldc = 1;
   double C_expected[] = { -0.271166, 0.0374796 };
   cblas_zsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyr2k(case 1374) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyr2k(case 1374) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 112;
   int N = 1;
   int K = 1;
   double alpha[2] = {0, 0.1};
   double beta[2] = {0, 0.1};
   double A[] = { 0.049, 0.133 };
   int lda = 1;
   double B[] = { -0.918, -0.147 };
   int ldb = 1;
   double C[] = { -0.688, -0.265 };
   int ldc = 1;
   double C_expected[] = { 0.0523594, -0.0738862 };
   cblas_zsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyr2k(case 1375) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyr2k(case 1375) imag");
     };
   };
  };


}
