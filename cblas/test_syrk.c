#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_syrk () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 1;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {1, 0};
   float A[] = { -0.392, -0.07 };
   int lda = 1;
   float C[] = { -0.044, 0.563 };
   int ldc = 1;
   float C_expected[] = { -0.0941172, 0.5614124 };
   cblas_csyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyrk(case 1312) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyrk(case 1312) imag");
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
   float beta[2] = {1, 0};
   float A[] = { -0.437, 0.787 };
   int lda = 1;
   float C[] = { 0.11, -0.826 };
   int ldc = 1;
   float C_expected[] = { 0.3073038, -0.6624886 };
   cblas_csyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyrk(case 1313) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyrk(case 1313) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 112;
   int N = 1;
   int K = 1;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {-1, 0};
   float A[] = { 0.25, 0.86 };
   int lda = 1;
   float C[] = { 0.33, 0.267 };
   int ldc = 1;
   float C_expected[] = { -0.16987, -0.46371 };
   cblas_csyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyrk(case 1314) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyrk(case 1314) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 112;
   int N = 1;
   int K = 1;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {-1, 0};
   float A[] = { 0.097, 0.65 };
   int lda = 1;
   float C[] = { -0.018, 0.424 };
   int ldc = 1;
   float C_expected[] = { 0.1293173, -0.5031391 };
   cblas_csyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyrk(case 1315) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyrk(case 1315) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 113;
   int N = 1;
   int K = 1;
   float alpha[2] = {0, 1};
   float beta[2] = {1, 0};
   float A[] = { -0.089, -0.847 };
   int lda = 1;
   float C[] = { -0.139, 0.509 };
   int ldc = 1;
   float C_expected[] = { -0.289766, -0.200488 };
   cblas_csyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyrk(case 1316) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyrk(case 1316) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 113;
   int N = 1;
   int K = 1;
   float alpha[2] = {0, 1};
   float beta[2] = {1, 0};
   float A[] = { -0.811, 0.032 };
   int lda = 1;
   float C[] = { -0.573, -0.663 };
   int ldc = 1;
   float C_expected[] = { -0.521096, -0.006303 };
   cblas_csyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyrk(case 1317) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyrk(case 1317) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 1;
   float alpha[2] = {0, 0.1};
   float beta[2] = {0, 1};
   float A[] = { -0.848, 0.179 };
   int lda = 1;
   float C[] = { -0.134, 0.608 };
   int ldc = 1;
   float C_expected[] = { -0.5776416, -0.0652937 };
   cblas_csyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyrk(case 1318) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyrk(case 1318) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 1;
   float alpha[2] = {0, 0.1};
   float beta[2] = {0, 1};
   float A[] = { 0.357, -0.576 };
   int lda = 1;
   float C[] = { 0.041, 0.099 };
   int ldc = 1;
   float C_expected[] = { -0.0578736, 0.0205673 };
   cblas_csyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyrk(case 1319) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyrk(case 1319) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 112;
   int N = 1;
   int K = 1;
   float alpha[2] = {0, 0};
   float beta[2] = {1, 0};
   float A[] = { 0.298, -0.61 };
   int lda = 1;
   float C[] = { 0.768, 0.203 };
   int ldc = 1;
   float C_expected[] = { 0.768, 0.203 };
   cblas_csyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyrk(case 1320) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyrk(case 1320) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 112;
   int N = 1;
   int K = 1;
   float alpha[2] = {0, 0};
   float beta[2] = {1, 0};
   float A[] = { -0.342, -0.083 };
   int lda = 1;
   float C[] = { -0.393, -0.984 };
   int ldc = 1;
   float C_expected[] = { -0.393, -0.984 };
   cblas_csyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyrk(case 1321) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyrk(case 1321) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 113;
   int N = 1;
   int K = 1;
   float alpha[2] = {0, 0};
   float beta[2] = {-0.3, 0.1};
   float A[] = { -0.037, 0.959 };
   int lda = 1;
   float C[] = { -0.079, -0.424 };
   int ldc = 1;
   float C_expected[] = { 0.0661, 0.1193 };
   cblas_csyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyrk(case 1322) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyrk(case 1322) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 113;
   int N = 1;
   int K = 1;
   float alpha[2] = {0, 0};
   float beta[2] = {-0.3, 0.1};
   float A[] = { -0.499, 0.034 };
   int lda = 1;
   float C[] = { 0.296, -0.599 };
   int ldc = 1;
   float C_expected[] = { -0.0289, 0.2093 };
   cblas_csyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyrk(case 1323) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyrk(case 1323) imag");
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
   double beta[2] = {0, 1};
   double A[] = { -0.166, -0.946 };
   int lda = 1;
   double C[] = { 0.004, 0.37 };
   int ldc = 1;
   double C_expected[] = { 0.49736, -0.310072 };
   cblas_zsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyrk(case 1324) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyrk(case 1324) imag");
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
   double beta[2] = {0, 1};
   double A[] = { -0.674, 0.715 };
   int lda = 1;
   double C[] = { -0.922, 0.012 };
   int ldc = 1;
   double C_expected[] = { 0.044949, 0.04182 };
   cblas_zsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyrk(case 1325) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyrk(case 1325) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 112;
   int N = 1;
   int K = 1;
   double alpha[2] = {-1, 0};
   double beta[2] = {1, 0};
   double A[] = { 0.857, -0.491 };
   int lda = 1;
   double C[] = { -0.692, 0.039 };
   int ldc = 1;
   double C_expected[] = { -1.185368, 0.880574 };
   cblas_zsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyrk(case 1326) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyrk(case 1326) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 112;
   int N = 1;
   int K = 1;
   double alpha[2] = {-1, 0};
   double beta[2] = {1, 0};
   double A[] = { -0.123, 0.981 };
   int lda = 1;
   double C[] = { -0.217, -0.429 };
   int ldc = 1;
   double C_expected[] = { 0.730232, -0.187674 };
   cblas_zsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyrk(case 1327) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyrk(case 1327) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 113;
   int N = 1;
   int K = 1;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {0, 1};
   double A[] = { -0.872, 0.108 };
   int lda = 1;
   double C[] = { 0.321, 0.98 };
   int ldc = 1;
   double C_expected[] = { -1.1857808, 0.4523776 };
   cblas_zsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyrk(case 1328) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyrk(case 1328) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 113;
   int N = 1;
   int K = 1;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {0, 1};
   double A[] = { 0.348, 0.447 };
   int lda = 1;
   double C[] = { 0.231, -0.973 };
   int ldc = 1;
   double C_expected[] = { 0.9655003, 0.1297959 };
   cblas_zsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyrk(case 1329) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyrk(case 1329) imag");
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
   double beta[2] = {0, 1};
   double A[] = { -0.744, -0.513 };
   int lda = 1;
   double C[] = { 0.282, -0.841 };
   int ldc = 1;
   double C_expected[] = { 0.550633, -0.481344 };
   cblas_zsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyrk(case 1330) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyrk(case 1330) imag");
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
   double beta[2] = {0, 1};
   double A[] = { 0.307, -0.949 };
   int lda = 1;
   double C[] = { 0.087, 0.196 };
   int ldc = 1;
   double C_expected[] = { 0.610352, 0.669686 };
   cblas_zsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyrk(case 1331) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyrk(case 1331) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 112;
   int N = 1;
   int K = 1;
   double alpha[2] = {0, 1};
   double beta[2] = {-0.3, 0.1};
   double A[] = { 0.293, 0.566 };
   int lda = 1;
   double C[] = { -0.373, 0.92 };
   int ldc = 1;
   double C_expected[] = { -0.311776, -0.547807 };
   cblas_zsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyrk(case 1332) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyrk(case 1332) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 112;
   int N = 1;
   int K = 1;
   double alpha[2] = {0, 1};
   double beta[2] = {-0.3, 0.1};
   double A[] = { 0.404, 0.861 };
   int lda = 1;
   double C[] = { 0.627, -0.434 };
   int ldc = 1;
   double C_expected[] = { -0.840388, -0.385205 };
   cblas_zsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyrk(case 1333) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyrk(case 1333) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 113;
   int N = 1;
   int K = 1;
   double alpha[2] = {0, 1};
   double beta[2] = {0, 0};
   double A[] = { -0.143, -0.036 };
   int lda = 1;
   double C[] = { 0.327, 0.125 };
   int ldc = 1;
   double C_expected[] = { -0.010296, 0.019153 };
   cblas_zsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyrk(case 1334) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyrk(case 1334) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 113;
   int N = 1;
   int K = 1;
   double alpha[2] = {0, 1};
   double beta[2] = {0, 0};
   double A[] = { 0.22, -0.858 };
   int lda = 1;
   double C[] = { -0.776, 0.667 };
   int ldc = 1;
   double C_expected[] = { 0.37752, -0.687764 };
   cblas_zsyrk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyrk(case 1335) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyrk(case 1335) imag");
     };
   };
  };


}
