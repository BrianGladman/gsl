#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_herk () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 1;
   float alpha = 0;
   float beta = 0;
   float A[] = { 0.23, 0.079 };
   int lda = 1;
   float C[] = { 0.606, -0.865 };
   int ldc = 1;
   float C_expected[] = { 0.000000000000e+00, 0.000000000000e+00 };
   cblas_cherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cherk(case 1336) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cherk(case 1336) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 1;
   float alpha = 0;
   float beta = 0;
   float A[] = { -0.126, -0.804 };
   int lda = 1;
   float C[] = { 0.028, 0.582 };
   int ldc = 1;
   float C_expected[] = { 0.000000000000e+00, 0.000000000000e+00 };
   cblas_cherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cherk(case 1337) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cherk(case 1337) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 112;
   int N = 1;
   int K = 1;
   float alpha = 1;
   float beta = -1;
   float A[] = { 0.451, 0.144 };
   int lda = 1;
   float C[] = { -0.381, 0.53 };
   int ldc = 1;
   float C_expected[] = { 0.605137, 0 };
   cblas_cherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cherk(case 1338) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cherk(case 1338) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 112;
   int N = 1;
   int K = 1;
   float alpha = 1;
   float beta = -1;
   float A[] = { -0.698, -0.42 };
   int lda = 1;
   float C[] = { 0.032, 0.355 };
   int ldc = 1;
   float C_expected[] = { 0.631604, 0 };
   cblas_cherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cherk(case 1339) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cherk(case 1339) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 113;
   int N = 1;
   int K = 1;
   float alpha = 0;
   float beta = 0;
   float A[] = { 0.173, -0.38 };
   int lda = 1;
   float C[] = { -0.019, 0.588 };
   int ldc = 1;
   float C_expected[] = { 0.000000000000e+00, 0.000000000000e+00 };
   cblas_cherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cherk(case 1340) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cherk(case 1340) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 113;
   int N = 1;
   int K = 1;
   float alpha = 0;
   float beta = 0;
   float A[] = { -0.133, -0.072 };
   int lda = 1;
   float C[] = { 0.814, -0.937 };
   int ldc = 1;
   float C_expected[] = { 0.000000000000e+00, 0.000000000000e+00 };
   cblas_cherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cherk(case 1341) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cherk(case 1341) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 1;
   float alpha = -0.3;
   float beta = -0.3;
   float A[] = { 0.832, 0.875 };
   int lda = 1;
   float C[] = { 0.052, 0.751 };
   int ldc = 1;
   float C_expected[] = { -0.4529547, 0 };
   cblas_cherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cherk(case 1342) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cherk(case 1342) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 1;
   float alpha = -0.3;
   float beta = -0.3;
   float A[] = { 0.344, -0.668 };
   int lda = 1;
   float C[] = { -0.912, -0.153 };
   int ldc = 1;
   float C_expected[] = { 0.104232, 0 };
   cblas_cherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cherk(case 1343) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cherk(case 1343) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 112;
   int N = 1;
   int K = 1;
   float alpha = 1;
   float beta = 0;
   float A[] = { -0.082, -0.061 };
   int lda = 1;
   float C[] = { -0.515, 0.611 };
   int ldc = 1;
   float C_expected[] = { 0.010445, 0 };
   cblas_cherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cherk(case 1344) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cherk(case 1344) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 112;
   int N = 1;
   int K = 1;
   float alpha = 1;
   float beta = 0;
   float A[] = { 0.712, -0.075 };
   int lda = 1;
   float C[] = { -0.809, -0.652 };
   int ldc = 1;
   float C_expected[] = { 0.512569, 0 };
   cblas_cherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cherk(case 1345) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cherk(case 1345) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 113;
   int N = 1;
   int K = 1;
   float alpha = 0;
   float beta = -0.3;
   float A[] = { -0.627, 0.757 };
   int lda = 1;
   float C[] = { -0.38, -0.235 };
   int ldc = 1;
   float C_expected[] = { 0.114, 0 };
   cblas_cherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cherk(case 1346) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cherk(case 1346) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 113;
   int N = 1;
   int K = 1;
   float alpha = 0;
   float beta = -0.3;
   float A[] = { 0.757, -0.579 };
   int lda = 1;
   float C[] = { -0.914, -0.426 };
   int ldc = 1;
   float C_expected[] = { 0.2742, 0 };
   cblas_cherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cherk(case 1347) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cherk(case 1347) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 1;
   double alpha = 1;
   double beta = -0.3;
   double A[] = { -0.291, 0.845 };
   int lda = 1;
   double C[] = { 0.035, 0.621 };
   int ldc = 1;
   double C_expected[] = { 0.788206, 0 };
   cblas_zherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zherk(case 1348) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zherk(case 1348) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 1;
   double alpha = 1;
   double beta = -0.3;
   double A[] = { -0.891, 0.657 };
   int lda = 1;
   double C[] = { -0.406, -0.65 };
   int ldc = 1;
   double C_expected[] = { 1.34733, 0 };
   cblas_zherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zherk(case 1349) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zherk(case 1349) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 112;
   int N = 1;
   int K = 1;
   double alpha = 1;
   double beta = -0.3;
   double A[] = { 0.164, 0.257 };
   int lda = 1;
   double C[] = { -0.832, -0.151 };
   int ldc = 1;
   double C_expected[] = { 0.342545, 0 };
   cblas_zherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zherk(case 1350) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zherk(case 1350) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 112;
   int N = 1;
   int K = 1;
   double alpha = 1;
   double beta = -0.3;
   double A[] = { -0.009, -0.342 };
   int lda = 1;
   double C[] = { 0.326, 0.319 };
   int ldc = 1;
   double C_expected[] = { 0.019245, 0 };
   cblas_zherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zherk(case 1351) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zherk(case 1351) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 113;
   int N = 1;
   int K = 1;
   double alpha = 0.1;
   double beta = 0.1;
   double A[] = { -0.277, -0.369 };
   int lda = 1;
   double C[] = { 0.238, -0.499 };
   int ldc = 1;
   double C_expected[] = { 0.045089, 0 };
   cblas_zherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zherk(case 1352) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zherk(case 1352) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 113;
   int N = 1;
   int K = 1;
   double alpha = 0.1;
   double beta = 0.1;
   double A[] = { -0.905, 0.182 };
   int lda = 1;
   double C[] = { 0.638, 0.224 };
   int ldc = 1;
   double C_expected[] = { 0.1490149, 0 };
   cblas_zherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zherk(case 1353) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zherk(case 1353) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 1;
   double alpha = 0.1;
   double beta = -1;
   double A[] = { 0.362, 0.035 };
   int lda = 1;
   double C[] = { -0.855, 0.136 };
   int ldc = 1;
   double C_expected[] = { 0.8682269, 0 };
   cblas_zherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zherk(case 1354) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zherk(case 1354) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 1;
   double alpha = 0.1;
   double beta = -1;
   double A[] = { -0.736, -0.686 };
   int lda = 1;
   double C[] = { 0.133, -0.278 };
   int ldc = 1;
   double C_expected[] = { -0.0317708, 0 };
   cblas_zherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zherk(case 1355) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zherk(case 1355) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 112;
   int N = 1;
   int K = 1;
   double alpha = 0;
   double beta = 1;
   double A[] = { -0.064, 0.455 };
   int lda = 1;
   double C[] = { -0.764, -0.257 };
   int ldc = 1;
   double C_expected[] = { -0.764, -0.257 };
   cblas_zherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zherk(case 1356) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zherk(case 1356) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 112;
   int N = 1;
   int K = 1;
   double alpha = 0;
   double beta = 1;
   double A[] = { 0.819, 0.175 };
   int lda = 1;
   double C[] = { -0.285, 0.852 };
   int ldc = 1;
   double C_expected[] = { -0.285, 0.852 };
   cblas_zherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zherk(case 1357) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zherk(case 1357) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 113;
   int N = 1;
   int K = 1;
   double alpha = 0;
   double beta = 0;
   double A[] = { 0.862, 0.823 };
   int lda = 1;
   double C[] = { -0.222, 0.373 };
   int ldc = 1;
   double C_expected[] = { 0.000000000000e+00, 0.000000000000e+00 };
   cblas_zherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zherk(case 1358) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zherk(case 1358) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 113;
   int N = 1;
   int K = 1;
   double alpha = 0;
   double beta = 0;
   double A[] = { -0.824, 0.684 };
   int lda = 1;
   double C[] = { 0.699, -0.615 };
   int ldc = 1;
   double C_expected[] = { 0.000000000000e+00, 0.000000000000e+00 };
   cblas_zherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zherk(case 1359) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zherk(case 1359) imag");
     };
   };
  };


}
