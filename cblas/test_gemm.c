#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_gemm () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int transA = 111;
   int transB = 111;
   int M = 1;
   int N = 1;
   int K = 2;
   float alpha = -0.3;
   float beta = 0.1;
   float A[] = { -0.813, -0.611 };
   int lda = 2;
   float B[] = { 0.736, 0.293 };
   int ldb = 1;
   float C[] = { 0.572 };
   int ldc = 1;
   float C_expected[] = { 0.2904173 };
   cblas_sgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "sgemm(case 1212)");
     }
   };
  };


  {
   int order = 102;
   int transA = 111;
   int transB = 111;
   int M = 1;
   int N = 1;
   int K = 2;
   float alpha = -0.3;
   float beta = 0.1;
   float A[] = { -0.717, 0.259 };
   int lda = 1;
   float B[] = { -0.399, 0.398 };
   int ldb = 2;
   float C[] = { -0.096 };
   int ldc = 1;
   float C_expected[] = { -0.1263495 };
   cblas_sgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "sgemm(case 1213)");
     }
   };
  };


  {
   int order = 101;
   int transA = 111;
   int transB = 112;
   int M = 1;
   int N = 1;
   int K = 2;
   float alpha = 1;
   float beta = -1;
   float A[] = { -0.359, 0.613 };
   int lda = 2;
   float B[] = { -0.534, -0.804 };
   int ldb = 2;
   float C[] = { 0.358 };
   int ldc = 1;
   float C_expected[] = { -0.659146 };
   cblas_sgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "sgemm(case 1214)");
     }
   };
  };


  {
   int order = 102;
   int transA = 111;
   int transB = 112;
   int M = 1;
   int N = 1;
   int K = 2;
   float alpha = 1;
   float beta = -1;
   float A[] = { -0.078, -0.345 };
   int lda = 1;
   float B[] = { 0.369, -0.169 };
   int ldb = 1;
   float C[] = { -0.448 };
   int ldc = 1;
   float C_expected[] = { 0.477523 };
   cblas_sgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "sgemm(case 1215)");
     }
   };
  };


  {
   int order = 101;
   int transA = 112;
   int transB = 111;
   int M = 1;
   int N = 1;
   int K = 2;
   float alpha = 0.1;
   float beta = -1;
   float A[] = { 0.377, 0.498 };
   int lda = 1;
   float B[] = { -0.725, -0.655 };
   int ldb = 1;
   float C[] = { -0.698 };
   int ldc = 1;
   float C_expected[] = { 0.6380485 };
   cblas_sgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "sgemm(case 1216)");
     }
   };
  };


  {
   int order = 102;
   int transA = 112;
   int transB = 111;
   int M = 1;
   int N = 1;
   int K = 2;
   float alpha = 0.1;
   float beta = -1;
   float A[] = { -0.553, 0.237 };
   int lda = 2;
   float B[] = { 0.232, -0.558 };
   int ldb = 2;
   float C[] = { 0.389 };
   int ldc = 1;
   float C_expected[] = { -0.4150542 };
   cblas_sgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "sgemm(case 1217)");
     }
   };
  };


  {
   int order = 101;
   int transA = 112;
   int transB = 112;
   int M = 1;
   int N = 1;
   int K = 2;
   float alpha = 1;
   float beta = 0.1;
   float A[] = { -0.391, 0.467 };
   int lda = 1;
   float B[] = { 0.889, 0.734 };
   int ldb = 2;
   float C[] = { -0.662 };
   int ldc = 1;
   float C_expected[] = { -0.071021 };
   cblas_sgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "sgemm(case 1218)");
     }
   };
  };


  {
   int order = 102;
   int transA = 112;
   int transB = 112;
   int M = 1;
   int N = 1;
   int K = 2;
   float alpha = 1;
   float beta = 0.1;
   float A[] = { 0.151, -0.728 };
   int lda = 2;
   float B[] = { -0.581, 0.949 };
   int ldb = 1;
   float C[] = { -0.955 };
   int ldc = 1;
   float C_expected[] = { -0.874103 };
   cblas_sgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "sgemm(case 1219)");
     }
   };
  };


  {
   int order = 101;
   int transA = 111;
   int transB = 111;
   int M = 1;
   int N = 1;
   int K = 2;
   double alpha = 0.1;
   double beta = -0.3;
   double A[] = { -0.708, 0.575 };
   int lda = 2;
   double B[] = { 0.845, -0.87 };
   int ldb = 1;
   double C[] = { -0.683 };
   int ldc = 1;
   double C_expected[] = { 0.095049 };
   cblas_dgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dgemm(case 1220)");
     }
   };
  };


  {
   int order = 102;
   int transA = 111;
   int transB = 111;
   int M = 1;
   int N = 1;
   int K = 2;
   double alpha = 0.1;
   double beta = -0.3;
   double A[] = { -0.575, 0.032 };
   int lda = 1;
   double B[] = { -0.262, -0.357 };
   int ldb = 2;
   double C[] = { 0.071 };
   int ldc = 1;
   double C_expected[] = { -0.0073774 };
   cblas_dgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dgemm(case 1221)");
     }
   };
  };


  {
   int order = 101;
   int transA = 111;
   int transB = 112;
   int M = 1;
   int N = 1;
   int K = 2;
   double alpha = 0.1;
   double beta = 0;
   double A[] = { -0.375, 0.394 };
   int lda = 2;
   double B[] = { -0.507, -0.919 };
   int ldb = 2;
   double C[] = { -0.994 };
   int ldc = 1;
   double C_expected[] = { -0.0171961 };
   cblas_dgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dgemm(case 1222)");
     }
   };
  };


  {
   int order = 102;
   int transA = 111;
   int transB = 112;
   int M = 1;
   int N = 1;
   int K = 2;
   double alpha = 0.1;
   double beta = 0;
   double A[] = { -0.358, -0.157 };
   int lda = 1;
   double B[] = { 0.334, 0.366 };
   int ldb = 1;
   double C[] = { 0.105 };
   int ldc = 1;
   double C_expected[] = { -0.0177034 };
   cblas_dgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dgemm(case 1223)");
     }
   };
  };


  {
   int order = 101;
   int transA = 112;
   int transB = 111;
   int M = 1;
   int N = 1;
   int K = 2;
   double alpha = 1;
   double beta = 0.1;
   double A[] = { -0.586, 0.507 };
   int lda = 1;
   double B[] = { 0.855, 0.836 };
   int ldb = 1;
   double C[] = { -0.322 };
   int ldc = 1;
   double C_expected[] = { -0.109378 };
   cblas_dgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dgemm(case 1224)");
     }
   };
  };


  {
   int order = 102;
   int transA = 112;
   int transB = 111;
   int M = 1;
   int N = 1;
   int K = 2;
   double alpha = 1;
   double beta = 0.1;
   double A[] = { -0.695, -0.234 };
   int lda = 2;
   double B[] = { -0.002, -0.338 };
   int ldb = 2;
   double C[] = { -0.209 };
   int ldc = 1;
   double C_expected[] = { 0.059582 };
   cblas_dgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dgemm(case 1225)");
     }
   };
  };


  {
   int order = 101;
   int transA = 112;
   int transB = 112;
   int M = 1;
   int N = 1;
   int K = 2;
   double alpha = 0;
   double beta = -1;
   double A[] = { 0.029, -0.662 };
   int lda = 1;
   double B[] = { 0.061, -0.865 };
   int ldb = 2;
   double C[] = { 0.121 };
   int ldc = 1;
   double C_expected[] = { -0.121 };
   cblas_dgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dgemm(case 1226)");
     }
   };
  };


  {
   int order = 102;
   int transA = 112;
   int transB = 112;
   int M = 1;
   int N = 1;
   int K = 2;
   double alpha = 0;
   double beta = -1;
   double A[] = { 0.581, -0.09 };
   int lda = 2;
   double B[] = { -0.788, -0.694 };
   int ldb = 1;
   double C[] = { 0.556 };
   int ldc = 1;
   double C_expected[] = { -0.556 };
   cblas_dgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dgemm(case 1227)");
     }
   };
  };


  {
   int order = 101;
   int transA = 111;
   int transB = 111;
   int M = 1;
   int N = 1;
   int K = 2;
   float alpha[2] = {0, 1};
   float beta[2] = {0, 1};
   float A[] = { -0.849, -0.76, -0.382, 0.437 };
   int lda = 2;
   float B[] = { 0.397, 0.288, -0.826, 0.889 };
   int ldb = 1;
   float C[] = { 0.356, -0.96 };
   int ldc = 1;
   float C_expected[] = { 2.206792, 0.164866 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1228) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1228) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 111;
   int transB = 111;
   int M = 1;
   int N = 1;
   int K = 2;
   float alpha[2] = {0, 1};
   float beta[2] = {0, 1};
   float A[] = { 0.576, -0.392, 0.4, -0.041 };
   int lda = 1;
   float B[] = { 0.511, -0.524, -0.906, 0.154 };
   int ldb = 2;
   float C[] = { -0.809, -0.832 };
   int ldc = 1;
   float C_expected[] = { 1.23539, -1.076158 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1229) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1229) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 111;
   int transB = 112;
   int M = 1;
   int N = 1;
   int K = 2;
   float alpha[2] = {0, 0.1};
   float beta[2] = {0, 1};
   float A[] = { -0.888, -0.749, 0.656, 0.571 };
   int lda = 2;
   float B[] = { 0.465, -0.515, 0.653, 0.35 };
   int ldb = 2;
   float C[] = { -0.297, 0.511 };
   int ldc = 1;
   float C_expected[] = { -0.5821498, -0.3540137 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1230) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1230) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 111;
   int transB = 112;
   int M = 1;
   int N = 1;
   int K = 2;
   float alpha[2] = {0, 0.1};
   float beta[2] = {0, 1};
   float A[] = { -0.933, -0.851, -0.529, 0.532 };
   int lda = 1;
   float B[] = { 0.477, -0.789, 0.469, 0.288 };
   int ldb = 1;
   float C[] = { -0.707, 0.345 };
   int ldc = 1;
   float C_expected[] = { -0.3877366, -0.8587797 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1231) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1231) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 111;
   int transB = 113;
   int M = 1;
   int N = 1;
   int K = 2;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {-1, 0};
   float A[] = { 0.297, 0.451, -0.707, -0.113 };
   int lda = 2;
   float B[] = { 0.983, 0.351, 0.142, -0.62 };
   int ldb = 2;
   float C[] = { 0.647, -0.994 };
   int ldc = 1;
   float C_expected[] = { -0.7614454, 1.0705818 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1232) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1232) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 111;
   int transB = 113;
   int M = 1;
   int N = 1;
   int K = 2;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {-1, 0};
   float A[] = { -0.991, -0.517, -0.909, 0.658 };
   int lda = 1;
   float B[] = { -0.99, -0.138, -0.071, -0.104 };
   int ldb = 1;
   float C[] = { 0.511, -0.186 };
   int ldc = 1;
   float C_expected[] = { -0.8489447, 0.2207089 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1233) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1233) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 112;
   int transB = 111;
   int M = 1;
   int N = 1;
   int K = 2;
   float alpha[2] = {-1, 0};
   float beta[2] = {0, 0.1};
   float A[] = { 0.155, 0.03, 0.122, -0.768 };
   int lda = 1;
   float B[] = { -0.741, -0.515, 0.626, 0.196 };
   int ldb = 1;
   float C[] = { -0.794, 0.414 };
   int ldc = 1;
   float C_expected[] = { -0.168895, 0.479511 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1234) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1234) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 112;
   int transB = 111;
   int M = 1;
   int N = 1;
   int K = 2;
   float alpha[2] = {-1, 0};
   float beta[2] = {0, 0.1};
   float A[] = { -0.849, 0.67, 0.146, 0.143 };
   int lda = 2;
   float B[] = { 0.073, -0.368, -0.9, -0.583 };
   int ldb = 2;
   float C[] = { -0.982, 0.445 };
   int ldc = 1;
   float C_expected[] = { -0.181052, -0.245724 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1235) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1235) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 112;
   int transB = 112;
   int M = 1;
   int N = 1;
   int K = 2;
   float alpha[2] = {-1, 0};
   float beta[2] = {0, 0};
   float A[] = { 0.355, 0.985, -0.989, 0.63 };
   int lda = 1;
   float B[] = { -0.528, 0.106, -0.319, 0.962 };
   int ldb = 2;
   float C[] = { -0.397, 0.004 };
   int ldc = 1;
   float C_expected[] = { 0.582419, 1.634838 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1236) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1236) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 112;
   int transB = 112;
   int M = 1;
   int N = 1;
   int K = 2;
   float alpha[2] = {-1, 0};
   float beta[2] = {0, 0};
   float A[] = { 0.874, 0.415, 0.636, 0.315 };
   int lda = 2;
   float B[] = { -0.199, -0.867, 0.386, -0.259 };
   int ldb = 1;
   float C[] = { -0.759, -0.053 };
   int ldc = 1;
   float C_expected[] = { -0.51296, 0.883477 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1237) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1237) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 112;
   int transB = 113;
   int M = 1;
   int N = 1;
   int K = 2;
   float alpha[2] = {1, 0};
   float beta[2] = {-0.3, 0.1};
   float A[] = { 0.987, -0.329, -0.138, 0.463 };
   int lda = 1;
   float B[] = { -0.429, 0.842, -0.899, 0.192 };
   int ldb = 2;
   float C[] = { -0.475, 0.916 };
   int ldc = 1;
   float C_expected[] = { -0.436583, -1.401954 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1238) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1238) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 112;
   int transB = 113;
   int M = 1;
   int N = 1;
   int K = 2;
   float alpha[2] = {1, 0};
   float beta[2] = {-0.3, 0.1};
   float A[] = { 0.542, 0.547, 0.461, 0.583 };
   int lda = 2;
   float B[] = { -0.302, 0.216, 0.434, 0.192 };
   int ldb = 1;
   float C[] = { -0.162, -0.108 };
   int ldc = 1;
   float C_expected[] = { 0.325878, -0.101556 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1239) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1239) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 113;
   int transB = 111;
   int M = 1;
   int N = 1;
   int K = 2;
   float alpha[2] = {0, 0.1};
   float beta[2] = {1, 0};
   float A[] = { -0.101, -0.216, 0.889, -0.623 };
   int lda = 1;
   float B[] = { 0.048, -0.345, 0.293, -0.851 };
   int ldb = 1;
   float C[] = { 0.693, -0.22 };
   int ldc = 1;
   float C_expected[] = { 0.7458787, -0.1339678 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1240) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1240) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 113;
   int transB = 111;
   int M = 1;
   int N = 1;
   int K = 2;
   float alpha[2] = {0, 0.1};
   float beta[2] = {1, 0};
   float A[] = { 0.439, -0.952, -0.484, -0.508 };
   int lda = 2;
   float B[] = { 0.381, -0.192, -0.889, -0.279 };
   int ldb = 2;
   float C[] = { -0.708, 0.298 };
   int ldc = 1;
   float C_expected[] = { -0.7041848, 0.3902051 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1241) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1241) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 113;
   int transB = 112;
   int M = 1;
   int N = 1;
   int K = 2;
   float alpha[2] = {1, 0};
   float beta[2] = {1, 0};
   float A[] = { 0.94, 0.984, -0.091, -0.276 };
   int lda = 1;
   float B[] = { -0.342, -0.665, -0.484, -0.2 };
   int ldb = 2;
   float C[] = { -0.688, -0.203 };
   int ldc = 1;
   float C_expected[] = { -1.564596, -0.606956 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1242) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1242) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 113;
   int transB = 112;
   int M = 1;
   int N = 1;
   int K = 2;
   float alpha[2] = {1, 0};
   float beta[2] = {1, 0};
   float A[] = { -0.501, 0.566, -0.562, 0.955 };
   int lda = 2;
   float B[] = { -0.086, 0.177, -0.856, 0.974 };
   int ldb = 1;
   float C[] = { 0.349, 0.023 };
   int ldc = 1;
   float C_expected[] = { 1.90351, 0.253091 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1243) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1243) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 113;
   int transB = 113;
   int M = 1;
   int N = 1;
   int K = 2;
   float alpha[2] = {0, 0.1};
   float beta[2] = {0, 0};
   float A[] = { -0.818, 0.188, 0.337, 0.6 };
   int lda = 1;
   float B[] = { -0.743, 0.449, -0.803, -0.681 };
   int ldb = 2;
   float C[] = { 0.77, -0.77 };
   int ldc = 1;
   float C_expected[] = { -0.1218263, 0.0661351 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1244) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1244) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 113;
   int transB = 113;
   int M = 1;
   int N = 1;
   int K = 2;
   float alpha[2] = {0, 0.1};
   float beta[2] = {0, 0};
   float A[] = { 0.553, 0.481, -0.851, 0.155 };
   int lda = 2;
   float B[] = { 0.824, 0.42, -0.987, 0.147 };
   int ldb = 1;
   float C[] = { -0.464, -0.029 };
   int ldc = 1;
   float C_expected[] = { 0.0350522, 0.1070804 };
   cblas_cgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cgemm(case 1245) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cgemm(case 1245) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 111;
   int transB = 111;
   int M = 1;
   int N = 1;
   int K = 2;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {0, 1};
   double A[] = { 0.456, 0.758, -0.45, -0.357 };
   int lda = 2;
   double B[] = { -0.191, -0.573, 0.63, 0.146 };
   int ldb = 1;
   double C[] = { 0.224, 0.501 };
   int ldc = 1;
   double C_expected[] = { -0.4660904, 0.4445888 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1246) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1246) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 111;
   int transB = 111;
   int M = 1;
   int N = 1;
   int K = 2;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {0, 1};
   double A[] = { -0.46, -0.236, 0.308, 0.069 };
   int lda = 1;
   double B[] = { -0.003, -0.814, 0.605, -0.755 };
   int ldb = 2;
   double C[] = { -0.889, 0.796 };
   int ldc = 1;
   double C_expected[] = { -0.8287486, -0.9395348 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1247) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1247) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 111;
   int transB = 112;
   int M = 1;
   int N = 1;
   int K = 2;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {-1, 0};
   double A[] = { -0.512, -0.225, 0.234, 0.043 };
   int lda = 2;
   double B[] = { -0.094, 0.491, 0.257, -0.787 };
   int ldb = 2;
   double C[] = { -0.27, 0.545 };
   int ldc = 1;
   double C_expected[] = { 0.2345603, -0.3987371 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1248) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1248) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 111;
   int transB = 112;
   int M = 1;
   int N = 1;
   int K = 2;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {-1, 0};
   double A[] = { 0.001, 0.411, -0.002, -0.912 };
   int lda = 1;
   double B[] = { -0.894, 0.811, -0.308, -0.116 };
   int ldb = 1;
   double C[] = { 0.732, -0.085 };
   int ldc = 1;
   double C_expected[] = { -0.5916332, 0.0667094 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1249) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1249) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 111;
   int transB = 113;
   int M = 1;
   int N = 1;
   int K = 2;
   double alpha[2] = {0, 0};
   double beta[2] = {-0.3, 0.1};
   double A[] = { -0.093, 0.465, 0.842, -0.058 };
   int lda = 2;
   double B[] = { 0.093, 0.724, 0.047, 0.991 };
   int ldb = 2;
   double C[] = { -0.438, -0.188 };
   int ldc = 1;
   double C_expected[] = { 0.1502, 0.0126 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1250) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1250) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 111;
   int transB = 113;
   int M = 1;
   int N = 1;
   int K = 2;
   double alpha[2] = {0, 0};
   double beta[2] = {-0.3, 0.1};
   double A[] = { -0.753, 0.988, -0.972, 0.164 };
   int lda = 1;
   double B[] = { 0.065, 0.694, -0.746, -0.082 };
   int ldb = 1;
   double C[] = { -0.06, 0.302 };
   int ldc = 1;
   double C_expected[] = { -0.0122, -0.0966 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1251) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1251) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 112;
   int transB = 111;
   int M = 1;
   int N = 1;
   int K = 2;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {0, 0};
   double A[] = { -0.124, 0.331, 0.753, 0.504 };
   int lda = 1;
   double B[] = { -0.888, 0.982, -0.361, 0.612 };
   int ldb = 1;
   double C[] = { -0.212, 0.189 };
   int ldc = 1;
   double C_expected[] = { 0.2522437, -0.0384799 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1252) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1252) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 112;
   int transB = 111;
   int M = 1;
   int N = 1;
   int K = 2;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {0, 0};
   double A[] = { 0.727, -0.673, 0.116, -0.277 };
   int lda = 2;
   double B[] = { -0.909, 0.064, -0.124, -0.81 };
   int ldb = 2;
   double C[] = { -0.046, 0.82 };
   int ldc = 1;
   double C_expected[] = { 0.1970902, -0.2652544 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1253) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1253) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 112;
   int transB = 112;
   int M = 1;
   int N = 1;
   int K = 2;
   double alpha[2] = {0, 0.1};
   double beta[2] = {0, 1};
   double A[] = { -0.374, -0.794, 0.792, -0.281 };
   int lda = 1;
   double B[] = { 0.246, -0.231, 0.853, 0.373 };
   int ldb = 2;
   double C[] = { -0.779, -0.103 };
   int ldc = 1;
   double C_expected[] = { 0.1083207, -0.7285029 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1254) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1254) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 112;
   int transB = 112;
   int M = 1;
   int N = 1;
   int K = 2;
   double alpha[2] = {0, 0.1};
   double beta[2] = {0, 1};
   double A[] = { 0.884, 0.527, -0.055, 0.762 };
   int lda = 2;
   double B[] = { 0.891, -0.717, 0.384, -0.338 };
   int ldb = 1;
   double C[] = { -0.957, -0.308 };
   int ldc = 1;
   double C_expected[] = { 0.2933073, -0.8168061 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1255) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1255) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 112;
   int transB = 113;
   int M = 1;
   int N = 1;
   int K = 2;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {-0.3, 0.1};
   double A[] = { 0.282, 0.002, 0.188, -0.363 };
   int lda = 1;
   double B[] = { -0.646, -0.02, -0.407, -0.874 };
   int ldb = 2;
   double C[] = { -0.649, 0.636 };
   int ldc = 1;
   double C_expected[] = { 0.0818997, -0.3447669 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1256) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1256) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 112;
   int transB = 113;
   int M = 1;
   int N = 1;
   int K = 2;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {-0.3, 0.1};
   double A[] = { -0.795, 0.084, 0.95, 0.816 };
   int lda = 2;
   double B[] = { 0.099, 0.191, -0.906, -0.968 };
   int ldb = 1;
   double C[] = { -0.328, 0.79 };
   int ldc = 1;
   double C_expected[] = { 0.4993282, -0.5432644 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1257) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1257) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 113;
   int transB = 111;
   int M = 1;
   int N = 1;
   int K = 2;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {0, 1};
   double A[] = { 0.326, 0.945, -0.084, 0.428 };
   int lda = 1;
   double B[] = { -0.587, 0.841, 0.14, -0.213 };
   int ldb = 1;
   double C[] = { -0.591, 0.707 };
   int ldc = 1;
   double C_expected[] = { -0.935823, -0.77701 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1258) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1258) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 113;
   int transB = 111;
   int M = 1;
   int N = 1;
   int K = 2;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {0, 1};
   double A[] = { 0.968, 0.243, -0.956, -0.274 };
   int lda = 2;
   double B[] = { 0.393, 0.356, -0.221, -0.062 };
   int ldb = 2;
   double C[] = { 0.885, -0.765 };
   int ldc = 1;
   double C_expected[] = { 0.5316585, 0.8801715 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1259) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1259) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 113;
   int transB = 112;
   int M = 1;
   int N = 1;
   int K = 2;
   double alpha[2] = {1, 0};
   double beta[2] = {-0.3, 0.1};
   double A[] = { -0.928, 0.362, 0.365, -0.612 };
   int lda = 1;
   double B[] = { 0.503, -0.718, -0.849, -0.923 };
   int ldb = 2;
   double C[] = { 0.501, 0.019 };
   int ldc = 1;
   double C_expected[] = { -0.623909, -0.327865 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1260) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1260) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 113;
   int transB = 112;
   int M = 1;
   int N = 1;
   int K = 2;
   double alpha[2] = {1, 0};
   double beta[2] = {-0.3, 0.1};
   double A[] = { 0.461, -0.876, 0.153, -0.519 };
   int lda = 2;
   double B[] = { 0.526, -0.223, -0.697, -0.077 };
   int ldb = 1;
   double C[] = { 0.904, -0.367 };
   int ldc = 1;
   double C_expected[] = { 0.136656, 0.184949 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1261) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1261) imag");
     };
   };
  };


  {
   int order = 101;
   int transA = 113;
   int transB = 113;
   int M = 1;
   int N = 1;
   int K = 2;
   double alpha[2] = {-1, 0};
   double beta[2] = {-1, 0};
   double A[] = { 0.869, -0.921, -0.859, -0.091 };
   int lda = 1;
   double B[] = { 0.008, -0.122, -0.321, 0.602 };
   int ldb = 2;
   double C[] = { -0.807, -0.748 };
   int ldc = 1;
   double C_expected[] = { 0.581889, 0.146707 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1262) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1262) imag");
     };
   };
  };


  {
   int order = 102;
   int transA = 113;
   int transB = 113;
   int M = 1;
   int N = 1;
   int K = 2;
   double alpha[2] = {-1, 0};
   double beta[2] = {-1, 0};
   double A[] = { 0.254, -0.182, -0.271, 0.263 };
   int lda = 2;
   double B[] = { -0.595, -0.042, -0.672, -0.255 };
   int ldb = 1;
   double C[] = { -0.364, -0.96 };
   int ldc = 1;
   double C_expected[] = { 0.273597, 0.949991 };
   cblas_zgemm(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zgemm(case 1263) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zgemm(case 1263) imag");
     };
   };
  };


}
