#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_symm () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int M = 1;
   int N = 1;
   float alpha = 0;
   float beta = 0;
   float A[] = { -0.595 };
   int lda = 1;
   float B[] = { 0.884 };
   int ldb = 1;
   float C[] = { 0.421 };
   int ldc = 1;
   float C_expected[] = { 0.000000000000e+00 };
   cblas_ssymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssymm(case 1264)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int M = 1;
   int N = 1;
   float alpha = 0;
   float beta = 0;
   float A[] = { -0.387 };
   int lda = 1;
   float B[] = { -0.498 };
   int ldb = 1;
   float C[] = { -0.407 };
   int ldc = 1;
   float C_expected[] = { 0.000000000000e+00 };
   cblas_ssymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssymm(case 1265)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int M = 1;
   int N = 1;
   float alpha = 1;
   float beta = 0;
   float A[] = { -0.679 };
   int lda = 1;
   float B[] = { -0.743 };
   int ldb = 1;
   float C[] = { 0.904 };
   int ldc = 1;
   float C_expected[] = { 0.504497 };
   cblas_ssymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssymm(case 1266)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int M = 1;
   int N = 1;
   float alpha = 1;
   float beta = 0;
   float A[] = { 0.157 };
   int lda = 1;
   float B[] = { -0.078 };
   int ldb = 1;
   float C[] = { 0.77 };
   int ldc = 1;
   float C_expected[] = { -0.012246 };
   cblas_ssymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssymm(case 1267)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int M = 1;
   int N = 1;
   float alpha = -1;
   float beta = 0;
   float A[] = { -0.406 };
   int lda = 1;
   float B[] = { 0.565 };
   int ldb = 1;
   float C[] = { -0.641 };
   int ldc = 1;
   float C_expected[] = { 0.22939 };
   cblas_ssymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssymm(case 1268)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int M = 1;
   int N = 1;
   float alpha = -1;
   float beta = 0;
   float A[] = { -0.428 };
   int lda = 1;
   float B[] = { -0.34 };
   int ldb = 1;
   float C[] = { -0.948 };
   int ldc = 1;
   float C_expected[] = { -0.14552 };
   cblas_ssymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssymm(case 1269)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int M = 1;
   int N = 1;
   float alpha = 0;
   float beta = 1;
   float A[] = { 0.544 };
   int lda = 1;
   float B[] = { -0.601 };
   int ldb = 1;
   float C[] = { 0.3 };
   int ldc = 1;
   float C_expected[] = { 0.3 };
   cblas_ssymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssymm(case 1270)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int M = 1;
   int N = 1;
   float alpha = 0;
   float beta = 1;
   float A[] = { -0.76 };
   int lda = 1;
   float B[] = { 0.728 };
   int ldb = 1;
   float C[] = { -0.852 };
   int ldc = 1;
   float C_expected[] = { -0.852 };
   cblas_ssymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssymm(case 1271)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int M = 1;
   int N = 1;
   double alpha = 0.1;
   double beta = -0.3;
   double A[] = { -0.935 };
   int lda = 1;
   double B[] = { 0.346 };
   int ldb = 1;
   double C[] = { 0.209 };
   int ldc = 1;
   double C_expected[] = { -0.095051 };
   cblas_dsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsymm(case 1272)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int M = 1;
   int N = 1;
   double alpha = 0.1;
   double beta = -0.3;
   double A[] = { -0.066 };
   int lda = 1;
   double B[] = { 0.975 };
   int ldb = 1;
   double C[] = { -0.412 };
   int ldc = 1;
   double C_expected[] = { 0.117165 };
   cblas_dsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsymm(case 1273)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int M = 1;
   int N = 1;
   double alpha = 1;
   double beta = 0;
   double A[] = { 0.886 };
   int lda = 1;
   double B[] = { 0.486 };
   int ldb = 1;
   double C[] = { -0.695 };
   int ldc = 1;
   double C_expected[] = { 0.430596 };
   cblas_dsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsymm(case 1274)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int M = 1;
   int N = 1;
   double alpha = 1;
   double beta = 0;
   double A[] = { -0.894 };
   int lda = 1;
   double B[] = { -0.288 };
   int ldb = 1;
   double C[] = { 0.629 };
   int ldc = 1;
   double C_expected[] = { 0.257472 };
   cblas_dsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsymm(case 1275)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int M = 1;
   int N = 1;
   double alpha = -0.3;
   double beta = -1;
   double A[] = { -0.785 };
   int lda = 1;
   double B[] = { -0.731 };
   int ldb = 1;
   double C[] = { 0.022 };
   int ldc = 1;
   double C_expected[] = { -0.1941505 };
   cblas_dsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsymm(case 1276)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int M = 1;
   int N = 1;
   double alpha = -0.3;
   double beta = -1;
   double A[] = { 0.173 };
   int lda = 1;
   double B[] = { 0.079 };
   int ldb = 1;
   double C[] = { -0.243 };
   int ldc = 1;
   double C_expected[] = { 0.2388999 };
   cblas_dsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsymm(case 1277)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int M = 1;
   int N = 1;
   double alpha = 1;
   double beta = 1;
   double A[] = { 0.153 };
   int lda = 1;
   double B[] = { 0.193 };
   int ldb = 1;
   double C[] = { 0.847 };
   int ldc = 1;
   double C_expected[] = { 0.876529 };
   cblas_dsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsymm(case 1278)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int M = 1;
   int N = 1;
   double alpha = 1;
   double beta = 1;
   double A[] = { -0.319 };
   int lda = 1;
   double B[] = { -0.24 };
   int ldb = 1;
   double C[] = { -0.166 };
   int ldc = 1;
   double C_expected[] = { -0.08944 };
   cblas_dsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsymm(case 1279)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int M = 1;
   int N = 1;
   float alpha[2] = {1, 0};
   float beta[2] = {-1, 0};
   float A[] = { 0.094, -0.103 };
   int lda = 1;
   float B[] = { -0.151, 0.954 };
   int ldb = 1;
   float C[] = { -0.093, -0.537 };
   int ldc = 1;
   float C_expected[] = { 0.177068, 0.642229 };
   cblas_csymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csymm(case 1280) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csymm(case 1280) imag");
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
   float beta[2] = {-1, 0};
   float A[] = { -0.478, -0.391 };
   int lda = 1;
   float B[] = { 0.031, 0.098 };
   int ldb = 1;
   float C[] = { 0.029, -0.256 };
   int ldc = 1;
   float C_expected[] = { -0.0055, 0.197035 };
   cblas_csymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csymm(case 1281) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csymm(case 1281) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 1};
   float beta[2] = {-0.3, 0.1};
   float A[] = { -0.782, -0.842 };
   int lda = 1;
   float B[] = { -0.168, -0.692 };
   int ldb = 1;
   float C[] = { 0.662, -0.311 };
   int ldc = 1;
   float C_expected[] = { -0.8501, -0.291788 };
   cblas_csymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csymm(case 1282) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csymm(case 1282) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 1};
   float beta[2] = {-0.3, 0.1};
   float A[] = { -0.687, 0.24 };
   int lda = 1;
   float B[] = { -0.345, 0.904 };
   int ldb = 1;
   float C[] = { 0.009, -0.981 };
   int ldc = 1;
   float C_expected[] = { 0.799248, 0.315255 };
   cblas_csymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csymm(case 1283) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csymm(case 1283) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 1};
   float beta[2] = {0, 0};
   float A[] = { 0.019, 0.861 };
   int lda = 1;
   float B[] = { 0.873, 0.398 };
   int ldb = 1;
   float C[] = { 0.113, -0.707 };
   int ldc = 1;
   float C_expected[] = { -0.759215, -0.326091 };
   cblas_csymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csymm(case 1284) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csymm(case 1284) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int M = 1;
   int N = 1;
   float alpha[2] = {0, 1};
   float beta[2] = {0, 0};
   float A[] = { 0.471, 0.878 };
   int lda = 1;
   float B[] = { -0.441, 0.979 };
   int ldb = 1;
   float C[] = { 0.214, -0.509 };
   int ldc = 1;
   float C_expected[] = { -0.073911, -1.067273 };
   cblas_csymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csymm(case 1285) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csymm(case 1285) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int M = 1;
   int N = 1;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {0, 0};
   float A[] = { -0.023, -0.075 };
   int lda = 1;
   float B[] = { 0.627, -0.312 };
   int ldb = 1;
   float C[] = { -0.911, 0.807 };
   int ldc = 1;
   float C_expected[] = { 0.0153312, 0.0081726 };
   cblas_csymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csymm(case 1286) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csymm(case 1286) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int M = 1;
   int N = 1;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {0, 0};
   float A[] = { 0.323, -0.578 };
   int lda = 1;
   float B[] = { 0.797, 0.545 };
   int ldb = 1;
   float C[] = { 0.031, 0.308 };
   int ldc = 1;
   float C_expected[] = { -0.1432692, 0.1426334 };
   cblas_csymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csymm(case 1287) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csymm(case 1287) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int M = 1;
   int N = 1;
   double alpha[2] = {0, 0.1};
   double beta[2] = {0, 0.1};
   double A[] = { 0.755, 0.268 };
   int lda = 1;
   double B[] = { -0.99, 0.02 };
   int ldb = 1;
   double C[] = { 0.606, 0.727 };
   int ldc = 1;
   double C_expected[] = { -0.047678, -0.014681 };
   cblas_zsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsymm(case 1288) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsymm(case 1288) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int M = 1;
   int N = 1;
   double alpha[2] = {0, 0.1};
   double beta[2] = {0, 0.1};
   double A[] = { -0.568, -0.888 };
   int lda = 1;
   double B[] = { 0.281, -0.779 };
   int ldb = 1;
   double C[] = { -0.801, 0.83 };
   int ldc = 1;
   double C_expected[] = { -0.1022944, -0.165236 };
   cblas_zsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsymm(case 1289) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsymm(case 1289) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int M = 1;
   int N = 1;
   double alpha[2] = {1, 0};
   double beta[2] = {0, 1};
   double A[] = { 0.166, 0.808 };
   int lda = 1;
   double B[] = { 0.723, 0.9 };
   int ldb = 1;
   double C[] = { 0.053, -0.757 };
   int ldc = 1;
   double C_expected[] = { 0.149818, 0.786584 };
   cblas_zsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsymm(case 1290) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsymm(case 1290) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int M = 1;
   int N = 1;
   double alpha[2] = {1, 0};
   double beta[2] = {0, 1};
   double A[] = { -0.909, -0.111 };
   int lda = 1;
   double B[] = { 0.339, 0.908 };
   int ldb = 1;
   double C[] = { 0.99, -0.578 };
   int ldc = 1;
   double C_expected[] = { 0.370637, 0.126999 };
   cblas_zsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsymm(case 1291) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsymm(case 1291) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int M = 1;
   int N = 1;
   double alpha[2] = {0, 0};
   double beta[2] = {1, 0};
   double A[] = { 0.993, -0.653 };
   int lda = 1;
   double B[] = { -0.502, 0.796 };
   int ldb = 1;
   double C[] = { 0.097, 0.236 };
   int ldc = 1;
   double C_expected[] = { 0.097, 0.236 };
   cblas_zsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsymm(case 1292) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsymm(case 1292) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int M = 1;
   int N = 1;
   double alpha[2] = {0, 0};
   double beta[2] = {1, 0};
   double A[] = { 0.339, -0.122 };
   int lda = 1;
   double B[] = { 0.735, 0.201 };
   int ldb = 1;
   double C[] = { -0.35, -0.269 };
   int ldc = 1;
   double C_expected[] = { -0.35, -0.269 };
   cblas_zsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsymm(case 1293) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsymm(case 1293) imag");
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
   double beta[2] = {0, 1};
   double A[] = { -0.053, -0.826 };
   int lda = 1;
   double B[] = { 0.67, -0.613 };
   int ldb = 1;
   double C[] = { -0.098, -0.737 };
   int ldc = 1;
   double C_expected[] = { 1.278848, 0.422931 };
   cblas_zsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsymm(case 1294) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsymm(case 1294) imag");
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
   double beta[2] = {0, 1};
   double A[] = { -0.398, -0.204 };
   int lda = 1;
   double B[] = { -0.934, 0.183 };
   int ldb = 1;
   double C[] = { -0.857, -0.927 };
   int ldc = 1;
   double C_expected[] = { 0.517936, -0.974702 };
   cblas_zsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsymm(case 1295) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsymm(case 1295) imag");
     };
   };
  };


}
