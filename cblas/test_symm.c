#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

#include "tests.h"

void
test_symm (void) {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int M = 1;
   int N = 2;
   float alpha = -0.3;
   float beta = -1;
   float A[] = { -0.581 };
   int lda = 1;
   float B[] = { 0.157, 0.451 };
   int ldb = 2;
   float C[] = { -0.869, -0.871 };
   int ldc = 2;
   float C_expected[] = { 0.896365, 0.949609 };
   cblas_ssymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssymm(case 1518)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int M = 1;
   int N = 2;
   float alpha = -0.3;
   float beta = -1;
   float A[] = { 0.874 };
   int lda = 1;
   float B[] = { 0.085, 0.069 };
   int ldb = 1;
   float C[] = { -0.495, -0.828 };
   int ldc = 1;
   float C_expected[] = { 0.472713, 0.809908 };
   cblas_ssymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssymm(case 1519)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int M = 1;
   int N = 2;
   float alpha = -1;
   float beta = 0;
   float A[] = { -0.671, -0.343, 0.6, 0.177 };
   int lda = 2;
   float B[] = { 0.043, 0.01 };
   int ldb = 2;
   float C[] = { 0.988, 0.478 };
   int ldc = 2;
   float C_expected[] = { 0.032283, 0.012979 };
   cblas_ssymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssymm(case 1520)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int M = 1;
   int N = 2;
   float alpha = -1;
   float beta = 0;
   float A[] = { 0.069, 0.096, 0.139, -0.044 };
   int lda = 2;
   float B[] = { -0.448, 0.07 };
   int ldb = 1;
   float C[] = { 0.361, 0.995 };
   int ldc = 1;
   float C_expected[] = { 0.021182, 0.065352 };
   cblas_ssymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssymm(case 1521)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int M = 1;
   int N = 2;
   float alpha = 0;
   float beta = -0.3;
   float A[] = { 0.745 };
   int lda = 1;
   float B[] = { -0.269, 0.448 };
   int ldb = 2;
   float C[] = { -0.986, 0.2 };
   int ldc = 2;
   float C_expected[] = { 0.2958, -0.06 };
   cblas_ssymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssymm(case 1522)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int M = 1;
   int N = 2;
   float alpha = 0;
   float beta = -0.3;
   float A[] = { 0.96 };
   int lda = 1;
   float B[] = { 0.392, -0.07 };
   int ldb = 1;
   float C[] = { -0.235, 0.554 };
   int ldc = 1;
   float C_expected[] = { 0.0705, -0.1662 };
   cblas_ssymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssymm(case 1523)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int M = 1;
   int N = 2;
   float alpha = -0.3;
   float beta = 0.1;
   float A[] = { -0.839, 0.498, -0.215, -0.314 };
   int lda = 2;
   float B[] = { -0.66, 0.593 };
   int ldb = 2;
   float C[] = { -0.806, 0.525 };
   int ldc = 2;
   float C_expected[] = { -0.208474, 0.0657906 };
   cblas_ssymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssymm(case 1524)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int M = 1;
   int N = 2;
   float alpha = -0.3;
   float beta = 0.1;
   float A[] = { 0.994, -0.117, -0.639, 0.925 };
   int lda = 2;
   float B[] = { -0.478, 0.147 };
   int ldb = 1;
   float C[] = { -0.814, 0.316 };
   int ldc = 1;
   float C_expected[] = { 0.0662993, -0.0259703 };
   cblas_ssymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssymm(case 1525)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int M = 1;
   int N = 2;
   double alpha = -0.3;
   double beta = 1;
   double A[] = { -0.981 };
   int lda = 1;
   double B[] = { -0.823, 0.83 };
   int ldb = 2;
   double C[] = { 0.991, 0.382 };
   int ldc = 2;
   double C_expected[] = { 0.7487911, 0.626269 };
   cblas_dsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsymm(case 1526)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int M = 1;
   int N = 2;
   double alpha = -0.3;
   double beta = 1;
   double A[] = { -0.248 };
   int lda = 1;
   double B[] = { 0.74, 0.068 };
   int ldb = 1;
   double C[] = { -0.905, 0.742 };
   int ldc = 1;
   double C_expected[] = { -0.849944, 0.7470592 };
   cblas_dsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsymm(case 1527)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int M = 1;
   int N = 2;
   double alpha = -1;
   double beta = 1;
   double A[] = { 0.591, -0.01, -0.192, -0.376 };
   int lda = 2;
   double B[] = { 0.561, 0.946 };
   int ldb = 2;
   double C[] = { 0.763, 0.189 };
   int ldc = 2;
   double C_expected[] = { 0.440909, 0.550306 };
   cblas_dsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsymm(case 1528)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int M = 1;
   int N = 2;
   double alpha = -1;
   double beta = 1;
   double A[] = { -0.786, 0.87, 0.222, -0.043 };
   int lda = 2;
   double B[] = { -0.503, -0.526 };
   int ldb = 1;
   double C[] = { -0.027, -0.391 };
   int ldc = 1;
   double C_expected[] = { -0.305586, -0.301952 };
   cblas_dsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsymm(case 1529)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int M = 1;
   int N = 2;
   double alpha = 0.1;
   double beta = 0.1;
   double A[] = { -0.468 };
   int lda = 1;
   double B[] = { -0.881, 0.692 };
   int ldb = 2;
   double C[] = { -0.812, -0.395 };
   int ldc = 2;
   double C_expected[] = { -0.0399692, -0.0718856 };
   cblas_dsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsymm(case 1530)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int M = 1;
   int N = 2;
   double alpha = 0.1;
   double beta = 0.1;
   double A[] = { 0.849 };
   int lda = 1;
   double B[] = { -0.887, 0.518 };
   int ldb = 1;
   double C[] = { 0.414, -0.251 };
   int ldc = 1;
   double C_expected[] = { -0.0339063, 0.0188782 };
   cblas_dsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsymm(case 1531)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int M = 1;
   int N = 2;
   double alpha = -1;
   double beta = 1;
   double A[] = { 0.457, 0.624, 0.807, 0.349 };
   int lda = 2;
   double B[] = { -0.609, 0.03 };
   int ldb = 2;
   double C[] = { 0.719, -0.624 };
   int ldc = 2;
   double C_expected[] = { 0.973103, -0.143007 };
   cblas_dsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsymm(case 1532)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int M = 1;
   int N = 2;
   double alpha = -1;
   double beta = 1;
   double A[] = { -0.133, -0.117, -0.163, 0.795 };
   int lda = 2;
   double B[] = { -0.882, 0.549 };
   int ldb = 1;
   double C[] = { 0.715, -0.327 };
   int ldc = 1;
   double C_expected[] = { 0.661927, -0.866649 };
   cblas_dsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsymm(case 1533)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int M = 1;
   int N = 2;
   float alpha[2] = {-1, 0};
   float beta[2] = {1, 0};
   float A[] = { 0.476, 0.816 };
   int lda = 1;
   float B[] = { 0.282, 0.852, -0.891, -0.588 };
   int ldb = 2;
   float C[] = { 0.9, 0.486, -0.78, -0.637 };
   int ldc = 2;
   float C_expected[] = { 1.461, -0.149664, -0.835692, 0.369944 };
   cblas_csymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csymm(case 1534) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csymm(case 1534) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int M = 1;
   int N = 2;
   float alpha[2] = {-1, 0};
   float beta[2] = {1, 0};
   float A[] = { 0.048, 0.172 };
   int lda = 1;
   float B[] = { 0.786, 0.783, 0.809, -0.569 };
   int ldb = 1;
   float C[] = { -0.227, -0.215, 0.881, 0.233 };
   int ldc = 1;
   float C_expected[] = { -0.130052, -0.387776, 0.7443, 0.121164 };
   cblas_csymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csymm(case 1535) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csymm(case 1535) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int M = 1;
   int N = 2;
   float alpha[2] = {0, 1};
   float beta[2] = {0, 1};
   float A[] = { -0.495, -0.012, 0.843, -0.986, -0.243, 0.833, 0.921, 0.004 };
   int lda = 2;
   float B[] = { 0.876, 0.612, 0.805, -0.57 };
   int ldb = 2;
   float C[] = { 0.938, -0.24, -0.874, -0.062 };
   int ldc = 2;
   float C_expected[] = { 1.82769, 0.628319, 0.93157, 1.21158 };
   cblas_csymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csymm(case 1536) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csymm(case 1536) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int M = 1;
   int N = 2;
   float alpha[2] = {0, 1};
   float beta[2] = {0, 1};
   float A[] = { -0.812, 0.83, 0.705, 0.15, -0.463, 0.901, -0.547, -0.483 };
   int lda = 2;
   float B[] = { -0.808, -0.664, 0.352, -0.102 };
   int ldb = 1;
   float C[] = { -0.64, 0.399, 0.896, -0.163 };
   int ldc = 1;
   float C_expected[] = { -0.631906, 0.496142, 0.697798, 1.62656 };
   cblas_csymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csymm(case 1537) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csymm(case 1537) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int M = 1;
   int N = 2;
   float alpha[2] = {-1, 0};
   float beta[2] = {0, 1};
   float A[] = { 0.342, -0.906 };
   int lda = 1;
   float B[] = { 0.676, 0.863, -0.517, -0.138 };
   int ldb = 2;
   float C[] = { 0.274, 0.388, -0.271, 0.205 };
   int ldc = 2;
   float C_expected[] = { -1.40107, 0.59131, 0.096842, -0.692206 };
   cblas_csymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csymm(case 1538) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csymm(case 1538) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int M = 1;
   int N = 2;
   float alpha[2] = {-1, 0};
   float beta[2] = {0, 1};
   float A[] = { 0.418, 0.354 };
   int lda = 1;
   float B[] = { -0.74, 0.018, 0.395, 0.248 };
   int ldb = 1;
   float C[] = { -0.162, 0.175, -0.853, 0.652 };
   int ldc = 1;
   float C_expected[] = { 0.140692, 0.092436, -0.729318, -1.09649 };
   cblas_csymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csymm(case 1539) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csymm(case 1539) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int M = 1;
   int N = 2;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {0, 0.1};
   float A[] = { 0.12, 0.496, 0.313, -0.136, 0.987, 0.532, 0.58, -0.687 };
   int lda = 2;
   float B[] = { -0.587, 0.278, 0.857, 0.136 };
   int ldb = 2;
   float C[] = { 0.162, 0.249, -0.665, 0.456 };
   int ldc = 2;
   float C_expected[] = { -0.22769, -0.0269913, 0.0502096, 0.0841558 };
   cblas_csymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csymm(case 1540) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csymm(case 1540) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int M = 1;
   int N = 2;
   float alpha[2] = {-0.3, 0.1};
   float beta[2] = {0, 0.1};
   float A[] = { 0.579, -0.859, 0.192, -0.737, 0.396, -0.498, 0.751, -0.379 };
   int lda = 2;
   float B[] = { 0.84, -0.755, -0.019, -0.063 };
   int ldb = 1;
   float C[] = { 0.04, 0.639, -0.876, -0.778 };
   int ldc = 1;
   float C_expected[] = { 0.115459, 0.329813, 0.288206, 0.110315 };
   cblas_csymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csymm(case 1541) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csymm(case 1541) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int M = 1;
   int N = 2;
   double alpha[2] = {0, 0};
   double beta[2] = {0, 0};
   double A[] = { 0.511, -0.486 };
   int lda = 1;
   double B[] = { 0.985, -0.923, -0.234, -0.756 };
   int ldb = 2;
   double C[] = { -0.16, 0.049, 0.618, -0.349 };
   int ldc = 2;
   double C_expected[] = { 0, 0, 0, 0 };
   cblas_zsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsymm(case 1542) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsymm(case 1542) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int M = 1;
   int N = 2;
   double alpha[2] = {0, 0};
   double beta[2] = {0, 0};
   double A[] = { 0.46, -0.816 };
   int lda = 1;
   double B[] = { 0.404, 0.113, -0.904, -0.627 };
   int ldb = 1;
   double C[] = { 0.114, 0.318, 0.636, -0.839 };
   int ldc = 1;
   double C_expected[] = { 0, 0, 0, 0 };
   cblas_zsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsymm(case 1543) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsymm(case 1543) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int M = 1;
   int N = 2;
   double alpha[2] = {-1, 0};
   double beta[2] = {-0.3, 0.1};
   double A[] = { -0.835, 0.344, 0.975, 0.634, 0.312, -0.659, -0.624, -0.175 };
   int lda = 2;
   double B[] = { -0.707, -0.846, 0.825, -0.661 };
   int ldb = 2;
   double C[] = { 0.352, -0.499, 0.267, 0.548 };
   int ldc = 2;
   double C_expected[] = { -2.160518, -0.156877, 0.648536, 0.867299 };
   cblas_zsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsymm(case 1544) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsymm(case 1544) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int M = 1;
   int N = 2;
   double alpha[2] = {-1, 0};
   double beta[2] = {-0.3, 0.1};
   double A[] = { -0.409, 0.013, -0.308, -0.317, -0.535, -0.697, -0.385, 0.119 };
   int lda = 2;
   double B[] = { 0.299, -0.233, 0.093, 0.664 };
   int ldb = 1;
   double C[] = { 0.699, 0.47, -0.347, -0.182 };
   int ldc = 1;
   double C_expected[] = { -0.550491, 0.249777, 0.559487, 0.348221 };
   cblas_zsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsymm(case 1545) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsymm(case 1545) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int M = 1;
   int N = 2;
   double alpha[2] = {1, 0};
   double beta[2] = {0, 1};
   double A[] = { -0.151, 0.635 };
   int lda = 1;
   double B[] = { 0.711, -0.869, 0.153, 0.647 };
   int ldb = 2;
   double C[] = { -0.299, 0.43, -0.307, 0.133 };
   int ldc = 2;
   double C_expected[] = { 0.014454, 0.283704, -0.566948, -0.307542 };
   cblas_zsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsymm(case 1546) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsymm(case 1546) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int M = 1;
   int N = 2;
   double alpha[2] = {1, 0};
   double beta[2] = {0, 1};
   double A[] = { 0.793, -0.543 };
   int lda = 1;
   double B[] = { 0.054, -0.045, 0.989, 0.453 };
   int ldb = 1;
   double C[] = { 0.443, -0.641, -0.809, -0.83 };
   int ldc = 1;
   double C_expected[] = { 0.659387, 0.377993, 1.860256, -0.986798 };
   cblas_zsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsymm(case 1547) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsymm(case 1547) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int M = 1;
   int N = 2;
   double alpha[2] = {1, 0};
   double beta[2] = {-1, 0};
   double A[] = { -0.432, -0.293, -0.819, 0.44, -0.818, -0.258, -0.836, 0.683 };
   int lda = 2;
   double B[] = { -0.259, -0.878, 0.161, 0.744 };
   int ldb = 2;
   double C[] = { 0.436, -0.655, -0.61, -0.875 };
   int ldc = 2;
   double C_expected[] = { -0.521112, 0.460053, -0.04741, 1.148005 };
   cblas_zsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsymm(case 1548) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsymm(case 1548) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int M = 1;
   int N = 2;
   double alpha[2] = {1, 0};
   double beta[2] = {-1, 0};
   double A[] = { -0.656, 0.378, -0.688, 0.676, 0.967, -0.804, 0.455, -0.425 };
   int lda = 2;
   double B[] = { 0.791, -0.947, -0.945, -0.444 };
   int ldb = 1;
   double C[] = { 0.014, -0.814, -0.091, -0.417 };
   int ldc = 1;
   double C_expected[] = { 0.775374, 1.400882, -0.431711, 1.802857 };
   cblas_zsymm(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsymm(case 1549) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsymm(case 1549) imag");
     };
   };
  };


}
