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
   int N = 2;
   int K = 1;
   float alpha = 1;
   float beta = 1;
   float A[] = { 0.934, 0.664, 0.426, 0.263 };
   int lda = 1;
   float C[] = { 0.251, -0.97, 0.76, -0.349, 0.152, -0.899, -0.17, 0.707 };
   int ldc = 2;
   float C_expected[] = { 1.56425, 0, 1.33252, -0.311778, 0.152, -0.899, 0.080645, 0 };
   cblas_cherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cherk(case 1606) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cherk(case 1606) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 111;
   int N = 2;
   int K = 1;
   float alpha = 1;
   float beta = 1;
   float A[] = { 0.16, 0.464, -0.623, 0.776 };
   int lda = 2;
   float C[] = { 0.771, -0.449, 0.776, 0.112, -0.134, 0.317, 0.547, -0.551 };
   int ldc = 2;
   float C_expected[] = { 1.0119, 0, 0.776, 0.112, 0.126384, -0.096232, 1.5373, 0 };
   cblas_cherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cherk(case 1607) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cherk(case 1607) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 113;
   int N = 2;
   int K = 1;
   float alpha = 0.1;
   float beta = 1;
   float A[] = { 0.787, 0.057, -0.49, 0.47 };
   int lda = 2;
   float C[] = { -0.758, 0.912, 0.992, -0.356, 0.584, 0.806, 0.965, 0.674 };
   int ldc = 2;
   float C_expected[] = { -0.695738, 0, 0.956116, -0.316218, 0.584, 0.806, 1.0111, 0 };
   cblas_cherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cherk(case 1608) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cherk(case 1608) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 113;
   int N = 2;
   int K = 1;
   float alpha = 0.1;
   float beta = 1;
   float A[] = { 0.961, -0.384, 0.165, 0.395 };
   int lda = 1;
   float C[] = { -0.186, 0.404, -0.873, 0.09, -0.451, -0.972, -0.203, -0.304 };
   int ldc = 2;
   float C_expected[] = { -0.0789023, 0, -0.873, 0.09, -0.450312, -0.927704, -0.184675, 0 };
   cblas_cherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cherk(case 1609) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cherk(case 1609) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 111;
   int N = 2;
   int K = 1;
   float alpha = 0;
   float beta = -0.3;
   float A[] = { 0.04, 0.608, 0.21, -0.44 };
   int lda = 1;
   float C[] = { 0.285, -0.943, 0.581, -0.56, 0.112, 0.529, 0.16, -0.913 };
   int ldc = 2;
   float C_expected[] = { -0.0855, 0, 0.581, -0.56, -0.0336, -0.1587, -0.048, 0 };
   cblas_cherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cherk(case 1610) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cherk(case 1610) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 111;
   int N = 2;
   int K = 1;
   float alpha = 0;
   float beta = -0.3;
   float A[] = { -0.984, -0.398, -0.379, 0.919 };
   int lda = 2;
   float C[] = { -0.44, -0.087, 0.156, -0.945, -0.943, -0.355, 0.577, 0.053 };
   int ldc = 2;
   float C_expected[] = { 0.132, 0, -0.0468, 0.2835, -0.943, -0.355, -0.1731, 0 };
   cblas_cherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cherk(case 1611) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cherk(case 1611) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 113;
   int N = 2;
   int K = 1;
   float alpha = 1;
   float beta = -1;
   float A[] = { 0.269, -0.428, -0.029, 0.964 };
   int lda = 2;
   float C[] = { 0.473, -0.932, -0.689, -0.072, -0.952, -0.862, 0.001, 0.282 };
   int ldc = 2;
   float C_expected[] = { -0.217455, 0, -0.689, -0.072, 0.531607, 0.615096, 0.929137, 0 };
   cblas_cherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cherk(case 1612) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cherk(case 1612) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 113;
   int N = 2;
   int K = 1;
   float alpha = 1;
   float beta = -1;
   float A[] = { -0.303, -0.037, -0.411, -0.243 };
   int lda = 1;
   float C[] = { 0.652, -0.227, -0.849, 0.87, -0.051, -0.535, 0.418, -0.681 };
   int ldc = 2;
   float C_expected[] = { -0.558822, 0, 0.982524, -0.928422, -0.051, -0.535, -0.19003, 0 };
   cblas_cherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cherk(case 1613) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cherk(case 1613) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 111;
   int N = 2;
   int K = 1;
   double alpha = 0.1;
   double beta = 1;
   double A[] = { -0.384, -0.851, 0.518, 0.492 };
   int lda = 1;
   double C[] = { -0.117, -0.194, -0.915, 0.069, 0.445, 0.089, 0.213, -0.889 };
   int ldc = 2;
   double C_expected[] = { -0.0298343, 0, -0.9767604, 0.043811, 0.445, 0.089, 0.2640388, 0 };
   cblas_zherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zherk(case 1614) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zherk(case 1614) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 111;
   int N = 2;
   int K = 1;
   double alpha = 0.1;
   double beta = 1;
   double A[] = { 0.13, 0.236, 0.788, 0.629 };
   int lda = 2;
   double C[] = { 0.021, -0.376, -0.804, 0.689, -0.912, 0.21, -0.581, 0.406 };
   int ldc = 2;
   double C_expected[] = { 0.0282596, 0, -0.804, 0.689, -0.8869116, 0.2204198, -0.4793415, 0 };
   cblas_zherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zherk(case 1615) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zherk(case 1615) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 113;
   int N = 2;
   int K = 1;
   double alpha = 0;
   double beta = 1;
   double A[] = { 0.593, 0.846, -0.144, 0.128 };
   int lda = 2;
   double C[] = { 0.02, 0.313, 0.222, 0.301, 0.412, -0.645, -0.411, -0.02 };
   int ldc = 2;
   double C_expected[] = { 0.02, 0.313, 0.222, 0.301, 0.412, -0.645, -0.411, -0.02 };
   cblas_zherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zherk(case 1616) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zherk(case 1616) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 113;
   int N = 2;
   int K = 1;
   double alpha = 0;
   double beta = 1;
   double A[] = { 0.857, 0.994, -0.933, 0.069 };
   int lda = 1;
   double C[] = { 0.253, -0.521, 0.937, -0.73, 0.24, 0.177, -0.27, -0.225 };
   int ldc = 2;
   double C_expected[] = { 0.253, -0.521, 0.937, -0.73, 0.24, 0.177, -0.27, -0.225 };
   cblas_zherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zherk(case 1617) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zherk(case 1617) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 111;
   int N = 2;
   int K = 1;
   double alpha = 0;
   double beta = 1;
   double A[] = { -0.343, -0.433, -0.381, -0.087 };
   int lda = 1;
   double C[] = { -0.695, 0.911, 0.719, -0.074, -0.621, -0.256, 0.216, -0.889 };
   int ldc = 2;
   double C_expected[] = { -0.695, 0.911, 0.719, -0.074, -0.621, -0.256, 0.216, -0.889 };
   cblas_zherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zherk(case 1618) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zherk(case 1618) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 111;
   int N = 2;
   int K = 1;
   double alpha = 0;
   double beta = 1;
   double A[] = { -0.887, 0.557, -0.43, 0.912 };
   int lda = 2;
   double C[] = { -0.083, 0.219, 0.417, 0.817, -0.294, -0.683, -0.633, 0.831 };
   int ldc = 2;
   double C_expected[] = { -0.083, 0.219, 0.417, 0.817, -0.294, -0.683, -0.633, 0.831 };
   cblas_zherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zherk(case 1619) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zherk(case 1619) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 113;
   int N = 2;
   int K = 1;
   double alpha = -0.3;
   double beta = -1;
   double A[] = { -0.531, 0.187, -0.777, -0.329 };
   int lda = 2;
   double C[] = { -0.173, 0.833, 0.155, -0.52, -0.99, 0.28, 0.455, 0.481 };
   int ldc = 2;
   double C_expected[] = { 0.077921, 0, 0.155, -0.52, 0.8846808, -0.1840006, -0.668591, 0 };
   cblas_zherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zherk(case 1620) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zherk(case 1620) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 113;
   int N = 2;
   int K = 1;
   double alpha = -0.3;
   double beta = -1;
   double A[] = { -0.287, 0.068, 0.917, -0.449 };
   int lda = 1;
   double C[] = { -0.248, -0.608, -0.124, -0.718, -0.037, -0.115, 0.998, -0.551 };
   int ldc = 2;
   double C_expected[] = { 0.2219021, 0, 0.2121133, 0.7379521, -0.037, -0.115, -1.310747, 0 };
   cblas_zherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zherk(case 1621) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zherk(case 1621) imag");
     };
   };
  };


}
