#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_tbmv () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { 0.439, -0.484, -0.952, -0.508, 0.381, -0.889, -0.192, -0.279, -0.155 };
   float X[] = { -0.089, -0.688, -0.203 };
   int incX = -1;
   float x_expected[] = { 0.017088, 0.315595, 0.243875 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 894)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { 0.439, -0.484, -0.952, -0.508, 0.381, -0.889, -0.192, -0.279, -0.155 };
   float X[] = { -0.089, -0.688, -0.203 };
   int incX = -1;
   float x_expected[] = { -0.089, -0.721909, 0.129992 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 895)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { 0.439, -0.484, -0.952, -0.508, 0.381, -0.889, -0.192, -0.279, -0.155 };
   float X[] = { -0.089, -0.688, -0.203 };
   int incX = -1;
   float x_expected[] = { 0.156927, -0.159004, 0.098252 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 896)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { 0.439, -0.484, -0.952, -0.508, 0.381, -0.889, -0.192, -0.279, -0.155 };
   float X[] = { -0.089, -0.688, -0.203 };
   int incX = -1;
   float x_expected[] = { 0.043096, -0.584876, -0.203 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 897)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { 0.439, -0.484, -0.952, -0.508, 0.381, -0.889, -0.192, -0.279, -0.155 };
   float X[] = { -0.089, -0.688, -0.203 };
   int incX = -1;
   float x_expected[] = { 0.024831, -0.24504, 0.447756 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 898)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { 0.439, -0.484, -0.952, -0.508, 0.381, -0.889, -0.192, -0.279, -0.155 };
   float X[] = { -0.089, -0.688, -0.203 };
   int incX = -1;
   float x_expected[] = { -0.089, -0.670912, 0.146504 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 899)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { 0.439, -0.484, -0.952, -0.508, 0.381, -0.889, -0.192, -0.279, -0.155 };
   float X[] = { -0.089, -0.688, -0.203 };
   int incX = -1;
   float x_expected[] = { -0.24504, 0.447756, -0.089117 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 900)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { 0.439, -0.484, -0.952, -0.508, 0.381, -0.889, -0.192, -0.279, -0.155 };
   float X[] = { -0.089, -0.688, -0.203 };
   int incX = -1;
   float x_expected[] = { -0.351128, -0.589748, -0.203 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 901)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { 0.94, -0.091, 0.984, -0.276, -0.342, -0.484, -0.665, -0.2, 0.349 };
   float X[] = { 0.023, -0.501, -0.562 };
   int incX = -1;
   float x_expected[] = { 0.156047, 0.189418, -0.52828 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 902)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { 0.94, -0.091, 0.984, -0.276, -0.342, -0.484, -0.665, -0.2, 0.349 };
   float X[] = { 0.023, -0.501, -0.562 };
   int incX = -1;
   float x_expected[] = { 0.194342, -0.449858, -0.562 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 903)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { 0.94, -0.091, 0.984, -0.276, -0.342, -0.484, -0.665, -0.2, 0.349 };
   float X[] = { 0.023, -0.501, -0.562 };
   int incX = -1;
   float x_expected[] = { -0.0046, 0.156047, 0.189418 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 904)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { 0.94, -0.091, 0.984, -0.276, -0.342, -0.484, -0.665, -0.2, 0.349 };
   float X[] = { 0.023, -0.501, -0.562 };
   int incX = -1;
   float x_expected[] = { 0.023, -0.516295, -0.423724 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 905)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { 0.94, -0.091, 0.984, -0.276, -0.342, -0.484, -0.665, -0.2, 0.349 };
   float X[] = { 0.023, -0.501, -0.562 };
   int incX = -1;
   float x_expected[] = { 0.328565, 0.326454, 0.051142 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 906)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { 0.94, -0.091, 0.984, -0.276, -0.342, -0.484, -0.665, -0.2, 0.349 };
   float X[] = { 0.023, -0.501, -0.562 };
   int incX = -1;
   float x_expected[] = { 0.356165, -0.345888, -0.562 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 907)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { 0.94, -0.091, 0.984, -0.276, -0.342, -0.484, -0.665, -0.2, 0.349 };
   float X[] = { 0.023, -0.501, -0.562 };
   int incX = -1;
   float x_expected[] = { -0.015295, 0.13041, -0.482689 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 908)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { 0.94, -0.091, 0.984, -0.276, -0.342, -0.484, -0.665, -0.2, 0.349 };
   float X[] = { 0.023, -0.501, -0.562 };
   int incX = -1;
   float x_expected[] = { 0.023, -0.508866, -0.516409 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 909)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { 0.566, 0.955, -0.086, -0.856, 0.177, 0.974, -0.652, -0.884, 0.77 };
   double X[] = { -0.77, -0.818, 0.337 };
   int incX = -1;
   double x_expected[] = { 0.50204, 0.563918, -0.590448 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 910)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { 0.566, 0.955, -0.086, -0.856, 0.177, 0.974, -0.652, -0.884, 0.77 };
   double X[] = { -0.77, -0.818, 0.337 };
   int incX = -1;
   double x_expected[] = { -0.77, -0.95429, -0.44419 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 911)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { 0.566, 0.955, -0.086, -0.856, 0.177, 0.974, -0.652, -0.884, 0.77 };
   double X[] = { -0.77, -0.818, 0.337 };
   int incX = -1;
   double x_expected[] = { 1.214016, -0.433258, 0.321835 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 912)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { 0.566, 0.955, -0.086, -0.856, 0.177, 0.974, -0.652, -0.884, 0.77 };
   double X[] = { -0.77, -0.818, 0.337 };
   int incX = -1;
   double x_expected[] = { -0.236664, -1.106472, 0.337 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 913)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { 0.566, 0.955, -0.086, -0.856, 0.177, 0.974, -0.652, -0.884, 0.77 };
   double X[] = { -0.77, -0.818, 0.337 };
   int incX = -1;
   double x_expected[] = { 0.68068, 0.357254, 1.022043 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 914)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { 0.566, 0.955, -0.086, -0.856, 0.177, 0.974, -0.652, -0.884, 0.77 };
   double X[] = { -0.77, -0.818, 0.337 };
   int incX = -1;
   double x_expected[] = { -0.77, -0.31596, 1.037208 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 915)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { 0.566, 0.955, -0.086, -0.856, 0.177, 0.974, -0.652, -0.884, 0.77 };
   double X[] = { -0.77, -0.818, 0.337 };
   int incX = -1;
   double x_expected[] = { 0.357254, 1.022043, 0.190742 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 916)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { 0.566, 0.955, -0.086, -0.856, 0.177, 0.974, -0.652, -0.884, 0.77 };
   double X[] = { -0.77, -0.818, 0.337 };
   int incX = -1;
   double x_expected[] = { -0.914786, -0.496165, 0.337 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 917)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { 0.188, 0.6, -0.743, -0.803, 0.449, -0.681, -0.464, -0.029, 0.553 };
   double X[] = { -0.851, 0.481, 0.155 };
   int incX = -1;
   double x_expected[] = { 0.610833, -0.293243, 0.02914 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 918)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { 0.188, 0.6, -0.743, -0.803, 0.449, -0.681, -0.464, -0.029, 0.553 };
   double X[] = { -0.851, 0.481, 0.155 };
   int incX = -1;
   double x_expected[] = { -0.635031, 0.574, 0.155 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 919)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { 0.188, 0.6, -0.743, -0.803, 0.449, -0.681, -0.464, -0.029, 0.553 };
   double X[] = { -0.851, 0.481, 0.155 };
   int incX = -1;
   double x_expected[] = { 0.024679, 0.610833, -0.293243 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 920)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { 0.188, 0.6, -0.743, -0.803, 0.449, -0.681, -0.464, -0.029, 0.553 };
   double X[] = { -0.851, 0.481, 0.155 };
   int incX = -1;
   double x_expected[] = { -0.851, 0.875864, -0.231243 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 921)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { 0.188, 0.6, -0.743, -0.803, 0.449, -0.681, -0.464, -0.029, 0.553 };
   double X[] = { -0.851, 0.481, 0.155 };
   int incX = -1;
   double x_expected[] = { -0.198505, 0.091504, 0.093 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 922)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { 0.188, 0.6, -0.743, -0.803, 0.449, -0.681, -0.464, -0.029, 0.553 };
   double X[] = { -0.851, 0.481, 0.155 };
   int incX = -1;
   double x_expected[] = { -1.074184, 0.356535, 0.155 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 923)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { 0.188, 0.6, -0.743, -0.803, 0.449, -0.681, -0.464, -0.029, 0.553 };
   double X[] = { -0.851, 0.481, 0.155 };
   int incX = -1;
   double x_expected[] = { 0.394864, -0.768342, 0.31774 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 924)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { 0.188, 0.6, -0.743, -0.803, 0.449, -0.681, -0.464, -0.029, 0.553 };
   double X[] = { -0.851, 0.481, 0.155 };
   int incX = -1;
   double x_expected[] = { -0.851, 0.098901, 0.4436 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 925)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { 0.824, -0.45, -0.987, 0.758, 0.42, -0.357, 0.147, -0.191, 0.88, 0.63, 0.155, -0.573, 0.224, 0.146, 0.501, -0.889, 0.456, 0.796 };
   float X[] = { -0.46, 0.069, 0.308, -0.003, -0.236, 0.605 };
   int incX = -1;
   float x_expected[] = { -0.113114, -0.051704, -0.403567, -0.288349, -0.223936, 0.841145 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 926) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 926) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { 0.824, -0.45, -0.987, 0.758, 0.42, -0.357, 0.147, -0.191, 0.88, 0.63, 0.155, -0.573, 0.224, 0.146, 0.501, -0.889, 0.456, 0.796 };
   float X[] = { -0.46, 0.069, 0.308, -0.003, -0.236, 0.605 };
   int incX = -1;
   float x_expected[] = { -0.46, 0.069, -0.14027, -0.23208, -0.537722, 0.841425 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 927) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 927) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { 0.824, -0.45, -0.987, 0.758, 0.42, -0.357, 0.147, -0.191, 0.88, 0.63, 0.155, -0.573, 0.224, 0.146, 0.501, -0.889, 0.456, 0.796 };
   float X[] = { -0.46, 0.069, 0.308, -0.003, -0.236, 0.605 };
   int incX = -1;
   float x_expected[] = { -0.099689, 0.487805, 0.353793, 0.325411, -0.225658, -0.776023 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 928) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 928) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { 0.824, -0.45, -0.987, 0.758, 0.42, -0.357, 0.147, -0.191, 0.88, 0.63, 0.155, -0.573, 0.224, 0.146, 0.501, -0.889, 0.456, 0.796 };
   float X[] = { -0.46, 0.069, 0.308, -0.003, -0.236, 0.605 };
   int incX = -1;
   float x_expected[] = { -0.39057, 0.113296, 0.388863, 0.131011, -0.236, 0.605 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 929) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 929) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { 0.824, -0.45, -0.987, 0.758, 0.42, -0.357, 0.147, -0.191, 0.88, 0.63, 0.155, -0.573, 0.224, 0.146, 0.501, -0.889, 0.456, 0.796 };
   float X[] = { -0.46, 0.069, 0.308, -0.003, -0.236, 0.605 };
   int incX = -1;
   float x_expected[] = { -0.169119, 0.443509, 0.159816, 0.139696, -0.180955, -0.835292 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 930) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 930) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { 0.824, -0.45, -0.987, 0.758, 0.42, -0.357, 0.147, -0.191, 0.88, 0.63, 0.155, -0.573, 0.224, 0.146, 0.501, -0.889, 0.456, 0.796 };
   float X[] = { -0.46, 0.069, 0.308, -0.003, -0.236, 0.605 };
   int incX = -1;
   float x_expected[] = { -0.46, 0.069, 0.194886, -0.054704, -0.191297, 0.545731 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 931) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 931) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { 0.824, -0.45, -0.987, 0.758, 0.42, -0.357, 0.147, -0.191, 0.88, 0.63, 0.155, -0.573, 0.224, 0.146, 0.501, -0.889, 0.456, 0.796 };
   float X[] = { -0.46, 0.069, 0.308, -0.003, -0.236, 0.605 };
   int incX = -1;
   float x_expected[] = { 0.159816, 0.139696, -0.180955, -0.835292, 0.077786, 0.60472 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 932) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 932) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { 0.824, -0.45, -0.987, 0.758, 0.42, -0.357, 0.147, -0.191, 0.88, 0.63, 0.155, -0.573, 0.224, 0.146, 0.501, -0.889, 0.456, 0.796 };
   float X[] = { -0.46, 0.069, 0.308, -0.003, -0.236, 0.605 };
   int incX = -1;
   float x_expected[] = { -0.18707, 0.2604, 0.082342, -0.779023, -0.236, 0.605 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 933) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 933) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { -0.814, 0.043, -0.755, -0.094, 0.876, 0.257, 0.406, 0.491, -0.27, -0.787, 0.545, 0.732, -0.512, -0.085, 0.234, 0.001, -0.225, -0.002 };
   float X[] = { 0.411, -0.308, -0.912, 0.811, -0.894, -0.116 };
   int incX = -1;
   float x_expected[] = { 0.647885, 0.621535, -0.104407, 0.05309, 0.732704, 0.055982 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 934) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 934) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { -0.814, 0.043, -0.755, -0.094, 0.876, 0.257, 0.406, 0.491, -0.27, -0.787, 0.545, 0.732, -0.512, -0.085, 0.234, 0.001, -0.225, -0.002 };
   float X[] = { 0.411, -0.308, -0.912, 0.811, -0.894, -0.116 };
   int incX = -1;
   float x_expected[] = { 1.2955, 0.190774, -0.247934, 0.982616, -0.894, -0.116 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 935) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 935) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { -0.814, 0.043, -0.755, -0.094, 0.876, 0.257, 0.406, 0.491, -0.27, -0.787, 0.545, 0.732, -0.512, -0.085, 0.234, 0.001, -0.225, -0.002 };
   float X[] = { 0.411, -0.308, -0.912, 0.811, -0.894, -0.116 };
   int incX = -1;
   float x_expected[] = { 0.096482, -0.071661, 0.647885, 0.621535, -0.104407, 0.05309 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 936) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 936) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { -0.814, 0.043, -0.755, -0.094, 0.876, 0.257, 0.406, 0.491, -0.27, -0.787, 0.545, 0.732, -0.512, -0.085, 0.234, 0.001, -0.225, -0.002 };
   float X[] = { 0.411, -0.308, -0.912, 0.811, -0.894, -0.116 };
   int incX = -1;
   float x_expected[] = { 0.411, -0.308, -1.14861, 0.933761, -1.66247, -0.234526 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 937) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 937) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { -0.814, 0.043, -0.755, -0.094, 0.876, 0.257, 0.406, 0.491, -0.27, -0.787, 0.545, 0.732, -0.512, -0.085, 0.234, 0.001, -0.225, -0.002 };
   float X[] = { 0.411, -0.308, -0.912, 0.811, -0.894, -0.116 };
   int incX = -1;
   float x_expected[] = { 0.632361, -0.409373, 0.578489, 0.012724, 0.664066, 0.171616 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 938) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 938) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { -0.814, 0.043, -0.755, -0.094, 0.876, 0.257, 0.406, 0.491, -0.27, -0.787, 0.545, 0.732, -0.512, -0.085, 0.234, 0.001, -0.225, -0.002 };
   float X[] = { 0.411, -0.308, -0.912, 0.811, -0.894, -0.116 };
   int incX = -1;
   float x_expected[] = { 0.946879, -0.645712, -1.21801, 0.32495, -0.894, -0.116 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 939) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 939) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { -0.814, 0.043, -0.755, -0.094, 0.876, 0.257, 0.406, 0.491, -0.27, -0.787, 0.545, 0.732, -0.512, -0.085, 0.234, 0.001, -0.225, -0.002 };
   float X[] = { 0.411, -0.308, -0.912, 0.811, -0.894, -0.116 };
   int incX = -1;
   float x_expected[] = { -0.236612, 0.122761, -1.12184, -0.358823, 1.4975, -0.470595 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 940) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 940) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { -0.814, 0.043, -0.755, -0.094, 0.876, 0.257, 0.406, 0.491, -0.27, -0.787, 0.545, 0.732, -0.512, -0.085, 0.234, 0.001, -0.225, -0.002 };
   float X[] = { 0.411, -0.308, -0.912, 0.811, -0.894, -0.116 };
   int incX = -1;
   float x_expected[] = { 0.411, -0.308, -1.26537, 0.570703, -0.129206, -0.642577 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 941) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 941) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 121;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { -0.675, 0.047, 0.695, 0.724, -0.438, 0.991, -0.188, -0.06, -0.093, 0.302, 0.842, -0.753, 0.465, -0.972, -0.058, 0.988, 0.093, 0.164 };
   float X[] = { 0.065, -0.082, -0.746, 0.775, 0.694, -0.954 };
   int incX = -1;
   float x_expected[] = { 0.413357, 0.178267, -0.114618, -1.35595, -0.513288, 0.611332 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 942) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 942) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 121;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { -0.675, 0.047, 0.695, 0.724, -0.438, 0.991, -0.188, -0.06, -0.093, 0.302, 0.842, -0.753, 0.465, -0.972, -0.058, 0.988, 0.093, 0.164 };
   float X[] = { 0.065, -0.082, -0.746, 0.775, 0.694, -0.954 };
   int incX = -1;
   float x_expected[] = { 0.368428, 0.071217, -0.954366, -0.390486, 0.694, -0.954 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 943) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 943) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 122;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { -0.675, 0.047, 0.695, 0.724, -0.438, 0.991, -0.188, -0.06, -0.093, 0.302, 0.842, -0.753, 0.465, -0.972, -0.058, 0.988, 0.093, 0.164 };
   float X[] = { 0.065, -0.082, -0.746, 0.775, 0.694, -0.954 };
   int incX = -1;
   float x_expected[] = { -0.084786, -0.059464, 0.413357, 0.178267, -0.114618, -1.35595 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 944) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 944) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 122;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { -0.675, 0.047, 0.695, 0.724, -0.438, 0.991, -0.188, -0.06, -0.093, 0.302, 0.842, -0.753, 0.465, -0.972, -0.058, 0.988, 0.093, 0.164 };
   float X[] = { 0.065, -0.082, -0.746, 0.775, 0.694, -0.954 };
   int incX = -1;
   float x_expected[] = { 0.065, -0.082, -0.636071, 0.80005, 0.787748, -1.14446 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 945) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 945) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 121;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { -0.675, 0.047, 0.695, 0.724, -0.438, 0.991, -0.188, -0.06, -0.093, 0.302, 0.842, -0.753, 0.465, -0.972, -0.058, 0.988, 0.093, 0.164 };
   float X[] = { 0.065, -0.082, -0.746, 0.775, 0.694, -0.954 };
   int incX = -1;
   float x_expected[] = { -1.18498, -0.424201, 0.230196, 0.374209, -0.208366, -1.16549 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 946) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 946) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 121;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { -0.675, 0.047, 0.695, 0.724, -0.438, 0.991, -0.188, -0.06, -0.093, 0.302, 0.842, -0.753, 0.465, -0.972, -0.058, 0.988, 0.093, 0.164 };
   float X[] = { 0.065, -0.082, -0.746, 0.775, 0.694, -0.954 };
   int incX = -1;
   float x_expected[] = { -1.03519, -0.446737, -0.819232, 0.995992, 0.694, -0.954 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 947) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 947) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 122;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { -0.675, 0.047, 0.695, 0.724, -0.438, 0.991, -0.188, -0.06, -0.093, 0.302, 0.842, -0.753, 0.465, -0.972, -0.058, 0.988, 0.093, 0.164 };
   float X[] = { 0.065, -0.082, -0.746, 0.775, 0.694, -0.954 };
   int incX = -1;
   float x_expected[] = { 0.109929, 0.02505, 0.062939, -0.202464, -0.470658, 1.69006 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 948) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 948) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 122;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { -0.675, 0.047, 0.695, 0.724, -0.438, 0.991, -0.188, -0.06, -0.093, 0.302, 0.842, -0.753, 0.465, -0.972, -0.058, 0.988, 0.093, 0.164 };
   float X[] = { 0.065, -0.082, -0.746, 0.775, 0.694, -0.954 };
   int incX = -1;
   float x_expected[] = { 0.065, -0.082, -0.776809, 0.762996, 0.73663, 0.124729 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 949) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 949) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { -0.212, 0.612, 0.189, -0.046, -0.124, 0.82, 0.753, 0.727, 0.331, 0.116, 0.504, -0.673, -0.888, -0.277, -0.361, -0.909, 0.982, -0.124 };
   double X[] = { 0.064, 0.169, -0.81, -0.779, -0.641, -0.103 };
   int incX = -1;
   double x_expected[] = { -0.010019, -0.1678, -0.042017, -1.112094, 0.010004, -0.480427 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 950) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 950) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { -0.212, 0.612, 0.189, -0.046, -0.124, 0.82, 0.753, 0.727, 0.331, 0.116, 0.504, -0.673, -0.888, -0.277, -0.361, -0.909, 0.982, -0.124 };
   double X[] = { 0.064, 0.169, -0.81, -0.779, -0.641, -0.103 };
   int incX = -1;
   double x_expected[] = { 0.064, 0.169, -0.80842, -0.715637, -0.829924, -0.212971 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 951) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 951) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { -0.212, 0.612, 0.189, -0.046, -0.124, 0.82, 0.753, 0.727, 0.331, 0.116, 0.504, -0.673, -0.888, -0.277, -0.361, -0.909, 0.982, -0.124 };
   double X[] = { 0.064, 0.169, -0.81, -0.779, -0.641, -0.103 };
   int incX = -1;
   double x_expected[] = { 0.634014, 0.796937, -0.585538, -0.895375, -0.125887, 0.010019 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 952) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 952) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { -0.212, 0.612, 0.189, -0.046, -0.124, 0.82, 0.753, 0.727, 0.331, 0.116, 0.504, -0.673, -0.888, -0.277, -0.361, -0.909, 0.982, -0.124 };
   double X[] = { 0.064, 0.169, -0.81, -0.779, -0.641, -0.103 };
   int incX = -1;
   double x_expected[] = { 0.567497, 1.085122, -1.217792, -1.322566, -0.641, -0.103 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 953) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 953) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { -0.212, 0.612, 0.189, -0.046, -0.124, 0.82, 0.753, 0.727, 0.331, 0.116, 0.504, -0.673, -0.888, -0.277, -0.361, -0.909, 0.982, -0.124 };
   double X[] = { 0.064, 0.169, -0.81, -0.779, -0.641, -0.103 };
   int incX = -1;
   double x_expected[] = { 0.130517, -0.119185, -0.187765, -0.519609, -0.169484, -1.165438 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 954) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 954) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { -0.212, 0.612, 0.189, -0.046, -0.124, 0.82, 0.753, 0.727, 0.331, 0.116, 0.504, -0.673, -0.888, -0.277, -0.361, -0.909, 0.982, -0.124 };
   double X[] = { 0.064, 0.169, -0.81, -0.779, -0.641, -0.103 };
   int incX = -1;
   double x_expected[] = { 0.064, 0.169, -0.820019, -0.9468, -0.684597, -1.278457 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 955) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 955) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { -0.212, 0.612, 0.189, -0.046, -0.124, 0.82, 0.753, 0.727, 0.331, 0.116, 0.504, -0.673, -0.888, -0.277, -0.361, -0.909, 0.982, -0.124 };
   double X[] = { 0.064, 0.169, -0.81, -0.779, -0.641, -0.103 };
   int incX = -1;
   double x_expected[] = { -0.187765, -0.519609, -0.169484, -1.165438, 0.198928, -0.370456 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 956) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 956) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { -0.212, 0.612, 0.189, -0.046, -0.124, 0.82, 0.753, 0.727, 0.331, 0.116, 0.504, -0.673, -0.888, -0.277, -0.361, -0.909, 0.982, -0.124 };
   double X[] = { 0.064, 0.169, -0.81, -0.779, -0.641, -0.103 };
   int incX = -1;
   double x_expected[] = { -0.113746, -0.182809, -0.935887, -0.768981, -0.641, -0.103 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 957) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 957) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { -0.374, -0.308, 0.792, 0.884, -0.794, -0.055, -0.281, 0.527, 0.246, 0.762, 0.853, 0.891, -0.231, 0.384, 0.373, -0.717, -0.957, -0.338 };
   double X[] = { 0.884, 0.636, 0.921, 0.282, -0.649, 0.188 };
   int incX = -1;
   double x_expected[] = { -0.436746, 0.963714, -1.087615, -0.018695, 0.30063, 0.12958 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 958) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 958) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { -0.374, -0.308, 0.792, 0.884, -0.794, -0.055, -0.281, 0.527, 0.246, 0.762, 0.853, 0.891, -0.231, 0.384, 0.373, -0.717, -0.957, -0.338 };
   double X[] = { 0.884, 0.636, 0.921, 0.282, -0.649, 0.188 };
   int incX = -1;
   double x_expected[] = { 0.895682, 1.407174, 0.2408, -0.14282, -0.649, 0.188 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 959) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 959) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { -0.374, -0.308, 0.792, 0.884, -0.794, -0.055, -0.281, 0.527, 0.246, 0.762, 0.853, 0.891, -0.231, 0.384, 0.373, -0.717, -0.957, -0.338 };
   double X[] = { 0.884, 0.636, 0.921, 0.282, -0.649, 0.188 };
   int incX = -1;
   double x_expected[] = { 0.785744, -0.3966, -0.436746, 0.963714, -1.087615, -0.018695 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 960) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 960) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { -0.374, -0.308, 0.792, 0.884, -0.794, -0.055, -0.281, 0.527, 0.246, 0.762, 0.853, 0.891, -0.231, 0.384, 0.373, -0.717, -0.957, -0.338 };
   double X[] = { 0.884, 0.636, 0.921, 0.282, -0.649, 0.188 };
   int incX = -1;
   double x_expected[] = { 0.884, 0.636, 0.472572, 0.47454, -1.056415, 0.594125 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 961) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 961) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { -0.374, -0.308, 0.792, 0.884, -0.794, -0.055, -0.281, 0.527, 0.246, 0.762, 0.853, 0.891, -0.231, 0.384, 0.373, -0.717, -0.957, -0.338 };
   double X[] = { 0.884, 0.636, 0.921, 0.282, -0.649, 0.188 };
   int incX = -1;
   double x_expected[] = { 0.464705, -0.108078, 0.094975, 0.376323, -0.6802, -0.42482 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 962) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 962) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { -0.374, -0.308, 0.792, 0.884, -0.794, -0.055, -0.281, 0.527, 0.246, 0.762, 0.853, 0.891, -0.231, 0.384, 0.373, -0.717, -0.957, -0.338 };
   double X[] = { 0.884, 0.636, 0.921, 0.282, -0.649, 0.188 };
   int incX = -1;
   double x_expected[] = { 0.562961, 0.924522, 1.004293, -0.112851, -0.649, 0.188 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 963) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 963) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { -0.374, -0.308, 0.792, 0.884, -0.794, -0.055, -0.281, 0.527, 0.246, 0.762, 0.853, 0.891, -0.231, 0.384, 0.373, -0.717, -0.957, -0.338 };
   double X[] = { 0.884, 0.636, 0.921, 0.282, -0.649, 0.188 };
   int incX = -1;
   double x_expected[] = { -0.448428, 0.19254, -0.674583, 1.236189, 0.780774, 1.167088 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 964) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 964) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { -0.374, -0.308, 0.792, 0.884, -0.794, -0.055, -0.281, 0.527, 0.246, 0.762, 0.853, 0.891, -0.231, 0.384, 0.373, -0.717, -0.957, -0.338 };
   double X[] = { 0.884, 0.636, 0.921, 0.282, -0.649, 0.188 };
   int incX = -1;
   double x_expected[] = { 0.884, 0.636, 0.653832, 1.112064, -0.168856, 1.225508 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 965) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 965) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 121;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { 0.002, 0.95, -0.363, 0.084, -0.646, 0.816, -0.407, 0.099, -0.02, -0.906, -0.874, 0.191, -0.328, -0.968, 0.79, 0.826, -0.795, 0.277 };
   double X[] = { -0.591, -0.084, 0.707, 0.945, 0.326, 0.428 };
   int incX = -1;
   double x_expected[] = { -0.59515, 0.077106, -0.27658, -0.637356, 0.407252, -0.308844 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 966) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 966) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 121;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { 0.002, 0.95, -0.363, 0.084, -0.646, 0.816, -0.407, 0.099, -0.02, -0.906, -0.874, 0.191, -0.328, -0.968, 0.79, 0.826, -0.795, 0.277 };
   double X[] = { -0.591, -0.084, 0.707, 0.945, 0.326, 0.428 };
   int incX = -1;
   double x_expected[] = { -1.46131, 0.537642, 0.624614, 0.762252, 0.326, 0.428 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 967) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 967) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 122;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { 0.002, 0.95, -0.363, 0.084, -0.646, 0.816, -0.407, 0.099, -0.02, -0.906, -0.874, 0.191, -0.328, -0.968, 0.79, 0.826, -0.795, 0.277 };
   double X[] = { -0.591, -0.084, 0.707, 0.945, 0.326, 0.428 };
   int incX = -1;
   double x_expected[] = { -0.536274, 0.421806, -0.59515, 0.077106, -0.27658, -0.637356 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 968) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 968) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 122;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { 0.002, 0.95, -0.363, 0.084, -0.646, 0.816, -0.407, 0.099, -0.02, -0.906, -0.874, 0.191, -0.328, -0.968, 0.79, 0.826, -0.795, 0.277 };
   double X[] = { -0.591, -0.084, 0.707, 0.945, 0.326, 0.428 };
   int incX = -1;
   double x_expected[] = { -0.591, -0.084, 0.98216, 0.400464, 0.131806, -0.026608 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 969) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 969) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 121;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { 0.002, 0.95, -0.363, 0.084, -0.646, 0.816, -0.407, 0.099, -0.02, -0.906, -0.874, 0.191, -0.328, -0.968, 0.79, 0.826, -0.795, 0.277 };
   double X[] = { -0.591, -0.084, 0.707, 0.945, 0.326, 0.428 };
   int incX = -1;
   double x_expected[] = { -1.68293, 0.796222, -0.96062, 0.415172, -0.082386, -0.182748 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 970) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 970) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 121;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { 0.002, 0.95, -0.363, 0.084, -0.646, 0.816, -0.407, 0.099, -0.02, -0.906, -0.874, 0.191, -0.328, -0.968, 0.79, 0.826, -0.795, 0.277 };
   double X[] = { -0.591, -0.084, 0.707, 0.945, 0.326, 0.428 };
   int incX = -1;
   double x_expected[] = { -1.737656, 0.290416, 0.61669, 0.73853, 0.326, 0.428 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 971) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 971) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 122;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { 0.002, 0.95, -0.363, 0.084, -0.646, 0.816, -0.407, 0.099, -0.02, -0.906, -0.874, 0.191, -0.328, -0.968, 0.79, 0.826, -0.795, 0.277 };
   double X[] = { -0.591, -0.084, 0.707, 0.945, 0.326, 0.428 };
   int incX = -1;
   double x_expected[] = { 0.27516, -0.544536, -0.10627, -0.988374, 0.229991, -0.711267 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 972) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 972) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 122;
   int diag = 132;
   int N = 3;
   int K = 1;
   int lda = 3;
   double A[] = { 0.002, 0.95, -0.363, 0.084, -0.646, 0.816, -0.407, 0.099, -0.02, -0.906, -0.874, 0.191, -0.328, -0.968, 0.79, 0.826, -0.795, 0.277 };
   double X[] = { -0.591, -0.084, 0.707, 0.945, 0.326, 0.428 };
   int incX = -1;
   double x_expected[] = { -0.591, -0.084, 0.794924, 0.411234, 0.148739, 0.025577 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 973) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 973) imag");
     };
   };
  };


}
