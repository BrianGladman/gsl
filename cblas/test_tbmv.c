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
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.39, -0.128, 0.28, -0.984 };
   float X[] = { 0.939, -0.065 };
   int incX = 1;
   float x_expected[] = { 0.27124, 0.06396 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 640)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.39, -0.128, 0.28, -0.984 };
   float X[] = { 0.939, -0.065 };
   int incX = 1;
   float x_expected[] = { 0.94732, -0.065 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 641)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.39, -0.128, 0.28, -0.984 };
   float X[] = { 0.939, -0.065 };
   int incX = 1;
   float x_expected[] = { 0.36621, 0.27124 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 642)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.39, -0.128, 0.28, -0.984 };
   float X[] = { 0.939, -0.065 };
   int incX = 1;
   float x_expected[] = { 0.939, 0.19792 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 643)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.39, -0.128, 0.28, -0.984 };
   float X[] = { 0.939, -0.065 };
   int incX = 1;
   float x_expected[] = { -0.138392, 0.06396 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 644)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.39, -0.128, 0.28, -0.984 };
   float X[] = { 0.939, -0.065 };
   int incX = 1;
   float x_expected[] = { 0.9208, -0.065 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 645)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.39, -0.128, 0.28, -0.984 };
   float X[] = { 0.939, -0.065 };
   int incX = 1;
   float x_expected[] = { 0.36621, -0.138392 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 646)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.39, -0.128, 0.28, -0.984 };
   float X[] = { 0.939, -0.065 };
   int incX = 1;
   float x_expected[] = { 0.939, -0.185192 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 647)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.56, 0.887, -0.887, 0.016 };
   float X[] = { -0.141, 0.815 };
   int incX = 1;
   float x_expected[] = { 0.125067, -0.112027 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 648)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.56, 0.887, -0.887, 0.016 };
   float X[] = { -0.141, 0.815 };
   int incX = 1;
   float x_expected[] = { -0.141, 0.689933 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 649)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.56, 0.887, -0.887, 0.016 };
   float X[] = { -0.141, 0.815 };
   int incX = 1;
   float x_expected[] = { -0.801865, 0.722905 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 650)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.56, 0.887, -0.887, 0.016 };
   float X[] = { -0.141, 0.815 };
   int incX = 1;
   float x_expected[] = { -0.863905, 0.815 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 651)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.56, 0.887, -0.887, 0.016 };
   float X[] = { -0.141, 0.815 };
   int incX = 1;
   float x_expected[] = { -0.125067, 0.138107 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 652)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.56, 0.887, -0.887, 0.016 };
   float X[] = { -0.141, 0.815 };
   int incX = 1;
   float x_expected[] = { -0.141, 0.940067 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 653)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.56, 0.887, -0.887, 0.016 };
   float X[] = { -0.141, 0.815 };
   int incX = 1;
   float x_expected[] = { 0.643945, -0.722905 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 654)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.56, 0.887, -0.887, 0.016 };
   float X[] = { -0.141, 0.815 };
   int incX = 1;
   float x_expected[] = { 0.581905, 0.815 };
   cblas_stbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbmv(case 655)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { -0.674, -0.645, 0.187, 0.833 };
   double X[] = { 0.508, -0.683 };
   int incX = 1;
   double x_expected[] = { 0.535531, -0.568939 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 656)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { -0.674, -0.645, 0.187, 0.833 };
   double X[] = { 0.508, -0.683 };
   int incX = 1;
   double x_expected[] = { 0.948535, -0.683 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 657)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { -0.674, -0.645, 0.187, 0.833 };
   double X[] = { 0.508, -0.683 };
   int incX = 1;
   double x_expected[] = { -0.342392, 0.535531 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 658)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { -0.674, -0.645, 0.187, 0.833 };
   double X[] = { 0.508, -0.683 };
   int incX = 1;
   double x_expected[] = { 0.508, -0.588004 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 659)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { -0.674, -0.645, 0.187, 0.833 };
   double X[] = { 0.508, -0.683 };
   int incX = 1;
   double x_expected[] = { -0.455381, -0.568939 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 660)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { -0.674, -0.645, 0.187, 0.833 };
   double X[] = { 0.508, -0.683 };
   int incX = 1;
   double x_expected[] = { 0.380279, -0.683 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 661)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { -0.674, -0.645, 0.187, 0.833 };
   double X[] = { 0.508, -0.683 };
   int incX = 1;
   double x_expected[] = { -0.342392, -0.455381 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 662)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { -0.674, -0.645, 0.187, 0.833 };
   double X[] = { 0.508, -0.683 };
   int incX = 1;
   double x_expected[] = { 0.508, -1.01066 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 663)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { 0.404, -0.054, 0.021, -0.274 };
   double X[] = { -0.006, 0.372 };
   int incX = 1;
   double x_expected[] = { -1.260000000000e-04, -0.101604 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 664)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { 0.404, -0.054, 0.021, -0.274 };
   double X[] = { -0.006, 0.372 };
   int incX = 1;
   double x_expected[] = { -0.006, 0.372324 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 665)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { 0.404, -0.054, 0.021, -0.274 };
   double X[] = { -0.006, 0.372 };
   int incX = 1;
   double x_expected[] = { 0.005388, -0.020088 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 666)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { 0.404, -0.054, 0.021, -0.274 };
   double X[] = { -0.006, 0.372 };
   int incX = 1;
   double x_expected[] = { 0.001812, 0.372 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 667)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { 0.404, -0.054, 0.021, -0.274 };
   double X[] = { -0.006, 0.372 };
   int incX = 1;
   double x_expected[] = { 3.240000000000e-04, -0.102054 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 668)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { 0.404, -0.054, 0.021, -0.274 };
   double X[] = { -0.006, 0.372 };
   int incX = 1;
   double x_expected[] = { -0.006, 0.371874 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 669)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { 0.404, -0.054, 0.021, -0.274 };
   double X[] = { -0.006, 0.372 };
   int incX = 1;
   double x_expected[] = { -0.022512, 0.007812 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 670)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { 0.404, -0.054, 0.021, -0.274 };
   double X[] = { -0.006, 0.372 };
   int incX = 1;
   double x_expected[] = { -0.026088, 0.372 };
   cblas_dtbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbmv(case 671)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.616, -0.54, -0.145, 0.352, -0.695, -0.805, 0.165, -0.39 };
   float X[] = { -0.223, -0.483, -0.118, -0.682 };
   int incX = 1;
   float x_expected[] = { 0.023344, 0.572554, -0.28545, -0.06651 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 672) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 672) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.616, -0.54, -0.145, 0.352, -0.695, -0.805, 0.165, -0.39 };
   float X[] = { -0.223, -0.483, -0.118, -0.682 };
   int incX = 1;
   float x_expected[] = { 0.034174, -0.425646, -0.118, -0.682 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 673) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 673) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.616, -0.54, -0.145, 0.352, -0.695, -0.805, 0.165, -0.39 };
   float X[] = { -0.223, -0.483, -0.118, -0.682 };
   int incX = 1;
   float x_expected[] = { -0.398188, -0.177108, 0.023344, 0.572554 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 674) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 674) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.616, -0.54, -0.145, 0.352, -0.695, -0.805, 0.165, -0.39 };
   float X[] = { -0.223, -0.483, -0.118, -0.682 };
   int incX = 1;
   float x_expected[] = { -0.223, -0.483, -0.35183, -0.1668 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 675) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 675) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.616, -0.54, -0.145, 0.352, -0.695, -0.805, 0.165, -0.39 };
   float X[] = { -0.223, -0.483, -0.118, -0.682 };
   int incX = 1;
   float x_expected[] = { -0.264649, 0.560519, -0.28545, -0.06651 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 676) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 676) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.616, -0.54, -0.145, 0.352, -0.695, -0.805, 0.165, -0.39 };
   float X[] = { -0.223, -0.483, -0.118, -0.682 };
   int incX = 1;
   float x_expected[] = { -0.69, 0.08598, -0.118, -0.682 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 677) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 677) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.616, -0.54, -0.145, 0.352, -0.695, -0.805, 0.165, -0.39 };
   float X[] = { -0.223, -0.483, -0.118, -0.682 };
   int incX = 1;
   float x_expected[] = { -0.398188, -0.177108, -0.264649, 0.560519 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 678) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 678) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.616, -0.54, -0.145, 0.352, -0.695, -0.805, 0.165, -0.39 };
   float X[] = { -0.223, -0.483, -0.118, -0.682 };
   int incX = 1;
   float x_expected[] = { -0.223, -0.483, 0.084351, -0.690461 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 679) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 679) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.112, 0.162, -0.42, 0.164, -0.209, -0.52, 0.385, 0.558 };
   float X[] = { -0.933, 0.78, -0.04, 0.453 };
   int incX = 1;
   float x_expected[] = { 0.600597, 0.32214, -0.004234, -0.328527 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 680) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 680) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.112, 0.162, -0.42, 0.164, -0.209, -0.52, 0.385, 0.558 };
   float X[] = { -0.933, 0.78, -0.04, 0.453 };
   int incX = 1;
   float x_expected[] = { -0.933, 0.78, 0.22394, -0.027612 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 681) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 681) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.112, 0.162, -0.42, 0.164, -0.209, -0.52, 0.385, 0.558 };
   float X[] = { -0.933, 0.78, -0.04, 0.453 };
   int incX = 1;
   float x_expected[] = { 0.013064, -0.137663, -0.057492, -0.19682 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 682) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 682) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.112, 0.162, -0.42, 0.164, -0.209, -0.52, 0.385, 0.558 };
   float X[] = { -0.933, 0.78, -0.04, 0.453 };
   int incX = 1;
   float x_expected[] = { -0.68908, 0.706123, -0.04, 0.453 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 683) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 683) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.112, 0.162, -0.42, 0.164, -0.209, -0.52, 0.385, 0.558 };
   float X[] = { -0.933, 0.78, -0.04, 0.453 };
   int incX = 1;
   float x_expected[] = { 0.26394, -0.480612, 0.332423, 0.474225 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 684) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 684) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.112, 0.162, -0.42, 0.164, -0.209, -0.52, 0.385, 0.558 };
   float X[] = { -0.933, 0.78, -0.04, 0.453 };
   int incX = 1;
   float x_expected[] = { -0.933, 0.78, 0.560597, 0.77514 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 685) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 685) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.112, 0.162, -0.42, 0.164, -0.209, -0.52, 0.385, 0.558 };
   float X[] = { -0.933, 0.78, -0.04, 0.453 };
   int incX = 1;
   float x_expected[] = { -0.288348, -0.260606, 0.24392, -0.073877 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 686) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 686) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.112, 0.162, -0.42, 0.164, -0.209, -0.52, 0.385, 0.558 };
   float X[] = { -0.933, 0.78, -0.04, 0.453 };
   int incX = 1;
   float x_expected[] = { -0.990492, 0.58318, -0.04, 0.453 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 687) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 687) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 121;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.808, -0.192, -0.223, 0.835, 0.275, -0.857, -0.76, 0.394 };
   float X[] = { 0.934, 0.336, -0.288, 0.381 };
   int incX = 1;
   float x_expected[] = { -0.031102, 0.892838, 0.441272, -1.030906 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 688) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 688) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 121;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.808, -0.192, -0.223, 0.835, 0.275, -0.857, -0.76, 0.394 };
   float X[] = { 0.934, 0.336, -0.288, 0.381 };
   int incX = 1;
   float x_expected[] = { 0.934, 0.336, -0.215722, -0.473818 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 689) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 689) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 122;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.808, -0.192, -0.223, 0.835, 0.275, -0.857, -0.76, 0.394 };
   float X[] = { 0.934, 0.336, -0.288, 0.381 };
   int incX = 1;
   float x_expected[] = { 0.284443, 0.308775, 0.382359, 0.155517 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 690) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 690) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 122;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.808, -0.192, -0.223, 0.835, 0.275, -0.857, -0.76, 0.394 };
   float X[] = { 0.934, 0.336, -0.288, 0.381 };
   int incX = 1;
   float x_expected[] = { 0.528283, 0.193959, -0.288, 0.381 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 691) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 691) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 121;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.808, -0.192, -0.223, 0.835, 0.275, -0.857, -0.76, 0.394 };
   float X[] = { 0.934, 0.336, -0.288, 0.381 };
   int incX = 1;
   float x_expected[] = { 0.072278, -0.854818, 0.337892, 0.71675 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 692) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 692) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 121;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.808, -0.192, -0.223, 0.835, 0.275, -0.857, -0.76, 0.394 };
   float X[] = { 0.934, 0.336, -0.288, 0.381 };
   int incX = 1;
   float x_expected[] = { 0.934, 0.336, -0.319102, 1.273838 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 693) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 693) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 122;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.808, -0.192, -0.223, 0.835, 0.275, -0.857, -0.76, 0.394 };
   float X[] = { 0.934, 0.336, -0.288, 0.381 };
   int incX = 1;
   float x_expected[] = { 1.072519, 0.606333, -0.405717, -0.142041 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 694) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 694) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 122;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.808, -0.192, -0.223, 0.835, 0.275, -0.857, -0.76, 0.394 };
   float X[] = { 0.934, 0.336, -0.288, 0.381 };
   int incX = 1;
   float x_expected[] = { 1.316359, 0.491517, -0.288, 0.381 };
   cblas_ctbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbmv(case 695) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbmv(case 695) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { 0.966, 0.754, -0.403, -0.337, -0.651, 0.941, -0.117, 0.511 };
   double X[] = { 0.118, -0.954, 0.027, -0.421 };
   int incX = 1;
   double x_expected[] = { 0.668138, 0.892656, 0.211972, 0.063054 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 696) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 696) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { 0.966, 0.754, -0.403, -0.337, -0.651, 0.941, -0.117, 0.511 };
   double X[] = { 0.118, -0.954, 0.027, -0.421 };
   int incX = 1;
   double x_expected[] = { -0.034758, -0.793436, 0.027, -0.421 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 697) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 697) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { 0.966, 0.754, -0.403, -0.337, -0.651, 0.941, -0.117, 0.511 };
   double X[] = { 0.118, -0.954, 0.027, -0.421 };
   int incX = 1;
   double x_expected[] = { 0.833304, -0.832592, 0.668138, 0.892656 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 698) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 698) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { 0.966, 0.754, -0.403, -0.337, -0.651, 0.941, -0.117, 0.511 };
   double X[] = { 0.118, -0.954, 0.027, -0.421 };
   int incX = 1;
   double x_expected[] = { 0.118, -0.954, 0.847896, 0.311092 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 699) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 699) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { 0.966, 0.754, -0.403, -0.337, -0.651, 0.941, -0.117, 0.511 };
   double X[] = { 0.118, -0.954, 0.027, -0.421 };
   int incX = 1;
   double x_expected[] = { 0.009532, 0.644174, 0.211972, 0.063054 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 700) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 700) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { 0.966, 0.754, -0.403, -0.337, -0.651, 0.941, -0.117, 0.511 };
   double X[] = { 0.118, -0.954, 0.027, -0.421 };
   int incX = 1;
   double x_expected[] = { 0.496584, -0.654522, 0.027, -0.421 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 701) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 701) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { 0.966, 0.754, -0.403, -0.337, -0.651, 0.941, -0.117, 0.511 };
   double X[] = { 0.118, -0.954, 0.027, -0.421 };
   int incX = 1;
   double x_expected[] = { 0.833304, -0.832592, 0.009532, 0.644174 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 702) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 702) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { 0.966, 0.754, -0.403, -0.337, -0.651, 0.941, -0.117, 0.511 };
   double X[] = { 0.118, -0.954, 0.027, -0.421 };
   int incX = 1;
   double x_expected[] = { 0.118, -0.954, -0.342052, -0.076304 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 703) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 703) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { 0.752, -0.632, 0.963, -0.188, 0.401, -0.339, 0.829, -0.826 };
   double X[] = { 0.25, 0.222, -0.578, 0.183 };
   int incX = 1;
   double x_expected[] = { 0.175508, 0.004272, -0.045518, 0.795921 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 704) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 704) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { 0.752, -0.632, 0.963, -0.188, 0.401, -0.339, 0.829, -0.826 };
   double X[] = { 0.25, 0.222, -0.578, 0.183 };
   int incX = 1;
   double x_expected[] = { 0.25, 0.222, -0.295514, 0.349786 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 705) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 705) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { 0.752, -0.632, 0.963, -0.188, 0.401, -0.339, 0.829, -0.826 };
   double X[] = { 0.25, 0.222, -0.578, 0.183 };
   int incX = 1;
   double x_expected[] = { 0.158563, 0.278269, -0.52221, 0.284893 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 706) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 706) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { 0.752, -0.632, 0.963, -0.188, 0.401, -0.339, 0.829, -0.826 };
   double X[] = { 0.25, 0.222, -0.578, 0.183 };
   int incX = 1;
   double x_expected[] = { 0.080259, 0.491325, -0.578, 0.183 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 707) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 707) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { 0.752, -0.632, 0.963, -0.188, 0.401, -0.339, 0.829, -0.826 };
   double X[] = { 0.25, 0.222, -0.578, 0.183 };
   int incX = 1;
   double x_expected[] = { 0.282486, 0.166786, -0.152496, 0.633407 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 708) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 708) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { 0.752, -0.632, 0.963, -0.188, 0.401, -0.339, 0.829, -0.826 };
   double X[] = { 0.25, 0.222, -0.578, 0.183 };
   int incX = 1;
   double x_expected[] = { 0.25, 0.222, -0.402492, 0.187272 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 709) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 709) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { 0.752, -0.632, 0.963, -0.188, 0.401, -0.339, 0.829, -0.826 };
   double X[] = { 0.25, 0.222, -0.578, 0.183 };
   int incX = 1;
   double x_expected[] = { -0.193906, 0.293837, -0.169741, 0.269325 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 710) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 710) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { 0.752, -0.632, 0.963, -0.188, 0.401, -0.339, 0.829, -0.826 };
   double X[] = { 0.25, 0.222, -0.578, 0.183 };
   int incX = 1;
   double x_expected[] = { -0.27221, 0.506893, -0.578, 0.183 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 711) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 711) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 121;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { -0.177, 0.908, 0.118, 0.826, 0.964, 0.553, -0.787, 0.991 };
   double X[] = { -0.935, -0.111, 0.672, -0.103 };
   int incX = 1;
   double x_expected[] = { -0.962723, 0.410051, -0.832953, 0.174321 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 712) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 712) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 121;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { -0.177, 0.908, 0.118, 0.826, 0.964, 0.553, -0.787, 0.991 };
   double X[] = { -0.935, -0.111, 0.672, -0.103 };
   int incX = 1;
   double x_expected[] = { -0.935, -0.111, 0.469984, 0.656212 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 713) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 713) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 122;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { -0.177, 0.908, 0.118, 0.826, 0.964, 0.553, -0.787, 0.991 };
   double X[] = { -0.935, -0.111, 0.672, -0.103 };
   int incX = 1;
   double x_expected[] = { 0.655556, 0.397719, -0.005782, -0.567226 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 714) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 714) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 122;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { -0.177, 0.908, 0.118, 0.826, 0.964, 0.553, -0.787, 0.991 };
   double X[] = { -0.935, -0.111, 0.672, -0.103 };
   int incX = 1;
   double x_expected[] = { -0.344151, -0.581908, 0.672, -0.103 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 715) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 715) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 121;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { -0.177, 0.908, 0.118, 0.826, 0.964, 0.553, -0.787, 0.991 };
   double X[] = { -0.935, -0.111, 0.672, -0.103 };
   int incX = 1;
   double x_expected[] = { -0.202016, 0.759212, -1.59366, -0.17484 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 716) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 716) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 121;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { -0.177, 0.908, 0.118, 0.826, 0.964, 0.553, -0.787, 0.991 };
   double X[] = { -0.935, -0.111, 0.672, -0.103 };
   int incX = 1;
   double x_expected[] = { -0.935, -0.111, -0.290723, 0.307051 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 717) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 717) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 122;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { -0.177, 0.908, 0.118, 0.826, 0.964, 0.553, -0.787, 0.991 };
   double X[] = { -0.935, -0.111, 0.672, -0.103 };
   int incX = 1;
   double x_expected[] = { 0.058925, 0.301401, 0.590849, -0.470908 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 718) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 718) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 122;
   int diag = 132;
   int N = 2;
   int K = 1;
   int lda = 2;
   double A[] = { -0.177, 0.908, 0.118, 0.826, 0.964, 0.553, -0.787, 0.991 };
   double X[] = { -0.935, -0.111, 0.672, -0.103 };
   int incX = 1;
   double x_expected[] = { -0.940782, -0.678226, 0.672, -0.103 };
   cblas_ztbmv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbmv(case 719) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbmv(case 719) imag");
     };
   };
  };


}
