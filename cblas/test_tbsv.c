#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_tbsv () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 2;
   int K = 1;
   int lda = 2;
   float A[] = { 0.799, 0.389, 0.471, -0.082 };
   float X[] = { 0.629, 0.801 };
   int incX = 1;
   float x_expected[] = { 9.40311221584, -9.76829268293 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 976)");
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
   float A[] = { 0.799, 0.389, 0.471, -0.082 };
   float X[] = { 0.629, 0.801 };
   int incX = 1;
   float x_expected[] = { 0.317411, 0.801 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 977)");
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
   float A[] = { 0.799, 0.389, 0.471, -0.082 };
   float X[] = { 0.629, 0.801 };
   int incX = 1;
   float x_expected[] = { 0.787234042553, 1.10594541377 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 978)");
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
   float A[] = { 0.799, 0.389, 0.471, -0.082 };
   float X[] = { 0.629, 0.801 };
   int incX = 1;
   float x_expected[] = { 0.629, 0.504741 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 979)");
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
   float A[] = { 0.799, 0.389, 0.471, -0.082 };
   float X[] = { 0.629, 0.801 };
   int incX = 1;
   float x_expected[] = { 13.4443852279, -9.76829268293 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 980)");
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
   float A[] = { 0.799, 0.389, 0.471, -0.082 };
   float X[] = { 0.629, 0.801 };
   int incX = 1;
   float x_expected[] = { 0.251729, 0.801 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 981)");
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
   float A[] = { 0.799, 0.389, 0.471, -0.082 };
   float X[] = { 0.629, 0.801 };
   int incX = 1;
   float x_expected[] = { 0.787234042553, 1.05045850838 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 982)");
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
   float A[] = { 0.799, 0.389, 0.471, -0.082 };
   float X[] = { 0.629, 0.801 };
   int incX = 1;
   float x_expected[] = { 0.629, 0.556319 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 983)");
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
   float A[] = { 0.778, -0.073, -0.976, -0.682 };
   float X[] = { -0.465, 0.29 };
   int incX = 1;
   float x_expected[] = { 0.47643442623, -0.476216588145 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 984)");
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
   float A[] = { 0.778, -0.073, -0.976, -0.682 };
   float X[] = { -0.465, 0.29 };
   int incX = 1;
   float x_expected[] = { -0.465, 0.256055 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 985)");
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
   float A[] = { 0.778, -0.073, -0.976, -0.682 };
   float X[] = { -0.465, 0.29 };
   int incX = 1;
   float x_expected[] = { -5.58131140613, -3.97260273973 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 986)");
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
   float A[] = { 0.778, -0.073, -0.976, -0.682 };
   float X[] = { -0.465, 0.29 };
   int incX = 1;
   float x_expected[] = { -0.18196, 0.29 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 987)");
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
   float A[] = { 0.778, -0.073, -0.976, -0.682 };
   float X[] = { -0.465, 0.29 };
   int incX = 1;
   float x_expected[] = { 6.3698630137, -9.54103563251 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 988)");
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
   float A[] = { 0.778, -0.073, -0.976, -0.682 };
   float X[] = { -0.465, 0.29 };
   int incX = 1;
   float x_expected[] = { -0.465, -0.16384 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 989)");
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
   float A[] = { 0.778, -0.073, -0.976, -0.682 };
   float X[] = { -0.465, 0.29 };
   int incX = 1;
   float x_expected[] = { -0.625566290193, -0.297131147541 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 990)");
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
   float A[] = { 0.778, -0.073, -0.976, -0.682 };
   float X[] = { -0.465, 0.29 };
   int incX = 1;
   float x_expected[] = { -0.44383, 0.29 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 991)");
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
   double A[] = { 0.145, 0.902, -0.474, 0.784 };
   double X[] = { 0.932, -0.724 };
   int incX = 1;
   double x_expected[] = { -3.72356410919, -0.923469387755 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 992)");
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
   double A[] = { 0.145, 0.902, -0.474, 0.784 };
   double X[] = { 0.932, -0.724 };
   int incX = 1;
   double x_expected[] = { 1.585048, -0.724 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 993)");
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
   double A[] = { 0.145, 0.902, -0.474, 0.784 };
   double X[] = { 0.932, -0.724 };
   int incX = 1;
   double x_expected[] = { 6.4275862069, 2.57502867192 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 994)");
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
   double A[] = { 0.145, 0.902, -0.474, 0.784 };
   double X[] = { 0.932, -0.724 };
   int incX = 1;
   double x_expected[] = { 0.932, -0.282232 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 995)");
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
   double A[] = { 0.145, 0.902, -0.474, 0.784 };
   double X[] = { 0.932, -0.724 };
   int incX = 1;
   double x_expected[] = { 0.54797728404, -0.923469387755 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 996)");
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
   double A[] = { 0.145, 0.902, -0.474, 0.784 };
   double X[] = { 0.932, -0.724 };
   int incX = 1;
   double x_expected[] = { 0.588824, -0.724 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 997)");
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
   double A[] = { 0.145, 0.902, -0.474, 0.784 };
   double X[] = { 0.932, -0.724 };
   int incX = 1;
   double x_expected[] = { 6.4275862069, 13.7588243853 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 998)");
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
   double A[] = { 0.145, 0.902, -0.474, 0.784 };
   double X[] = { 0.932, -0.724 };
   int incX = 1;
   double x_expected[] = { 0.932, -1.564664 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 999)");
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
   double A[] = { 0.334, -0.317, 0.348, 0.07 };
   double X[] = { -0.627, -0.226 };
   int incX = 1;
   double x_expected[] = { -1.80172413793, -11.3878078818 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 1000)");
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
   double A[] = { 0.334, -0.317, 0.348, 0.07 };
   double X[] = { -0.627, -0.226 };
   int incX = 1;
   double x_expected[] = { -0.627, -0.424759 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 1001)");
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
   double A[] = { 0.334, -0.317, 0.348, 0.07 };
   double X[] = { -0.627, -0.226 };
   int incX = 1;
   double x_expected[] = { -2.62006271369, 0.712933753943 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 1002)");
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
   double A[] = { 0.334, -0.317, 0.348, 0.07 };
   double X[] = { -0.627, -0.226 };
   int incX = 1;
   double x_expected[] = { -0.548352, -0.226 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 1003)");
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
   double A[] = { 0.334, -0.317, 0.348, 0.07 };
   double X[] = { -0.627, -0.226 };
   int incX = 1;
   double x_expected[] = { 1.97791798107, -13.0616493916 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 1004)");
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
   double A[] = { 0.334, -0.317, 0.348, 0.07 };
   double X[] = { -0.627, -0.226 };
   int incX = 1;
   double x_expected[] = { -0.627, -0.007804 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 1005)");
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
   double A[] = { 0.334, -0.317, 0.348, 0.07 };
   double X[] = { -0.627, -0.226 };
   int incX = 1;
   double x_expected[] = { -2.49361621584, -0.649425287356 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 1006)");
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
   double A[] = { 0.334, -0.317, 0.348, 0.07 };
   double X[] = { -0.627, -0.226 };
   int incX = 1;
   double x_expected[] = { -0.698642, -0.226 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 1007)");
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
   float A[] = { 0.677, 0.423, 0.249, -0.143, -0.135, -0.182, -0.689, -0.076 };
   float X[] = { -0.855, 0.389, 0.599, -0.473 };
   int incX = 1;
   float x_expected[] = { 1.726022562, -2.95209538305, -0.784111035032, 0.772993379771 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1008) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1008) imag");
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
   float A[] = { 0.677, 0.423, 0.249, -0.143, -0.135, -0.182, -0.689, -0.076 };
   float X[] = { -0.855, 0.389, 0.599, -0.473 };
   int incX = 1;
   float x_expected[] = { -0.936512, 0.592434, 0.599, -0.473 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1009) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1009) imag");
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
   float A[] = { 0.677, 0.423, 0.249, -0.143, -0.135, -0.182, -0.689, -0.076 };
   float X[] = { -0.855, 0.389, 0.599, -0.473 };
   int incX = 1;
   float x_expected[] = { -0.650110316387, 0.980792708762, 1.80078284525, -0.808839819394 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1010) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1010) imag");
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
   float A[] = { 0.677, 0.423, 0.249, -0.143, -0.135, -0.182, -0.689, -0.076 };
   float X[] = { -0.855, 0.389, 0.599, -0.473 };
   int incX = 1;
   float x_expected[] = { -0.855, 0.389, 0.412777, -0.576095 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1011) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1011) imag");
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
   float A[] = { 0.677, 0.423, 0.249, -0.143, -0.135, -0.182, -0.689, -0.076 };
   float X[] = { -0.855, 0.389, 0.599, -0.473 };
   int incX = 1;
   float x_expected[] = { -3.93481831202, -0.851538637451, -0.784111035032, 0.772993379771 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1012) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1012) imag");
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
   float A[] = { 0.677, 0.423, 0.249, -0.143, -0.135, -0.182, -0.689, -0.076 };
   float X[] = { -0.855, 0.389, 0.599, -0.473 };
   int incX = 1;
   float x_expected[] = { -0.688049, 0.434163, 0.599, -0.473 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1013) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1013) imag");
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
   float A[] = { 0.677, 0.423, 0.249, -0.143, -0.135, -0.182, -0.689, -0.076 };
   float X[] = { -0.855, 0.389, 0.599, -0.473 };
   int incX = 1;
   float x_expected[] = { -0.650110316387, 0.980792708762, 1.23992833409, 4.32974965126 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1014) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1014) imag");
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
   float A[] = { 0.677, 0.423, 0.249, -0.143, -0.135, -0.182, -0.689, -0.076 };
   float X[] = { -0.855, 0.389, 0.599, -0.473 };
   int incX = 1;
   float x_expected[] = { -0.855, 0.389, 0.756268, -0.692126 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1015) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1015) imag");
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
   float A[] = { 0.493, 0.112, -0.902, 0.128, -0.889, -0.277, -0.256, -0.777 };
   float X[] = { 0.488, -0.633, 0.029, 0.84 };
   int incX = 1;
   float x_expected[] = { -0.298126982296, 0.804928204832, -1.81008488252, -0.772537268646 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1016) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1016) imag");
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
   float A[] = { 0.493, 0.112, -0.902, 0.128, -0.889, -0.277, -0.256, -0.777 };
   float X[] = { 0.488, -0.633, 0.029, 0.84 };
   int incX = 1;
   float x_expected[] = { 0.488, -0.633, 0.388152, 0.20657 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1017) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1017) imag");
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
   float A[] = { 0.493, 0.112, -0.902, 0.128, -0.889, -0.277, -0.256, -0.777 };
   float X[] = { 0.488, -0.633, 0.029, 0.84 };
   int incX = 1;
   float x_expected[] = { 0.976663064999, -3.1049885701, 0.0980279232953, -0.917353021971 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1018) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1018) imag");
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
   float A[] = { 0.493, 0.112, -0.902, 0.128, -0.889, -0.277, -0.256, -0.777 };
   float X[] = { 0.488, -0.633, 0.029, 0.84 };
   int incX = 1;
   float x_expected[] = { 0.281101, 0.121793, 0.029, 0.84 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1019) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1019) imag");
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
   float A[] = { 0.493, 0.112, -0.902, 0.128, -0.889, -0.277, -0.256, -0.777 };
   float X[] = { 0.488, -0.633, 0.029, 0.84 };
   int incX = 1;
   float x_expected[] = { -0.627960886181, 0.61266186981, -1.1382452424, -1.27458860722 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1020) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1020) imag");
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
   float A[] = { 0.493, 0.112, -0.902, 0.128, -0.889, -0.277, -0.256, -0.777 };
   float X[] = { 0.488, -0.633, 0.029, 0.84 };
   int incX = 1;
   float x_expected[] = { 0.488, -0.633, 0.638173, 0.412439 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1021) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1021) imag");
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
   float A[] = { 0.493, 0.112, -0.902, 0.128, -0.889, -0.277, -0.256, -0.777 };
   float X[] = { 0.488, -0.633, 0.029, 0.84 };
   int incX = 1;
   float x_expected[] = { -0.38511691793, -2.67792189938, -0.298092382216, -0.852000461334 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1022) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1022) imag");
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
   float A[] = { 0.493, 0.112, -0.902, 0.128, -0.889, -0.277, -0.256, -0.777 };
   float X[] = { 0.488, -0.633, 0.029, 0.84 };
   int incX = 1;
   float x_expected[] = { 0.621678, 0.120968, 0.029, 0.84 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1023) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1023) imag");
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
   float A[] = { 0.874, 0.452, 0.322, -0.066, -0.477, -0.153, 0.416, -0.619 };
   float X[] = { 0.174, 0.543, 0.777, 0.614 };
   int incX = 1;
   float x_expected[] = { 0.000322788896062, -1.13826124381, 1.61603372146, -0.0476587912069 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1024) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1024) imag");
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
   float A[] = { 0.874, 0.452, 0.322, -0.066, -0.477, -0.153, 0.416, -0.619 };
   float X[] = { 0.174, 0.543, 0.777, 0.614 };
   int incX = 1;
   float x_expected[] = { 0.174, 0.543, 0.75681, 0.42767 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1025) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1025) imag");
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
   float A[] = { 0.874, 0.452, 0.322, -0.066, -0.477, -0.153, 0.416, -0.619 };
   float X[] = { 0.174, 0.543, 0.777, 0.614 };
   int incX = 1;
   float x_expected[] = { 1.13982704022, 1.47938123676, 2.69083672714, 1.35529433543 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1026) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1026) imag");
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
   float A[] = { 0.874, 0.452, 0.322, -0.066, -0.477, -0.153, 0.416, -0.619 };
   float X[] = { 0.174, 0.543, 0.777, 0.614 };
   int incX = 1;
   float x_expected[] = { 0.638571, 0.716997, 0.777, 0.614 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1027) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1027) imag");
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
   float A[] = { 0.874, 0.452, 0.322, -0.066, -0.477, -0.153, 0.416, -0.619 };
   float X[] = { 0.174, 0.543, 0.777, 0.614 };
   int incX = 1;
   float x_expected[] = { 0.850296186598, 1.51205109219, 2.39868026598, -0.672187572639 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1028) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1028) imag");
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
   float A[] = { 0.874, 0.452, 0.322, -0.066, -0.477, -0.153, 0.416, -0.619 };
   float X[] = { 0.174, 0.543, 0.777, 0.614 };
   int incX = 1;
   float x_expected[] = { 0.174, 0.543, 0.943077, 0.846389 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1029) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1029) imag");
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
   float A[] = { 0.874, 0.452, 0.322, -0.066, -0.477, -0.153, 0.416, -0.619 };
   float X[] = { 0.174, 0.543, 0.777, 0.614 };
   int incX = 1;
   float x_expected[] = { -0.154330935169, 1.22926535115, -1.10261100351, -1.64087942041 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1030) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1030) imag");
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
   float A[] = { 0.874, 0.452, 0.322, -0.066, -0.477, -0.153, 0.416, -0.619 };
   float X[] = { 0.174, 0.543, 0.777, 0.614 };
   int incX = 1;
   float x_expected[] = { -0.03567, 0.29401, 0.777, 0.614 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1031) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1031) imag");
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
   double A[] = { -0.446, 0.222, -0.138, 0.165, -0.767, -0.063, 0.725, -0.047 };
   double X[] = { 0.267, 0.102, -0.88, -0.278 };
   int incX = 1;
   double x_expected[] = { -0.0607220214106, -0.299895695129, -1.18395935086, -0.460201502745 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1032) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1032) imag");
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
   double A[] = { -0.446, 0.222, -0.138, 0.165, -0.767, -0.063, 0.725, -0.047 };
   double X[] = { 0.267, 0.102, -0.88, -0.278 };
   int incX = 1;
   double x_expected[] = { 0.09969, 0.208836, -0.88, -0.278 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1033) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1033) imag");
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
   double A[] = { -0.446, 0.222, -0.138, 0.165, -0.767, -0.063, 0.725, -0.047 };
   double X[] = { 0.267, 0.102, -0.88, -0.278 };
   int incX = 1;
   double x_expected[] = { -0.388549556809, -0.422103142627, 1.20098588277, 5.97387248645 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1034) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1034) imag");
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
   double A[] = { -0.446, 0.222, -0.138, 0.165, -0.767, -0.063, 0.725, -0.047 };
   double X[] = { 0.267, 0.102, -0.88, -0.278 };
   int incX = 1;
   double x_expected[] = { 0.267, 0.102, -0.681637, -0.182945 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1035) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1035) imag");
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
   double A[] = { -0.446, 0.222, -0.138, 0.165, -0.767, -0.063, 0.725, -0.047 };
   double X[] = { 0.267, 0.102, -0.88, -0.278 };
   int incX = 1;
   double x_expected[] = { 0.664641789402, 3.15383976059, -1.18395935086, -0.460201502745 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1036) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1036) imag");
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
   double A[] = { -0.446, 0.222, -0.138, 0.165, -0.767, -0.063, 0.725, -0.047 };
   double X[] = { 0.267, 0.102, -0.88, -0.278 };
   int incX = 1;
   double x_expected[] = { -0.390446, -0.166666, -0.88, -0.278 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1037) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1037) imag");
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
   double A[] = { -0.446, 0.222, -0.138, 0.165, -0.767, -0.063, 0.725, -0.047 };
   double X[] = { 0.267, 0.102, -0.88, -0.278 };
   int incX = 1;
   double x_expected[] = { -0.388549556809, -0.422103142627, 1.32822261866, 0.245712557801 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1038) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1038) imag");
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
   double A[] = { -0.446, 0.222, -0.138, 0.165, -0.767, -0.063, 0.725, -0.047 };
   double X[] = { 0.267, 0.102, -0.88, -0.278 };
   int incX = 1;
   double x_expected[] = { 0.267, 0.102, -0.826324, -0.307979 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1039) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1039) imag");
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
   double A[] = { 0.403, -0.838, -0.096, -0.337, -0.007, -0.688, -0.657, 0.29 };
   double X[] = { 0.442, 0.313, -0.001, -0.073 };
   int incX = 1;
   double x_expected[] = { -0.461430566147, 0.637747072728, 0.237426654378, 0.359409606666 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1040) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1040) imag");
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
   double A[] = { 0.403, -0.838, -0.096, -0.337, -0.007, -0.688, -0.657, 0.29 };
   double X[] = { 0.442, 0.313, -0.001, -0.073 };
   int incX = 1;
   double x_expected[] = { 0.442, 0.313, -0.064049, 0.106002 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1041) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1041) imag");
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
   double A[] = { 0.403, -0.838, -0.096, -0.337, -0.007, -0.688, -0.657, 0.29 };
   double X[] = { 0.442, 0.313, -0.001, -0.073 };
   int incX = 1;
   double x_expected[] = { -0.248596169936, 0.604072420401, 0.201140204422, 0.0543307407257 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1042) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1042) imag");
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
   double A[] = { 0.403, -0.838, -0.096, -0.337, -0.007, -0.688, -0.657, 0.29 };
   double X[] = { 0.442, 0.313, -0.001, -0.073 };
   int incX = 1;
   double x_expected[] = { 0.492217, 0.311801, -0.001, -0.073 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1043) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1043) imag");
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
   double A[] = { 0.403, -0.838, -0.096, -0.337, -0.007, -0.688, -0.657, 0.29 };
   double X[] = { 0.442, 0.313, -0.001, -0.073 };
   int incX = 1;
   double x_expected[] = { -1.20465040518, 0.968408193183, 0.357494061751, 1.52008203854 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1044) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1044) imag");
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
   double A[] = { 0.403, -0.838, -0.096, -0.337, -0.007, -0.688, -0.657, 0.29 };
   double X[] = { 0.442, 0.313, -0.001, -0.073 };
   int incX = 1;
   double x_expected[] = { 0.442, 0.313, -0.21325, 0.233287 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1045) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1045) imag");
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
   double A[] = { 0.403, -0.838, -0.096, -0.337, -0.007, -0.688, -0.657, 0.29 };
   double X[] = { 0.442, 0.313, -0.001, -0.073 };
   int incX = 1;
   double x_expected[] = { -0.127158157121, 0.600903522868, 0.106108455343, -0.000373896529944 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1046) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1046) imag");
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
   double A[] = { 0.403, -0.838, -0.096, -0.337, -0.007, -0.688, -0.657, 0.29 };
   double X[] = { 0.442, 0.313, -0.001, -0.073 };
   int incX = 1;
   double x_expected[] = { 0.466505, 0.305655, -0.001, -0.073 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1047) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1047) imag");
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
   double A[] = { 0.991, 0.571, -0.654, 0.932, -0.12, -0.179, 0.416, -0.724 };
   double X[] = { 0.492, 0.187, -0.965, -0.338 };
   int incX = 1;
   double x_expected[] = { -0.550526474451, -2.37953532439, -1.96711301421, -2.36326190877 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1048) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1048) imag");
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
   double A[] = { 0.991, 0.571, -0.654, 0.932, -0.12, -0.179, 0.416, -0.724 };
   double X[] = { 0.492, 0.187, -0.965, -0.338 };
   int incX = 1;
   double x_expected[] = { 0.492, 0.187, -0.817516, 0.242842 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1049) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1049) imag");
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
   double A[] = { 0.991, 0.571, -0.654, 0.932, -0.12, -0.179, 0.416, -0.724 };
   double X[] = { 0.492, 0.187, -0.965, -0.338 };
   int incX = 1;
   double x_expected[] = { 0.370926759835, 0.207230499326, 0.729844022402, -0.523263958529 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1050) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1050) imag");
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
   double A[] = { 0.991, 0.571, -0.654, 0.932, -0.12, -0.179, 0.416, -0.724 };
   double X[] = { 0.492, 0.187, -0.965, -0.338 };
   int incX = 1;
   double x_expected[] = { 0.315698, 0.319175, -0.965, -0.338 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1051) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1051) imag");
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
   double A[] = { 0.991, 0.571, -0.654, 0.932, -0.12, -0.179, 0.416, -0.724 };
   double X[] = { 0.492, 0.187, -0.965, -0.338 };
   int incX = 1;
   double x_expected[] = { -0.382655784748, 0.259381026583, -0.822988539455, 0.8592904116 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1052) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1052) imag");
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
   double A[] = { 0.991, 0.571, -0.654, 0.932, -0.12, -0.179, 0.416, -0.724 };
   double X[] = { 0.492, 0.187, -0.965, -0.338 };
   int incX = 1;
   double x_expected[] = { 0.492, 0.187, -0.872487, -0.403628 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1053) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1053) imag");
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
   double A[] = { 0.991, 0.571, -0.654, 0.932, -0.12, -0.179, 0.416, -0.724 };
   double X[] = { 0.492, 0.187, -0.965, -0.338 };
   int incX = 1;
   double x_expected[] = { -4.15729238485, 1.94413182536, 1.19071510088, 4.59281669215 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1054) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1054) imag");
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
   double A[] = { 0.991, 0.571, -0.654, 0.932, -0.12, -0.179, 0.416, -0.724 };
   double X[] = { 0.492, 0.187, -0.965, -0.338 };
   int incX = 1;
   double x_expected[] = { 0.175906, -0.933432, -0.965, -0.338 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1055) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1055) imag");
     };
   };
  };


}
