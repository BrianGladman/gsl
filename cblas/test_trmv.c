#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_trmv () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { 0.987 };
   float X[] = { -0.138 };
   int incX = -1;
   float x_expected[] = { -0.136206 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 814)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { 0.987 };
   float X[] = { -0.138 };
   int incX = -1;
   float x_expected[] = { -0.138 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 815)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { 0.987 };
   float X[] = { -0.138 };
   int incX = -1;
   float x_expected[] = { -0.136206 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 816)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { 0.987 };
   float X[] = { -0.138 };
   int incX = -1;
   float x_expected[] = { -0.138 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 817)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { 0.987 };
   float X[] = { -0.138 };
   int incX = -1;
   float x_expected[] = { -0.136206 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 818)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { 0.987 };
   float X[] = { -0.138 };
   int incX = -1;
   float x_expected[] = { -0.138 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 819)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { 0.987 };
   float X[] = { -0.138 };
   int incX = -1;
   float x_expected[] = { -0.136206 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 820)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { 0.987 };
   float X[] = { -0.138 };
   int incX = -1;
   float x_expected[] = { -0.138 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 821)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { -0.329 };
   float X[] = { 0.463 };
   int incX = -1;
   float x_expected[] = { -0.152327 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 822)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { -0.329 };
   float X[] = { 0.463 };
   int incX = -1;
   float x_expected[] = { 0.463 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 823)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { -0.329 };
   float X[] = { 0.463 };
   int incX = -1;
   float x_expected[] = { -0.152327 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 824)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { -0.329 };
   float X[] = { 0.463 };
   int incX = -1;
   float x_expected[] = { 0.463 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 825)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { -0.329 };
   float X[] = { 0.463 };
   int incX = -1;
   float x_expected[] = { -0.152327 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 826)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { -0.329 };
   float X[] = { 0.463 };
   int incX = -1;
   float x_expected[] = { 0.463 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 827)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { -0.329 };
   float X[] = { 0.463 };
   int incX = -1;
   float x_expected[] = { -0.152327 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 828)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { -0.329 };
   float X[] = { 0.463 };
   int incX = -1;
   float x_expected[] = { 0.463 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 829)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { -0.429 };
   double X[] = { -0.899 };
   int incX = -1;
   double x_expected[] = { 0.385671 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 830)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { -0.429 };
   double X[] = { -0.899 };
   int incX = -1;
   double x_expected[] = { -0.899 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 831)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { -0.429 };
   double X[] = { -0.899 };
   int incX = -1;
   double x_expected[] = { 0.385671 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 832)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { -0.429 };
   double X[] = { -0.899 };
   int incX = -1;
   double x_expected[] = { -0.899 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 833)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { -0.429 };
   double X[] = { -0.899 };
   int incX = -1;
   double x_expected[] = { 0.385671 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 834)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { -0.429 };
   double X[] = { -0.899 };
   int incX = -1;
   double x_expected[] = { -0.899 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 835)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { -0.429 };
   double X[] = { -0.899 };
   int incX = -1;
   double x_expected[] = { 0.385671 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 836)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { -0.429 };
   double X[] = { -0.899 };
   int incX = -1;
   double x_expected[] = { -0.899 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 837)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { 0.842 };
   double X[] = { 0.192 };
   int incX = -1;
   double x_expected[] = { 0.161664 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 838)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { 0.842 };
   double X[] = { 0.192 };
   int incX = -1;
   double x_expected[] = { 0.192 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 839)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { 0.842 };
   double X[] = { 0.192 };
   int incX = -1;
   double x_expected[] = { 0.161664 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 840)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { 0.842 };
   double X[] = { 0.192 };
   int incX = -1;
   double x_expected[] = { 0.192 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 841)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { 0.842 };
   double X[] = { 0.192 };
   int incX = -1;
   double x_expected[] = { 0.161664 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 842)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { 0.842 };
   double X[] = { 0.192 };
   int incX = -1;
   double x_expected[] = { 0.192 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 843)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { 0.842 };
   double X[] = { 0.192 };
   int incX = -1;
   double x_expected[] = { 0.161664 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 844)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { 0.842 };
   double X[] = { 0.192 };
   int incX = -1;
   double x_expected[] = { 0.192 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 845)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { -0.162, -0.108 };
   float X[] = { 0.542, 0.461 };
   int incX = -1;
   float x_expected[] = { -0.038016, -0.133218 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 846) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 846) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { -0.162, -0.108 };
   float X[] = { 0.542, 0.461 };
   int incX = -1;
   float x_expected[] = { 0.542, 0.461 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 847) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 847) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { -0.162, -0.108 };
   float X[] = { 0.542, 0.461 };
   int incX = -1;
   float x_expected[] = { -0.038016, -0.133218 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 848) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 848) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { -0.162, -0.108 };
   float X[] = { 0.542, 0.461 };
   int incX = -1;
   float x_expected[] = { 0.542, 0.461 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 849) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 849) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { -0.162, -0.108 };
   float X[] = { 0.542, 0.461 };
   int incX = -1;
   float x_expected[] = { -0.038016, -0.133218 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 850) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 850) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { -0.162, -0.108 };
   float X[] = { 0.542, 0.461 };
   int incX = -1;
   float x_expected[] = { 0.542, 0.461 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 851) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 851) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { -0.162, -0.108 };
   float X[] = { 0.542, 0.461 };
   int incX = -1;
   float x_expected[] = { -0.038016, -0.133218 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 852) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 852) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { -0.162, -0.108 };
   float X[] = { 0.542, 0.461 };
   int incX = -1;
   float x_expected[] = { 0.542, 0.461 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 853) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 853) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { 0.547, 0.583 };
   float X[] = { -0.302, 0.434 };
   int incX = -1;
   float x_expected[] = { -0.418216, 0.061332 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 854) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 854) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { 0.547, 0.583 };
   float X[] = { -0.302, 0.434 };
   int incX = -1;
   float x_expected[] = { -0.302, 0.434 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 855) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 855) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { 0.547, 0.583 };
   float X[] = { -0.302, 0.434 };
   int incX = -1;
   float x_expected[] = { -0.418216, 0.061332 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 856) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 856) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { 0.547, 0.583 };
   float X[] = { -0.302, 0.434 };
   int incX = -1;
   float x_expected[] = { -0.302, 0.434 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 857) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 857) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { 0.547, 0.583 };
   float X[] = { -0.302, 0.434 };
   int incX = -1;
   float x_expected[] = { -0.418216, 0.061332 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 858) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 858) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { 0.547, 0.583 };
   float X[] = { -0.302, 0.434 };
   int incX = -1;
   float x_expected[] = { -0.302, 0.434 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 859) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 859) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { 0.547, 0.583 };
   float X[] = { -0.302, 0.434 };
   int incX = -1;
   float x_expected[] = { -0.418216, 0.061332 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 860) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 860) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { 0.547, 0.583 };
   float X[] = { -0.302, 0.434 };
   int incX = -1;
   float x_expected[] = { -0.302, 0.434 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 861) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 861) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { 0.216, 0.192 };
   float X[] = { -0.564, -0.297 };
   int incX = -1;
   float x_expected[] = { -0.178848, 0.044136 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 862) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 862) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { 0.216, 0.192 };
   float X[] = { -0.564, -0.297 };
   int incX = -1;
   float x_expected[] = { -0.564, -0.297 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 863) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 863) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { 0.216, 0.192 };
   float X[] = { -0.564, -0.297 };
   int incX = -1;
   float x_expected[] = { -0.178848, 0.044136 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 864) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 864) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { 0.216, 0.192 };
   float X[] = { -0.564, -0.297 };
   int incX = -1;
   float x_expected[] = { -0.564, -0.297 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 865) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 865) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { 0.216, 0.192 };
   float X[] = { -0.564, -0.297 };
   int incX = -1;
   float x_expected[] = { -0.178848, 0.044136 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 866) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 866) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { 0.216, 0.192 };
   float X[] = { -0.564, -0.297 };
   int incX = -1;
   float x_expected[] = { -0.564, -0.297 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 867) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 867) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { 0.216, 0.192 };
   float X[] = { -0.564, -0.297 };
   int incX = -1;
   float x_expected[] = { -0.178848, 0.044136 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 868) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 868) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   int lda = 1;
   float A[] = { 0.216, 0.192 };
   float X[] = { -0.564, -0.297 };
   int incX = -1;
   float x_expected[] = { -0.564, -0.297 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 869) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 869) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { 0.693, -0.22 };
   double X[] = { -0.101, 0.889 };
   int incX = -1;
   double x_expected[] = { 0.125587, 0.638297 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 870) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 870) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { 0.693, -0.22 };
   double X[] = { -0.101, 0.889 };
   int incX = -1;
   double x_expected[] = { -0.101, 0.889 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 871) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 871) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { 0.693, -0.22 };
   double X[] = { -0.101, 0.889 };
   int incX = -1;
   double x_expected[] = { 0.125587, 0.638297 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 872) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 872) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { 0.693, -0.22 };
   double X[] = { -0.101, 0.889 };
   int incX = -1;
   double x_expected[] = { -0.101, 0.889 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 873) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 873) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { 0.693, -0.22 };
   double X[] = { -0.101, 0.889 };
   int incX = -1;
   double x_expected[] = { 0.125587, 0.638297 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 874) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 874) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { 0.693, -0.22 };
   double X[] = { -0.101, 0.889 };
   int incX = -1;
   double x_expected[] = { -0.101, 0.889 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 875) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 875) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { 0.693, -0.22 };
   double X[] = { -0.101, 0.889 };
   int incX = -1;
   double x_expected[] = { 0.125587, 0.638297 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 876) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 876) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { 0.693, -0.22 };
   double X[] = { -0.101, 0.889 };
   int incX = -1;
   double x_expected[] = { -0.101, 0.889 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 877) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 877) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { -0.216, -0.623 };
   double X[] = { 0.048, 0.293 };
   int incX = -1;
   double x_expected[] = { 0.172171, -0.093192 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 878) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 878) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { -0.216, -0.623 };
   double X[] = { 0.048, 0.293 };
   int incX = -1;
   double x_expected[] = { 0.048, 0.293 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 879) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 879) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { -0.216, -0.623 };
   double X[] = { 0.048, 0.293 };
   int incX = -1;
   double x_expected[] = { 0.172171, -0.093192 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 880) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 880) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { -0.216, -0.623 };
   double X[] = { 0.048, 0.293 };
   int incX = -1;
   double x_expected[] = { 0.048, 0.293 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 881) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 881) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { -0.216, -0.623 };
   double X[] = { 0.048, 0.293 };
   int incX = -1;
   double x_expected[] = { 0.172171, -0.093192 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 882) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 882) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { -0.216, -0.623 };
   double X[] = { 0.048, 0.293 };
   int incX = -1;
   double x_expected[] = { 0.048, 0.293 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 883) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 883) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { -0.216, -0.623 };
   double X[] = { 0.048, 0.293 };
   int incX = -1;
   double x_expected[] = { 0.172171, -0.093192 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 884) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 884) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { -0.216, -0.623 };
   double X[] = { 0.048, 0.293 };
   int incX = -1;
   double x_expected[] = { 0.048, 0.293 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 885) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 885) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { -0.345, -0.851 };
   double X[] = { -0.708, 0.298 };
   int incX = -1;
   double x_expected[] = { -0.009338, -0.705318 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 886) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 886) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { -0.345, -0.851 };
   double X[] = { -0.708, 0.298 };
   int incX = -1;
   double x_expected[] = { -0.708, 0.298 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 887) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 887) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { -0.345, -0.851 };
   double X[] = { -0.708, 0.298 };
   int incX = -1;
   double x_expected[] = { -0.009338, -0.705318 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 888) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 888) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { -0.345, -0.851 };
   double X[] = { -0.708, 0.298 };
   int incX = -1;
   double x_expected[] = { -0.708, 0.298 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 889) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 889) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { -0.345, -0.851 };
   double X[] = { -0.708, 0.298 };
   int incX = -1;
   double x_expected[] = { -0.009338, -0.705318 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 890) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 890) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { -0.345, -0.851 };
   double X[] = { -0.708, 0.298 };
   int incX = -1;
   double x_expected[] = { -0.708, 0.298 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 891) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 891) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   int lda = 1;
   double A[] = { -0.345, -0.851 };
   double X[] = { -0.708, 0.298 };
   int incX = -1;
   double x_expected[] = { -0.009338, -0.705318 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 892) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 892) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   int lda = 1;
   double A[] = { -0.345, -0.851 };
   double X[] = { -0.708, 0.298 };
   int incX = -1;
   double x_expected[] = { -0.708, 0.298 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 893) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 893) imag");
     };
   };
  };


}
