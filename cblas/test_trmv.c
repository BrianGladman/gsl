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
   float A[] = { 0.882 };
   float X[] = { 0.872 };
   int incX = 1;
   float x_expected[] = { 0.769104 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 560)");
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
   float A[] = { 0.882 };
   float X[] = { 0.872 };
   int incX = 1;
   float x_expected[] = { 0.872 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 561)");
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
   float A[] = { 0.882 };
   float X[] = { 0.872 };
   int incX = 1;
   float x_expected[] = { 0.769104 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 562)");
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
   float A[] = { 0.882 };
   float X[] = { 0.872 };
   int incX = 1;
   float x_expected[] = { 0.872 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 563)");
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
   float A[] = { 0.882 };
   float X[] = { 0.872 };
   int incX = 1;
   float x_expected[] = { 0.769104 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 564)");
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
   float A[] = { 0.882 };
   float X[] = { 0.872 };
   int incX = 1;
   float x_expected[] = { 0.872 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 565)");
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
   float A[] = { 0.882 };
   float X[] = { 0.872 };
   int incX = 1;
   float x_expected[] = { 0.769104 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 566)");
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
   float A[] = { 0.882 };
   float X[] = { 0.872 };
   int incX = 1;
   float x_expected[] = { 0.872 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 567)");
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
   float A[] = { 0.637 };
   float X[] = { 0.526 };
   int incX = 1;
   float x_expected[] = { 0.335062 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 568)");
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
   float A[] = { 0.637 };
   float X[] = { 0.526 };
   int incX = 1;
   float x_expected[] = { 0.526 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 569)");
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
   float A[] = { 0.637 };
   float X[] = { 0.526 };
   int incX = 1;
   float x_expected[] = { 0.335062 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 570)");
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
   float A[] = { 0.637 };
   float X[] = { 0.526 };
   int incX = 1;
   float x_expected[] = { 0.526 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 571)");
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
   float A[] = { 0.637 };
   float X[] = { 0.526 };
   int incX = 1;
   float x_expected[] = { 0.335062 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 572)");
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
   float A[] = { 0.637 };
   float X[] = { 0.526 };
   int incX = 1;
   float x_expected[] = { 0.526 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 573)");
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
   float A[] = { 0.637 };
   float X[] = { 0.526 };
   int incX = 1;
   float x_expected[] = { 0.335062 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 574)");
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
   float A[] = { 0.637 };
   float X[] = { 0.526 };
   int incX = 1;
   float x_expected[] = { 0.526 };
   cblas_strmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strmv(case 575)");
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
   double A[] = { -0.694 };
   double X[] = { 0.65 };
   int incX = 1;
   double x_expected[] = { -0.4511 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 576)");
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
   double A[] = { -0.694 };
   double X[] = { 0.65 };
   int incX = 1;
   double x_expected[] = { 0.65 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 577)");
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
   double A[] = { -0.694 };
   double X[] = { 0.65 };
   int incX = 1;
   double x_expected[] = { -0.4511 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 578)");
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
   double A[] = { -0.694 };
   double X[] = { 0.65 };
   int incX = 1;
   double x_expected[] = { 0.65 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 579)");
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
   double A[] = { -0.694 };
   double X[] = { 0.65 };
   int incX = 1;
   double x_expected[] = { -0.4511 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 580)");
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
   double A[] = { -0.694 };
   double X[] = { 0.65 };
   int incX = 1;
   double x_expected[] = { 0.65 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 581)");
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
   double A[] = { -0.694 };
   double X[] = { 0.65 };
   int incX = 1;
   double x_expected[] = { -0.4511 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 582)");
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
   double A[] = { -0.694 };
   double X[] = { 0.65 };
   int incX = 1;
   double x_expected[] = { 0.65 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 583)");
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
   double A[] = { 0.222 };
   double X[] = { -0.117 };
   int incX = 1;
   double x_expected[] = { -0.025974 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 584)");
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
   double A[] = { 0.222 };
   double X[] = { -0.117 };
   int incX = 1;
   double x_expected[] = { -0.117 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 585)");
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
   double A[] = { 0.222 };
   double X[] = { -0.117 };
   int incX = 1;
   double x_expected[] = { -0.025974 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 586)");
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
   double A[] = { 0.222 };
   double X[] = { -0.117 };
   int incX = 1;
   double x_expected[] = { -0.117 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 587)");
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
   double A[] = { 0.222 };
   double X[] = { -0.117 };
   int incX = 1;
   double x_expected[] = { -0.025974 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 588)");
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
   double A[] = { 0.222 };
   double X[] = { -0.117 };
   int incX = 1;
   double x_expected[] = { -0.117 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 589)");
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
   double A[] = { 0.222 };
   double X[] = { -0.117 };
   int incX = 1;
   double x_expected[] = { -0.025974 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 590)");
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
   double A[] = { 0.222 };
   double X[] = { -0.117 };
   int incX = 1;
   double x_expected[] = { -0.117 };
   cblas_dtrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrmv(case 591)");
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
   float A[] = { -0.574, -0.031 };
   float X[] = { 0.95, -0.484 };
   int incX = 1;
   float x_expected[] = { -0.560304, 0.248366 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 592) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 592) imag");
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
   float A[] = { -0.574, -0.031 };
   float X[] = { 0.95, -0.484 };
   int incX = 1;
   float x_expected[] = { 0.95, -0.484 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 593) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 593) imag");
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
   float A[] = { -0.574, -0.031 };
   float X[] = { 0.95, -0.484 };
   int incX = 1;
   float x_expected[] = { -0.560304, 0.248366 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 594) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 594) imag");
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
   float A[] = { -0.574, -0.031 };
   float X[] = { 0.95, -0.484 };
   int incX = 1;
   float x_expected[] = { 0.95, -0.484 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 595) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 595) imag");
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
   float A[] = { -0.574, -0.031 };
   float X[] = { 0.95, -0.484 };
   int incX = 1;
   float x_expected[] = { -0.560304, 0.248366 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 596) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 596) imag");
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
   float A[] = { -0.574, -0.031 };
   float X[] = { 0.95, -0.484 };
   int incX = 1;
   float x_expected[] = { 0.95, -0.484 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 597) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 597) imag");
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
   float A[] = { -0.574, -0.031 };
   float X[] = { 0.95, -0.484 };
   int incX = 1;
   float x_expected[] = { -0.560304, 0.248366 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 598) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 598) imag");
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
   float A[] = { -0.574, -0.031 };
   float X[] = { 0.95, -0.484 };
   int incX = 1;
   float x_expected[] = { 0.95, -0.484 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 599) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 599) imag");
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
   float A[] = { -0.912, 0.627 };
   float X[] = { 0.252, 0.133 };
   int incX = 1;
   float x_expected[] = { -0.313215, 0.036708 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 600) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 600) imag");
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
   float A[] = { -0.912, 0.627 };
   float X[] = { 0.252, 0.133 };
   int incX = 1;
   float x_expected[] = { 0.252, 0.133 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 601) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 601) imag");
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
   float A[] = { -0.912, 0.627 };
   float X[] = { 0.252, 0.133 };
   int incX = 1;
   float x_expected[] = { -0.313215, 0.036708 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 602) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 602) imag");
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
   float A[] = { -0.912, 0.627 };
   float X[] = { 0.252, 0.133 };
   int incX = 1;
   float x_expected[] = { 0.252, 0.133 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 603) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 603) imag");
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
   float A[] = { -0.912, 0.627 };
   float X[] = { 0.252, 0.133 };
   int incX = 1;
   float x_expected[] = { -0.313215, 0.036708 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 604) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 604) imag");
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
   float A[] = { -0.912, 0.627 };
   float X[] = { 0.252, 0.133 };
   int incX = 1;
   float x_expected[] = { 0.252, 0.133 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 605) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 605) imag");
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
   float A[] = { -0.912, 0.627 };
   float X[] = { 0.252, 0.133 };
   int incX = 1;
   float x_expected[] = { -0.313215, 0.036708 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 606) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 606) imag");
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
   float A[] = { -0.912, 0.627 };
   float X[] = { 0.252, 0.133 };
   int incX = 1;
   float x_expected[] = { 0.252, 0.133 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 607) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 607) imag");
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
   float A[] = { 0.046, 0.027 };
   float X[] = { 0.453, 0.305 };
   int incX = 1;
   float x_expected[] = { 0.029073, 0.001799 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 608) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 608) imag");
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
   float A[] = { 0.046, 0.027 };
   float X[] = { 0.453, 0.305 };
   int incX = 1;
   float x_expected[] = { 0.453, 0.305 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 609) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 609) imag");
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
   float A[] = { 0.046, 0.027 };
   float X[] = { 0.453, 0.305 };
   int incX = 1;
   float x_expected[] = { 0.029073, 0.001799 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 610) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 610) imag");
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
   float A[] = { 0.046, 0.027 };
   float X[] = { 0.453, 0.305 };
   int incX = 1;
   float x_expected[] = { 0.453, 0.305 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 611) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 611) imag");
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
   float A[] = { 0.046, 0.027 };
   float X[] = { 0.453, 0.305 };
   int incX = 1;
   float x_expected[] = { 0.029073, 0.001799 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 612) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 612) imag");
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
   float A[] = { 0.046, 0.027 };
   float X[] = { 0.453, 0.305 };
   int incX = 1;
   float x_expected[] = { 0.453, 0.305 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 613) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 613) imag");
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
   float A[] = { 0.046, 0.027 };
   float X[] = { 0.453, 0.305 };
   int incX = 1;
   float x_expected[] = { 0.029073, 0.001799 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 614) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 614) imag");
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
   float A[] = { 0.046, 0.027 };
   float X[] = { 0.453, 0.305 };
   int incX = 1;
   float x_expected[] = { 0.453, 0.305 };
   cblas_ctrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrmv(case 615) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrmv(case 615) imag");
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
   double A[] = { -0.797, -0.418 };
   double X[] = { -0.82, 0.616 };
   int incX = 1;
   double x_expected[] = { 0.911028, -0.148192 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 616) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 616) imag");
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
   double A[] = { -0.797, -0.418 };
   double X[] = { -0.82, 0.616 };
   int incX = 1;
   double x_expected[] = { -0.82, 0.616 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 617) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 617) imag");
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
   double A[] = { -0.797, -0.418 };
   double X[] = { -0.82, 0.616 };
   int incX = 1;
   double x_expected[] = { 0.911028, -0.148192 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 618) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 618) imag");
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
   double A[] = { -0.797, -0.418 };
   double X[] = { -0.82, 0.616 };
   int incX = 1;
   double x_expected[] = { -0.82, 0.616 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 619) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 619) imag");
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
   double A[] = { -0.797, -0.418 };
   double X[] = { -0.82, 0.616 };
   int incX = 1;
   double x_expected[] = { 0.911028, -0.148192 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 620) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 620) imag");
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
   double A[] = { -0.797, -0.418 };
   double X[] = { -0.82, 0.616 };
   int incX = 1;
   double x_expected[] = { -0.82, 0.616 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 621) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 621) imag");
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
   double A[] = { -0.797, -0.418 };
   double X[] = { -0.82, 0.616 };
   int incX = 1;
   double x_expected[] = { 0.911028, -0.148192 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 622) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 622) imag");
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
   double A[] = { -0.797, -0.418 };
   double X[] = { -0.82, 0.616 };
   int incX = 1;
   double x_expected[] = { -0.82, 0.616 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 623) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 623) imag");
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
   double A[] = { 0.905, -0.611 };
   double X[] = { 0.439, 0.5 };
   int incX = 1;
   double x_expected[] = { 0.702795, 0.184271 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 624) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 624) imag");
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
   double A[] = { 0.905, -0.611 };
   double X[] = { 0.439, 0.5 };
   int incX = 1;
   double x_expected[] = { 0.439, 0.5 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 625) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 625) imag");
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
   double A[] = { 0.905, -0.611 };
   double X[] = { 0.439, 0.5 };
   int incX = 1;
   double x_expected[] = { 0.702795, 0.184271 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 626) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 626) imag");
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
   double A[] = { 0.905, -0.611 };
   double X[] = { 0.439, 0.5 };
   int incX = 1;
   double x_expected[] = { 0.439, 0.5 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 627) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 627) imag");
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
   double A[] = { 0.905, -0.611 };
   double X[] = { 0.439, 0.5 };
   int incX = 1;
   double x_expected[] = { 0.702795, 0.184271 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 628) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 628) imag");
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
   double A[] = { 0.905, -0.611 };
   double X[] = { 0.439, 0.5 };
   int incX = 1;
   double x_expected[] = { 0.439, 0.5 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 629) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 629) imag");
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
   double A[] = { 0.905, -0.611 };
   double X[] = { 0.439, 0.5 };
   int incX = 1;
   double x_expected[] = { 0.702795, 0.184271 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 630) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 630) imag");
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
   double A[] = { 0.905, -0.611 };
   double X[] = { 0.439, 0.5 };
   int incX = 1;
   double x_expected[] = { 0.439, 0.5 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 631) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 631) imag");
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
   double A[] = { -0.536, 0.949 };
   double X[] = { 0.029, -0.591 };
   int incX = 1;
   double x_expected[] = { -0.576403, 0.289255 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 632) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 632) imag");
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
   double A[] = { -0.536, 0.949 };
   double X[] = { 0.029, -0.591 };
   int incX = 1;
   double x_expected[] = { 0.029, -0.591 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 633) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 633) imag");
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
   double A[] = { -0.536, 0.949 };
   double X[] = { 0.029, -0.591 };
   int incX = 1;
   double x_expected[] = { -0.576403, 0.289255 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 634) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 634) imag");
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
   double A[] = { -0.536, 0.949 };
   double X[] = { 0.029, -0.591 };
   int incX = 1;
   double x_expected[] = { 0.029, -0.591 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 635) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 635) imag");
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
   double A[] = { -0.536, 0.949 };
   double X[] = { 0.029, -0.591 };
   int incX = 1;
   double x_expected[] = { -0.576403, 0.289255 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 636) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 636) imag");
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
   double A[] = { -0.536, 0.949 };
   double X[] = { 0.029, -0.591 };
   int incX = 1;
   double x_expected[] = { 0.029, -0.591 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 637) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 637) imag");
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
   double A[] = { -0.536, 0.949 };
   double X[] = { 0.029, -0.591 };
   int incX = 1;
   double x_expected[] = { -0.576403, 0.289255 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 638) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 638) imag");
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
   double A[] = { -0.536, 0.949 };
   double X[] = { 0.029, -0.591 };
   int incX = 1;
   double x_expected[] = { 0.029, -0.591 };
   cblas_ztrmv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrmv(case 639) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrmv(case 639) imag");
     };
   };
  };


}
