#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_tpmv () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   float A[] = { -0.23 };
   float X[] = { 0.628 };
   int incX = 1;
   float x_expected[] = { -0.14444 };
   cblas_stpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpmv(case 720)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   float A[] = { -0.23 };
   float X[] = { 0.628 };
   int incX = 1;
   float x_expected[] = { 0.628 };
   cblas_stpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpmv(case 721)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   float A[] = { -0.23 };
   float X[] = { 0.628 };
   int incX = 1;
   float x_expected[] = { -0.14444 };
   cblas_stpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpmv(case 722)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   float A[] = { -0.23 };
   float X[] = { 0.628 };
   int incX = 1;
   float x_expected[] = { 0.628 };
   cblas_stpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpmv(case 723)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   float A[] = { -0.23 };
   float X[] = { 0.628 };
   int incX = 1;
   float x_expected[] = { -0.14444 };
   cblas_stpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpmv(case 724)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   float A[] = { -0.23 };
   float X[] = { 0.628 };
   int incX = 1;
   float x_expected[] = { 0.628 };
   cblas_stpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpmv(case 725)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   float A[] = { -0.23 };
   float X[] = { 0.628 };
   int incX = 1;
   float x_expected[] = { -0.14444 };
   cblas_stpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpmv(case 726)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   float A[] = { -0.23 };
   float X[] = { 0.628 };
   int incX = 1;
   float x_expected[] = { 0.628 };
   cblas_stpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpmv(case 727)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   float A[] = { -0.426 };
   float X[] = { -0.535 };
   int incX = 1;
   float x_expected[] = { 0.22791 };
   cblas_stpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpmv(case 728)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   float A[] = { -0.426 };
   float X[] = { -0.535 };
   int incX = 1;
   float x_expected[] = { -0.535 };
   cblas_stpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpmv(case 729)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   float A[] = { -0.426 };
   float X[] = { -0.535 };
   int incX = 1;
   float x_expected[] = { 0.22791 };
   cblas_stpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpmv(case 730)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   float A[] = { -0.426 };
   float X[] = { -0.535 };
   int incX = 1;
   float x_expected[] = { -0.535 };
   cblas_stpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpmv(case 731)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   float A[] = { -0.426 };
   float X[] = { -0.535 };
   int incX = 1;
   float x_expected[] = { 0.22791 };
   cblas_stpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpmv(case 732)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   float A[] = { -0.426 };
   float X[] = { -0.535 };
   int incX = 1;
   float x_expected[] = { -0.535 };
   cblas_stpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpmv(case 733)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   float A[] = { -0.426 };
   float X[] = { -0.535 };
   int incX = 1;
   float x_expected[] = { 0.22791 };
   cblas_stpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpmv(case 734)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   float A[] = { -0.426 };
   float X[] = { -0.535 };
   int incX = 1;
   float x_expected[] = { -0.535 };
   cblas_stpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpmv(case 735)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   double A[] = { -0.291 };
   double X[] = { 0.461 };
   int incX = 1;
   double x_expected[] = { -0.134151 };
   cblas_dtpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpmv(case 736)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   double A[] = { -0.291 };
   double X[] = { 0.461 };
   int incX = 1;
   double x_expected[] = { 0.461 };
   cblas_dtpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpmv(case 737)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   double A[] = { -0.291 };
   double X[] = { 0.461 };
   int incX = 1;
   double x_expected[] = { -0.134151 };
   cblas_dtpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpmv(case 738)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   double A[] = { -0.291 };
   double X[] = { 0.461 };
   int incX = 1;
   double x_expected[] = { 0.461 };
   cblas_dtpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpmv(case 739)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   double A[] = { -0.291 };
   double X[] = { 0.461 };
   int incX = 1;
   double x_expected[] = { -0.134151 };
   cblas_dtpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpmv(case 740)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   double A[] = { -0.291 };
   double X[] = { 0.461 };
   int incX = 1;
   double x_expected[] = { 0.461 };
   cblas_dtpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpmv(case 741)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   double A[] = { -0.291 };
   double X[] = { 0.461 };
   int incX = 1;
   double x_expected[] = { -0.134151 };
   cblas_dtpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpmv(case 742)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   double A[] = { -0.291 };
   double X[] = { 0.461 };
   int incX = 1;
   double x_expected[] = { 0.461 };
   cblas_dtpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpmv(case 743)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   double A[] = { -0.162 };
   double X[] = { -0.449 };
   int incX = 1;
   double x_expected[] = { 0.072738 };
   cblas_dtpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpmv(case 744)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   double A[] = { -0.162 };
   double X[] = { -0.449 };
   int incX = 1;
   double x_expected[] = { -0.449 };
   cblas_dtpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpmv(case 745)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   double A[] = { -0.162 };
   double X[] = { -0.449 };
   int incX = 1;
   double x_expected[] = { 0.072738 };
   cblas_dtpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpmv(case 746)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   double A[] = { -0.162 };
   double X[] = { -0.449 };
   int incX = 1;
   double x_expected[] = { -0.449 };
   cblas_dtpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpmv(case 747)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   double A[] = { -0.162 };
   double X[] = { -0.449 };
   int incX = 1;
   double x_expected[] = { 0.072738 };
   cblas_dtpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpmv(case 748)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   double A[] = { -0.162 };
   double X[] = { -0.449 };
   int incX = 1;
   double x_expected[] = { -0.449 };
   cblas_dtpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpmv(case 749)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   double A[] = { -0.162 };
   double X[] = { -0.449 };
   int incX = 1;
   double x_expected[] = { 0.072738 };
   cblas_dtpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpmv(case 750)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   double A[] = { -0.162 };
   double X[] = { -0.449 };
   int incX = 1;
   double x_expected[] = { -0.449 };
   cblas_dtpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpmv(case 751)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   float A[] = { -0.623, 0.351 };
   float X[] = { 0.393, -0.642 };
   int incX = 1;
   float x_expected[] = { -0.019497, 0.537909 };
   cblas_ctpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpmv(case 752) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpmv(case 752) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   float A[] = { -0.623, 0.351 };
   float X[] = { 0.393, -0.642 };
   int incX = 1;
   float x_expected[] = { 0.393, -0.642 };
   cblas_ctpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpmv(case 753) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpmv(case 753) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   float A[] = { -0.623, 0.351 };
   float X[] = { 0.393, -0.642 };
   int incX = 1;
   float x_expected[] = { -0.019497, 0.537909 };
   cblas_ctpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpmv(case 754) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpmv(case 754) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   float A[] = { -0.623, 0.351 };
   float X[] = { 0.393, -0.642 };
   int incX = 1;
   float x_expected[] = { 0.393, -0.642 };
   cblas_ctpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpmv(case 755) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpmv(case 755) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   float A[] = { -0.623, 0.351 };
   float X[] = { 0.393, -0.642 };
   int incX = 1;
   float x_expected[] = { -0.019497, 0.537909 };
   cblas_ctpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpmv(case 756) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpmv(case 756) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   float A[] = { -0.623, 0.351 };
   float X[] = { 0.393, -0.642 };
   int incX = 1;
   float x_expected[] = { 0.393, -0.642 };
   cblas_ctpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpmv(case 757) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpmv(case 757) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   float A[] = { -0.623, 0.351 };
   float X[] = { 0.393, -0.642 };
   int incX = 1;
   float x_expected[] = { -0.019497, 0.537909 };
   cblas_ctpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpmv(case 758) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpmv(case 758) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   float A[] = { -0.623, 0.351 };
   float X[] = { 0.393, -0.642 };
   int incX = 1;
   float x_expected[] = { 0.393, -0.642 };
   cblas_ctpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpmv(case 759) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpmv(case 759) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   float A[] = { -0.48, -0.88 };
   float X[] = { -0.546, -0.883 };
   int incX = 1;
   float x_expected[] = { -0.51496, 0.90432 };
   cblas_ctpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpmv(case 760) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpmv(case 760) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   float A[] = { -0.48, -0.88 };
   float X[] = { -0.546, -0.883 };
   int incX = 1;
   float x_expected[] = { -0.546, -0.883 };
   cblas_ctpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpmv(case 761) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpmv(case 761) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   float A[] = { -0.48, -0.88 };
   float X[] = { -0.546, -0.883 };
   int incX = 1;
   float x_expected[] = { -0.51496, 0.90432 };
   cblas_ctpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpmv(case 762) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpmv(case 762) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   float A[] = { -0.48, -0.88 };
   float X[] = { -0.546, -0.883 };
   int incX = 1;
   float x_expected[] = { -0.546, -0.883 };
   cblas_ctpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpmv(case 763) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpmv(case 763) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   float A[] = { -0.48, -0.88 };
   float X[] = { -0.546, -0.883 };
   int incX = 1;
   float x_expected[] = { -0.51496, 0.90432 };
   cblas_ctpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpmv(case 764) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpmv(case 764) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   float A[] = { -0.48, -0.88 };
   float X[] = { -0.546, -0.883 };
   int incX = 1;
   float x_expected[] = { -0.546, -0.883 };
   cblas_ctpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpmv(case 765) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpmv(case 765) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   float A[] = { -0.48, -0.88 };
   float X[] = { -0.546, -0.883 };
   int incX = 1;
   float x_expected[] = { -0.51496, 0.90432 };
   cblas_ctpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpmv(case 766) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpmv(case 766) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   float A[] = { -0.48, -0.88 };
   float X[] = { -0.546, -0.883 };
   int incX = 1;
   float x_expected[] = { -0.546, -0.883 };
   cblas_ctpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpmv(case 767) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpmv(case 767) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   float A[] = { 0.753, -0.326 };
   float X[] = { 0.553, -0.975 };
   int incX = 1;
   float x_expected[] = { 0.734259, -0.553897 };
   cblas_ctpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpmv(case 768) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpmv(case 768) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   float A[] = { 0.753, -0.326 };
   float X[] = { 0.553, -0.975 };
   int incX = 1;
   float x_expected[] = { 0.553, -0.975 };
   cblas_ctpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpmv(case 769) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpmv(case 769) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   float A[] = { 0.753, -0.326 };
   float X[] = { 0.553, -0.975 };
   int incX = 1;
   float x_expected[] = { 0.734259, -0.553897 };
   cblas_ctpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpmv(case 770) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpmv(case 770) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   float A[] = { 0.753, -0.326 };
   float X[] = { 0.553, -0.975 };
   int incX = 1;
   float x_expected[] = { 0.553, -0.975 };
   cblas_ctpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpmv(case 771) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpmv(case 771) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   float A[] = { 0.753, -0.326 };
   float X[] = { 0.553, -0.975 };
   int incX = 1;
   float x_expected[] = { 0.734259, -0.553897 };
   cblas_ctpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpmv(case 772) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpmv(case 772) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   float A[] = { 0.753, -0.326 };
   float X[] = { 0.553, -0.975 };
   int incX = 1;
   float x_expected[] = { 0.553, -0.975 };
   cblas_ctpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpmv(case 773) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpmv(case 773) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   float A[] = { 0.753, -0.326 };
   float X[] = { 0.553, -0.975 };
   int incX = 1;
   float x_expected[] = { 0.734259, -0.553897 };
   cblas_ctpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpmv(case 774) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpmv(case 774) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   float A[] = { 0.753, -0.326 };
   float X[] = { 0.553, -0.975 };
   int incX = 1;
   float x_expected[] = { 0.553, -0.975 };
   cblas_ctpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpmv(case 775) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpmv(case 775) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   double A[] = { -0.515, 0.819 };
   double X[] = { -0.202, 0.576 };
   int incX = 1;
   double x_expected[] = { -0.367714, -0.462078 };
   cblas_ztpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpmv(case 776) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpmv(case 776) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   double A[] = { -0.515, 0.819 };
   double X[] = { -0.202, 0.576 };
   int incX = 1;
   double x_expected[] = { -0.202, 0.576 };
   cblas_ztpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpmv(case 777) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpmv(case 777) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   double A[] = { -0.515, 0.819 };
   double X[] = { -0.202, 0.576 };
   int incX = 1;
   double x_expected[] = { -0.367714, -0.462078 };
   cblas_ztpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpmv(case 778) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpmv(case 778) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   double A[] = { -0.515, 0.819 };
   double X[] = { -0.202, 0.576 };
   int incX = 1;
   double x_expected[] = { -0.202, 0.576 };
   cblas_ztpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpmv(case 779) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpmv(case 779) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   double A[] = { -0.515, 0.819 };
   double X[] = { -0.202, 0.576 };
   int incX = 1;
   double x_expected[] = { -0.367714, -0.462078 };
   cblas_ztpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpmv(case 780) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpmv(case 780) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   double A[] = { -0.515, 0.819 };
   double X[] = { -0.202, 0.576 };
   int incX = 1;
   double x_expected[] = { -0.202, 0.576 };
   cblas_ztpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpmv(case 781) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpmv(case 781) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   double A[] = { -0.515, 0.819 };
   double X[] = { -0.202, 0.576 };
   int incX = 1;
   double x_expected[] = { -0.367714, -0.462078 };
   cblas_ztpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpmv(case 782) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpmv(case 782) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   double A[] = { -0.515, 0.819 };
   double X[] = { -0.202, 0.576 };
   int incX = 1;
   double x_expected[] = { -0.202, 0.576 };
   cblas_ztpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpmv(case 783) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpmv(case 783) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   double A[] = { 0.99, -0.144 };
   double X[] = { -0.917, 0.717 };
   int incX = 1;
   double x_expected[] = { -0.804582, 0.841878 };
   cblas_ztpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpmv(case 784) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpmv(case 784) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   double A[] = { 0.99, -0.144 };
   double X[] = { -0.917, 0.717 };
   int incX = 1;
   double x_expected[] = { -0.917, 0.717 };
   cblas_ztpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpmv(case 785) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpmv(case 785) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   double A[] = { 0.99, -0.144 };
   double X[] = { -0.917, 0.717 };
   int incX = 1;
   double x_expected[] = { -0.804582, 0.841878 };
   cblas_ztpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpmv(case 786) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpmv(case 786) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   double A[] = { 0.99, -0.144 };
   double X[] = { -0.917, 0.717 };
   int incX = 1;
   double x_expected[] = { -0.917, 0.717 };
   cblas_ztpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpmv(case 787) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpmv(case 787) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   double A[] = { 0.99, -0.144 };
   double X[] = { -0.917, 0.717 };
   int incX = 1;
   double x_expected[] = { -0.804582, 0.841878 };
   cblas_ztpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpmv(case 788) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpmv(case 788) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   double A[] = { 0.99, -0.144 };
   double X[] = { -0.917, 0.717 };
   int incX = 1;
   double x_expected[] = { -0.917, 0.717 };
   cblas_ztpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpmv(case 789) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpmv(case 789) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   double A[] = { 0.99, -0.144 };
   double X[] = { -0.917, 0.717 };
   int incX = 1;
   double x_expected[] = { -0.804582, 0.841878 };
   cblas_ztpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpmv(case 790) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpmv(case 790) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   double A[] = { 0.99, -0.144 };
   double X[] = { -0.917, 0.717 };
   int incX = 1;
   double x_expected[] = { -0.917, 0.717 };
   cblas_ztpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpmv(case 791) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpmv(case 791) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   double A[] = { -0.857, -0.679 };
   double X[] = { 0.175, -0.67 };
   int incX = 1;
   double x_expected[] = { 0.304955, 0.693015 };
   cblas_ztpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpmv(case 792) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpmv(case 792) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   double A[] = { -0.857, -0.679 };
   double X[] = { 0.175, -0.67 };
   int incX = 1;
   double x_expected[] = { 0.175, -0.67 };
   cblas_ztpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpmv(case 793) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpmv(case 793) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   double A[] = { -0.857, -0.679 };
   double X[] = { 0.175, -0.67 };
   int incX = 1;
   double x_expected[] = { 0.304955, 0.693015 };
   cblas_ztpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpmv(case 794) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpmv(case 794) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   double A[] = { -0.857, -0.679 };
   double X[] = { 0.175, -0.67 };
   int incX = 1;
   double x_expected[] = { 0.175, -0.67 };
   cblas_ztpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpmv(case 795) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpmv(case 795) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   double A[] = { -0.857, -0.679 };
   double X[] = { 0.175, -0.67 };
   int incX = 1;
   double x_expected[] = { 0.304955, 0.693015 };
   cblas_ztpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpmv(case 796) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpmv(case 796) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   double A[] = { -0.857, -0.679 };
   double X[] = { 0.175, -0.67 };
   int incX = 1;
   double x_expected[] = { 0.175, -0.67 };
   cblas_ztpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpmv(case 797) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpmv(case 797) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   double A[] = { -0.857, -0.679 };
   double X[] = { 0.175, -0.67 };
   int incX = 1;
   double x_expected[] = { 0.304955, 0.693015 };
   cblas_ztpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpmv(case 798) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpmv(case 798) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   double A[] = { -0.857, -0.679 };
   double X[] = { 0.175, -0.67 };
   int incX = 1;
   double x_expected[] = { 0.175, -0.67 };
   cblas_ztpmv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpmv(case 799) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpmv(case 799) imag");
     };
   };
  };


}
