#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_tpsv () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   float A[] = { -0.976 };
   float X[] = { -0.052 };
   int incX = 1;
   float x_expected[] = { 0.0532786885246 };
   cblas_stpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpsv(case 1056)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   float A[] = { -0.976 };
   float X[] = { -0.052 };
   int incX = 1;
   float x_expected[] = { -0.052 };
   cblas_stpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpsv(case 1057)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   float A[] = { -0.976 };
   float X[] = { -0.052 };
   int incX = 1;
   float x_expected[] = { 0.0532786885246 };
   cblas_stpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpsv(case 1058)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   float A[] = { -0.976 };
   float X[] = { -0.052 };
   int incX = 1;
   float x_expected[] = { -0.052 };
   cblas_stpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpsv(case 1059)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   float A[] = { -0.976 };
   float X[] = { -0.052 };
   int incX = 1;
   float x_expected[] = { 0.0532786885246 };
   cblas_stpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpsv(case 1060)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   float A[] = { -0.976 };
   float X[] = { -0.052 };
   int incX = 1;
   float x_expected[] = { -0.052 };
   cblas_stpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpsv(case 1061)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   float A[] = { -0.976 };
   float X[] = { -0.052 };
   int incX = 1;
   float x_expected[] = { 0.0532786885246 };
   cblas_stpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpsv(case 1062)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   float A[] = { -0.976 };
   float X[] = { -0.052 };
   int incX = 1;
   float x_expected[] = { -0.052 };
   cblas_stpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpsv(case 1063)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   float A[] = { -0.976 };
   float X[] = { -0.052 };
   int incX = 1;
   float x_expected[] = { 0.0532786885246 };
   cblas_stpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpsv(case 1064)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   float A[] = { -0.976 };
   float X[] = { -0.052 };
   int incX = 1;
   float x_expected[] = { -0.052 };
   cblas_stpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpsv(case 1065)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   float A[] = { -0.976 };
   float X[] = { -0.052 };
   int incX = 1;
   float x_expected[] = { 0.0532786885246 };
   cblas_stpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpsv(case 1066)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   float A[] = { -0.976 };
   float X[] = { -0.052 };
   int incX = 1;
   float x_expected[] = { -0.052 };
   cblas_stpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpsv(case 1067)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   float A[] = { -0.976 };
   float X[] = { -0.052 };
   int incX = 1;
   float x_expected[] = { 0.0532786885246 };
   cblas_stpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpsv(case 1068)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   float A[] = { -0.976 };
   float X[] = { -0.052 };
   int incX = 1;
   float x_expected[] = { -0.052 };
   cblas_stpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpsv(case 1069)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   float A[] = { -0.976 };
   float X[] = { -0.052 };
   int incX = 1;
   float x_expected[] = { 0.0532786885246 };
   cblas_stpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpsv(case 1070)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   float A[] = { -0.976 };
   float X[] = { -0.052 };
   int incX = 1;
   float x_expected[] = { -0.052 };
   cblas_stpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stpsv(case 1071)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   double A[] = { -0.101 };
   double X[] = { 0.8 };
   int incX = 1;
   double x_expected[] = { -7.92079207921 };
   cblas_dtpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpsv(case 1072)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   double A[] = { -0.101 };
   double X[] = { 0.8 };
   int incX = 1;
   double x_expected[] = { 0.8 };
   cblas_dtpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpsv(case 1073)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   double A[] = { -0.101 };
   double X[] = { 0.8 };
   int incX = 1;
   double x_expected[] = { -7.92079207921 };
   cblas_dtpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpsv(case 1074)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   double A[] = { -0.101 };
   double X[] = { 0.8 };
   int incX = 1;
   double x_expected[] = { 0.8 };
   cblas_dtpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpsv(case 1075)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   double A[] = { -0.101 };
   double X[] = { 0.8 };
   int incX = 1;
   double x_expected[] = { -7.92079207921 };
   cblas_dtpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpsv(case 1076)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   double A[] = { -0.101 };
   double X[] = { 0.8 };
   int incX = 1;
   double x_expected[] = { 0.8 };
   cblas_dtpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpsv(case 1077)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   double A[] = { -0.101 };
   double X[] = { 0.8 };
   int incX = 1;
   double x_expected[] = { -7.92079207921 };
   cblas_dtpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpsv(case 1078)");
     }
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   double A[] = { -0.101 };
   double X[] = { 0.8 };
   int incX = 1;
   double x_expected[] = { 0.8 };
   cblas_dtpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpsv(case 1079)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   double A[] = { -0.101 };
   double X[] = { 0.8 };
   int incX = 1;
   double x_expected[] = { -7.92079207921 };
   cblas_dtpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpsv(case 1080)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   double A[] = { -0.101 };
   double X[] = { 0.8 };
   int incX = 1;
   double x_expected[] = { 0.8 };
   cblas_dtpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpsv(case 1081)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   double A[] = { -0.101 };
   double X[] = { 0.8 };
   int incX = 1;
   double x_expected[] = { -7.92079207921 };
   cblas_dtpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpsv(case 1082)");
     }
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   double A[] = { -0.101 };
   double X[] = { 0.8 };
   int incX = 1;
   double x_expected[] = { 0.8 };
   cblas_dtpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpsv(case 1083)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   double A[] = { -0.101 };
   double X[] = { 0.8 };
   int incX = 1;
   double x_expected[] = { -7.92079207921 };
   cblas_dtpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpsv(case 1084)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   double A[] = { -0.101 };
   double X[] = { 0.8 };
   int incX = 1;
   double x_expected[] = { 0.8 };
   cblas_dtpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpsv(case 1085)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   double A[] = { -0.101 };
   double X[] = { 0.8 };
   int incX = 1;
   double x_expected[] = { -7.92079207921 };
   cblas_dtpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpsv(case 1086)");
     }
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   double A[] = { -0.101 };
   double X[] = { 0.8 };
   int incX = 1;
   double x_expected[] = { 0.8 };
   cblas_dtpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtpsv(case 1087)");
     }
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   float A[] = { 0.026, -0.335 };
   float X[] = { -0.669, 0.746 };
   int incX = 1;
   float x_expected[] = { -2.36759638976, -1.81326117572 };
   cblas_ctpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpsv(case 1088) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpsv(case 1088) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   float A[] = { 0.026, -0.335 };
   float X[] = { -0.669, 0.746 };
   int incX = 1;
   float x_expected[] = { -0.669, 0.746 };
   cblas_ctpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpsv(case 1089) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpsv(case 1089) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   float A[] = { 0.026, -0.335 };
   float X[] = { -0.669, 0.746 };
   int incX = 1;
   float x_expected[] = { -2.36759638976, -1.81326117572 };
   cblas_ctpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpsv(case 1090) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpsv(case 1090) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   float A[] = { 0.026, -0.335 };
   float X[] = { -0.669, 0.746 };
   int incX = 1;
   float x_expected[] = { -0.669, 0.746 };
   cblas_ctpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpsv(case 1091) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpsv(case 1091) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   float A[] = { 0.026, -0.335 };
   float X[] = { -0.669, 0.746 };
   int incX = 1;
   float x_expected[] = { -2.36759638976, -1.81326117572 };
   cblas_ctpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpsv(case 1092) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpsv(case 1092) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   float A[] = { 0.026, -0.335 };
   float X[] = { -0.669, 0.746 };
   int incX = 1;
   float x_expected[] = { -0.669, 0.746 };
   cblas_ctpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpsv(case 1093) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpsv(case 1093) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   float A[] = { 0.026, -0.335 };
   float X[] = { -0.669, 0.746 };
   int incX = 1;
   float x_expected[] = { -2.36759638976, -1.81326117572 };
   cblas_ctpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpsv(case 1094) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpsv(case 1094) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   float A[] = { 0.026, -0.335 };
   float X[] = { -0.669, 0.746 };
   int incX = 1;
   float x_expected[] = { -0.669, 0.746 };
   cblas_ctpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpsv(case 1095) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpsv(case 1095) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   float A[] = { 0.026, -0.335 };
   float X[] = { -0.669, 0.746 };
   int incX = 1;
   float x_expected[] = { -2.36759638976, -1.81326117572 };
   cblas_ctpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpsv(case 1096) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpsv(case 1096) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   float A[] = { 0.026, -0.335 };
   float X[] = { -0.669, 0.746 };
   int incX = 1;
   float x_expected[] = { -0.669, 0.746 };
   cblas_ctpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpsv(case 1097) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpsv(case 1097) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   float A[] = { 0.026, -0.335 };
   float X[] = { -0.669, 0.746 };
   int incX = 1;
   float x_expected[] = { -2.36759638976, -1.81326117572 };
   cblas_ctpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpsv(case 1098) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpsv(case 1098) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   float A[] = { 0.026, -0.335 };
   float X[] = { -0.669, 0.746 };
   int incX = 1;
   float x_expected[] = { -0.669, 0.746 };
   cblas_ctpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpsv(case 1099) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpsv(case 1099) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   float A[] = { 0.026, -0.335 };
   float X[] = { -0.669, 0.746 };
   int incX = 1;
   float x_expected[] = { -2.36759638976, -1.81326117572 };
   cblas_ctpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpsv(case 1100) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpsv(case 1100) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   float A[] = { 0.026, -0.335 };
   float X[] = { -0.669, 0.746 };
   int incX = 1;
   float x_expected[] = { -0.669, 0.746 };
   cblas_ctpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpsv(case 1101) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpsv(case 1101) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   float A[] = { 0.026, -0.335 };
   float X[] = { -0.669, 0.746 };
   int incX = 1;
   float x_expected[] = { -2.36759638976, -1.81326117572 };
   cblas_ctpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpsv(case 1102) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpsv(case 1102) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   float A[] = { 0.026, -0.335 };
   float X[] = { -0.669, 0.746 };
   int incX = 1;
   float x_expected[] = { -0.669, 0.746 };
   cblas_ctpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpsv(case 1103) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpsv(case 1103) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   float A[] = { 0.026, -0.335 };
   float X[] = { -0.669, 0.746 };
   int incX = 1;
   float x_expected[] = { 2.05946802951, 2.15685423513 };
   cblas_ctpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpsv(case 1104) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpsv(case 1104) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   float A[] = { 0.026, -0.335 };
   float X[] = { -0.669, 0.746 };
   int incX = 1;
   float x_expected[] = { -0.669, 0.746 };
   cblas_ctpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpsv(case 1105) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpsv(case 1105) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   float A[] = { 0.026, -0.335 };
   float X[] = { -0.669, 0.746 };
   int incX = 1;
   float x_expected[] = { 2.05946802951, 2.15685423513 };
   cblas_ctpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpsv(case 1106) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpsv(case 1106) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   float A[] = { 0.026, -0.335 };
   float X[] = { -0.669, 0.746 };
   int incX = 1;
   float x_expected[] = { -0.669, 0.746 };
   cblas_ctpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpsv(case 1107) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpsv(case 1107) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   float A[] = { 0.026, -0.335 };
   float X[] = { -0.669, 0.746 };
   int incX = 1;
   float x_expected[] = { 2.05946802951, 2.15685423513 };
   cblas_ctpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpsv(case 1108) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpsv(case 1108) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   float A[] = { 0.026, -0.335 };
   float X[] = { -0.669, 0.746 };
   int incX = 1;
   float x_expected[] = { -0.669, 0.746 };
   cblas_ctpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpsv(case 1109) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpsv(case 1109) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   float A[] = { 0.026, -0.335 };
   float X[] = { -0.669, 0.746 };
   int incX = 1;
   float x_expected[] = { 2.05946802951, 2.15685423513 };
   cblas_ctpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpsv(case 1110) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpsv(case 1110) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   float A[] = { 0.026, -0.335 };
   float X[] = { -0.669, 0.746 };
   int incX = 1;
   float x_expected[] = { -0.669, 0.746 };
   cblas_ctpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctpsv(case 1111) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctpsv(case 1111) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   double A[] = { 0.262, -0.449 };
   double X[] = { -0.954, -0.093 };
   int incX = 1;
   double x_expected[] = { -0.770378730411, -1.67519103036 };
   cblas_ztpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpsv(case 1112) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpsv(case 1112) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   double A[] = { 0.262, -0.449 };
   double X[] = { -0.954, -0.093 };
   int incX = 1;
   double x_expected[] = { -0.954, -0.093 };
   cblas_ztpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpsv(case 1113) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpsv(case 1113) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   double A[] = { 0.262, -0.449 };
   double X[] = { -0.954, -0.093 };
   int incX = 1;
   double x_expected[] = { -0.770378730411, -1.67519103036 };
   cblas_ztpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpsv(case 1114) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpsv(case 1114) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   double A[] = { 0.262, -0.449 };
   double X[] = { -0.954, -0.093 };
   int incX = 1;
   double x_expected[] = { -0.954, -0.093 };
   cblas_ztpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpsv(case 1115) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpsv(case 1115) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   double A[] = { 0.262, -0.449 };
   double X[] = { -0.954, -0.093 };
   int incX = 1;
   double x_expected[] = { -0.770378730411, -1.67519103036 };
   cblas_ztpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpsv(case 1116) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpsv(case 1116) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   double A[] = { 0.262, -0.449 };
   double X[] = { -0.954, -0.093 };
   int incX = 1;
   double x_expected[] = { -0.954, -0.093 };
   cblas_ztpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpsv(case 1117) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpsv(case 1117) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   double A[] = { 0.262, -0.449 };
   double X[] = { -0.954, -0.093 };
   int incX = 1;
   double x_expected[] = { -0.770378730411, -1.67519103036 };
   cblas_ztpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpsv(case 1118) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpsv(case 1118) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 111;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   double A[] = { 0.262, -0.449 };
   double X[] = { -0.954, -0.093 };
   int incX = 1;
   double x_expected[] = { -0.954, -0.093 };
   cblas_ztpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpsv(case 1119) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpsv(case 1119) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   double A[] = { 0.262, -0.449 };
   double X[] = { -0.954, -0.093 };
   int incX = 1;
   double x_expected[] = { -0.770378730411, -1.67519103036 };
   cblas_ztpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpsv(case 1120) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpsv(case 1120) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   double A[] = { 0.262, -0.449 };
   double X[] = { -0.954, -0.093 };
   int incX = 1;
   double x_expected[] = { -0.954, -0.093 };
   cblas_ztpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpsv(case 1121) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpsv(case 1121) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   double A[] = { 0.262, -0.449 };
   double X[] = { -0.954, -0.093 };
   int incX = 1;
   double x_expected[] = { -0.770378730411, -1.67519103036 };
   cblas_ztpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpsv(case 1122) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpsv(case 1122) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   double A[] = { 0.262, -0.449 };
   double X[] = { -0.954, -0.093 };
   int incX = 1;
   double x_expected[] = { -0.954, -0.093 };
   cblas_ztpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpsv(case 1123) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpsv(case 1123) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   double A[] = { 0.262, -0.449 };
   double X[] = { -0.954, -0.093 };
   int incX = 1;
   double x_expected[] = { -0.770378730411, -1.67519103036 };
   cblas_ztpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpsv(case 1124) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpsv(case 1124) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   double A[] = { 0.262, -0.449 };
   double X[] = { -0.954, -0.093 };
   int incX = 1;
   double x_expected[] = { -0.954, -0.093 };
   cblas_ztpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpsv(case 1125) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpsv(case 1125) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   double A[] = { 0.262, -0.449 };
   double X[] = { -0.954, -0.093 };
   int incX = 1;
   double x_expected[] = { -0.770378730411, -1.67519103036 };
   cblas_ztpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpsv(case 1126) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpsv(case 1126) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 112;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   double A[] = { 0.262, -0.449 };
   double X[] = { -0.954, -0.093 };
   int incX = 1;
   double x_expected[] = { -0.954, -0.093 };
   cblas_ztpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpsv(case 1127) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpsv(case 1127) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   double A[] = { 0.262, -0.449 };
   double X[] = { -0.954, -0.093 };
   int incX = 1;
   double x_expected[] = { -1.07940942478, 1.49486576995 };
   cblas_ztpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpsv(case 1128) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpsv(case 1128) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   double A[] = { 0.262, -0.449 };
   double X[] = { -0.954, -0.093 };
   int incX = 1;
   double x_expected[] = { -0.954, -0.093 };
   cblas_ztpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpsv(case 1129) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpsv(case 1129) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   double A[] = { 0.262, -0.449 };
   double X[] = { -0.954, -0.093 };
   int incX = 1;
   double x_expected[] = { -1.07940942478, 1.49486576995 };
   cblas_ztpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpsv(case 1130) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpsv(case 1130) imag");
     };
   };
  };


  {
   int order = 101;
   int trans = 113;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   double A[] = { 0.262, -0.449 };
   double X[] = { -0.954, -0.093 };
   int incX = 1;
   double x_expected[] = { -0.954, -0.093 };
   cblas_ztpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpsv(case 1131) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpsv(case 1131) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   double A[] = { 0.262, -0.449 };
   double X[] = { -0.954, -0.093 };
   int incX = 1;
   double x_expected[] = { -1.07940942478, 1.49486576995 };
   cblas_ztpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpsv(case 1132) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpsv(case 1132) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 121;
   int diag = 132;
   int N = 1;
   double A[] = { 0.262, -0.449 };
   double X[] = { -0.954, -0.093 };
   int incX = 1;
   double x_expected[] = { -0.954, -0.093 };
   cblas_ztpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpsv(case 1133) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpsv(case 1133) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 122;
   int diag = 131;
   int N = 1;
   double A[] = { 0.262, -0.449 };
   double X[] = { -0.954, -0.093 };
   int incX = 1;
   double x_expected[] = { -1.07940942478, 1.49486576995 };
   cblas_ztpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpsv(case 1134) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpsv(case 1134) imag");
     };
   };
  };


  {
   int order = 102;
   int trans = 113;
   int uplo = 122;
   int diag = 132;
   int N = 1;
   double A[] = { 0.262, -0.449 };
   double X[] = { -0.954, -0.093 };
   int incX = 1;
   double x_expected[] = { -0.954, -0.093 };
   cblas_ztpsv(order, uplo, trans, diag, N, A, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztpsv(case 1135) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztpsv(case 1135) imag");
     };
   };
  };


}
