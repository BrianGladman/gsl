#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_trsv () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 1;
   int lda = 1;
   float A[] = { -0.186 };
   float X[] = { -0.636 };
   int incX = 1;
   float x_expected[] = { 3.41935483871 };
   cblas_strsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strsv(case 896)");
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
   float A[] = { -0.186 };
   float X[] = { -0.636 };
   int incX = 1;
   float x_expected[] = { -0.636 };
   cblas_strsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strsv(case 897)");
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
   float A[] = { -0.186 };
   float X[] = { -0.636 };
   int incX = 1;
   float x_expected[] = { 3.41935483871 };
   cblas_strsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strsv(case 898)");
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
   float A[] = { -0.186 };
   float X[] = { -0.636 };
   int incX = 1;
   float x_expected[] = { -0.636 };
   cblas_strsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strsv(case 899)");
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
   float A[] = { -0.186 };
   float X[] = { -0.636 };
   int incX = 1;
   float x_expected[] = { 3.41935483871 };
   cblas_strsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strsv(case 900)");
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
   float A[] = { -0.186 };
   float X[] = { -0.636 };
   int incX = 1;
   float x_expected[] = { -0.636 };
   cblas_strsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strsv(case 901)");
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
   float A[] = { -0.186 };
   float X[] = { -0.636 };
   int incX = 1;
   float x_expected[] = { 3.41935483871 };
   cblas_strsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strsv(case 902)");
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
   float A[] = { -0.186 };
   float X[] = { -0.636 };
   int incX = 1;
   float x_expected[] = { -0.636 };
   cblas_strsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strsv(case 903)");
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
   float A[] = { 0.609 };
   float X[] = { 0.951 };
   int incX = 1;
   float x_expected[] = { 1.56157635468 };
   cblas_strsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strsv(case 904)");
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
   float A[] = { 0.609 };
   float X[] = { 0.951 };
   int incX = 1;
   float x_expected[] = { 0.951 };
   cblas_strsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strsv(case 905)");
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
   float A[] = { 0.609 };
   float X[] = { 0.951 };
   int incX = 1;
   float x_expected[] = { 1.56157635468 };
   cblas_strsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strsv(case 906)");
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
   float A[] = { 0.609 };
   float X[] = { 0.951 };
   int incX = 1;
   float x_expected[] = { 0.951 };
   cblas_strsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strsv(case 907)");
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
   float A[] = { 0.609 };
   float X[] = { 0.951 };
   int incX = 1;
   float x_expected[] = { 1.56157635468 };
   cblas_strsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strsv(case 908)");
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
   float A[] = { 0.609 };
   float X[] = { 0.951 };
   int incX = 1;
   float x_expected[] = { 0.951 };
   cblas_strsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strsv(case 909)");
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
   float A[] = { 0.609 };
   float X[] = { 0.951 };
   int incX = 1;
   float x_expected[] = { 1.56157635468 };
   cblas_strsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strsv(case 910)");
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
   float A[] = { 0.609 };
   float X[] = { 0.951 };
   int incX = 1;
   float x_expected[] = { 0.951 };
   cblas_strsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "strsv(case 911)");
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
   double A[] = { -0.488 };
   double X[] = { 0.624 };
   int incX = 1;
   double x_expected[] = { -1.27868852459 };
   cblas_dtrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrsv(case 912)");
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
   double A[] = { -0.488 };
   double X[] = { 0.624 };
   int incX = 1;
   double x_expected[] = { 0.624 };
   cblas_dtrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrsv(case 913)");
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
   double A[] = { -0.488 };
   double X[] = { 0.624 };
   int incX = 1;
   double x_expected[] = { -1.27868852459 };
   cblas_dtrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrsv(case 914)");
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
   double A[] = { -0.488 };
   double X[] = { 0.624 };
   int incX = 1;
   double x_expected[] = { 0.624 };
   cblas_dtrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrsv(case 915)");
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
   double A[] = { -0.488 };
   double X[] = { 0.624 };
   int incX = 1;
   double x_expected[] = { -1.27868852459 };
   cblas_dtrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrsv(case 916)");
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
   double A[] = { -0.488 };
   double X[] = { 0.624 };
   int incX = 1;
   double x_expected[] = { 0.624 };
   cblas_dtrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrsv(case 917)");
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
   double A[] = { -0.488 };
   double X[] = { 0.624 };
   int incX = 1;
   double x_expected[] = { -1.27868852459 };
   cblas_dtrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrsv(case 918)");
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
   double A[] = { -0.488 };
   double X[] = { 0.624 };
   int incX = 1;
   double x_expected[] = { 0.624 };
   cblas_dtrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrsv(case 919)");
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
   double A[] = { -0.547 };
   double X[] = { -0.684 };
   int incX = 1;
   double x_expected[] = { 1.25045703839 };
   cblas_dtrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrsv(case 920)");
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
   double A[] = { -0.547 };
   double X[] = { -0.684 };
   int incX = 1;
   double x_expected[] = { -0.684 };
   cblas_dtrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrsv(case 921)");
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
   double A[] = { -0.547 };
   double X[] = { -0.684 };
   int incX = 1;
   double x_expected[] = { 1.25045703839 };
   cblas_dtrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrsv(case 922)");
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
   double A[] = { -0.547 };
   double X[] = { -0.684 };
   int incX = 1;
   double x_expected[] = { -0.684 };
   cblas_dtrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrsv(case 923)");
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
   double A[] = { -0.547 };
   double X[] = { -0.684 };
   int incX = 1;
   double x_expected[] = { 1.25045703839 };
   cblas_dtrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrsv(case 924)");
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
   double A[] = { -0.547 };
   double X[] = { -0.684 };
   int incX = 1;
   double x_expected[] = { -0.684 };
   cblas_dtrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrsv(case 925)");
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
   double A[] = { -0.547 };
   double X[] = { -0.684 };
   int incX = 1;
   double x_expected[] = { 1.25045703839 };
   cblas_dtrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrsv(case 926)");
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
   double A[] = { -0.547 };
   double X[] = { -0.684 };
   int incX = 1;
   double x_expected[] = { -0.684 };
   cblas_dtrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtrsv(case 927)");
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
   float A[] = { 0.314, -0.42 };
   float X[] = { -0.928, -0.722 };
   int incX = 1;
   float x_expected[] = { 0.0430842630438, -2.24173442523 };
   cblas_ctrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrsv(case 928) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrsv(case 928) imag");
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
   float A[] = { 0.314, -0.42 };
   float X[] = { -0.928, -0.722 };
   int incX = 1;
   float x_expected[] = { -0.928, -0.722 };
   cblas_ctrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrsv(case 929) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrsv(case 929) imag");
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
   float A[] = { 0.314, -0.42 };
   float X[] = { -0.928, -0.722 };
   int incX = 1;
   float x_expected[] = { 0.0430842630438, -2.24173442523 };
   cblas_ctrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrsv(case 930) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrsv(case 930) imag");
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
   float A[] = { 0.314, -0.42 };
   float X[] = { -0.928, -0.722 };
   int incX = 1;
   float x_expected[] = { -0.928, -0.722 };
   cblas_ctrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrsv(case 931) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrsv(case 931) imag");
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
   float A[] = { 0.314, -0.42 };
   float X[] = { -0.928, -0.722 };
   int incX = 1;
   float x_expected[] = { 0.0430842630438, -2.24173442523 };
   cblas_ctrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrsv(case 932) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrsv(case 932) imag");
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
   float A[] = { 0.314, -0.42 };
   float X[] = { -0.928, -0.722 };
   int incX = 1;
   float x_expected[] = { -0.928, -0.722 };
   cblas_ctrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrsv(case 933) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrsv(case 933) imag");
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
   float A[] = { 0.314, -0.42 };
   float X[] = { -0.928, -0.722 };
   int incX = 1;
   float x_expected[] = { 0.0430842630438, -2.24173442523 };
   cblas_ctrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrsv(case 934) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrsv(case 934) imag");
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
   float A[] = { 0.314, -0.42 };
   float X[] = { -0.928, -0.722 };
   int incX = 1;
   float x_expected[] = { -0.928, -0.722 };
   cblas_ctrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrsv(case 935) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrsv(case 935) imag");
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
   float A[] = { -0.804, -0.587 };
   float X[] = { -0.797, 0.347 };
   int incX = 1;
   float x_expected[] = { 0.441075293773, -0.753620892344 };
   cblas_ctrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrsv(case 936) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrsv(case 936) imag");
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
   float A[] = { -0.804, -0.587 };
   float X[] = { -0.797, 0.347 };
   int incX = 1;
   float x_expected[] = { -0.797, 0.347 };
   cblas_ctrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrsv(case 937) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrsv(case 937) imag");
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
   float A[] = { -0.804, -0.587 };
   float X[] = { -0.797, 0.347 };
   int incX = 1;
   float x_expected[] = { 0.441075293773, -0.753620892344 };
   cblas_ctrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrsv(case 938) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrsv(case 938) imag");
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
   float A[] = { -0.804, -0.587 };
   float X[] = { -0.797, 0.347 };
   int incX = 1;
   float x_expected[] = { -0.797, 0.347 };
   cblas_ctrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrsv(case 939) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrsv(case 939) imag");
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
   float A[] = { -0.804, -0.587 };
   float X[] = { -0.797, 0.347 };
   int incX = 1;
   float x_expected[] = { 0.441075293773, -0.753620892344 };
   cblas_ctrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrsv(case 940) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrsv(case 940) imag");
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
   float A[] = { -0.804, -0.587 };
   float X[] = { -0.797, 0.347 };
   int incX = 1;
   float x_expected[] = { -0.797, 0.347 };
   cblas_ctrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrsv(case 941) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrsv(case 941) imag");
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
   float A[] = { -0.804, -0.587 };
   float X[] = { -0.797, 0.347 };
   int incX = 1;
   float x_expected[] = { 0.441075293773, -0.753620892344 };
   cblas_ctrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrsv(case 942) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrsv(case 942) imag");
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
   float A[] = { -0.804, -0.587 };
   float X[] = { -0.797, 0.347 };
   int incX = 1;
   float x_expected[] = { -0.797, 0.347 };
   cblas_ctrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrsv(case 943) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrsv(case 943) imag");
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
   float A[] = { -0.753, -0.144 };
   float X[] = { 0.26, -0.696 };
   int incX = 1;
   float x_expected[] = { -0.503626572748, 0.827991731108 };
   cblas_ctrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrsv(case 944) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrsv(case 944) imag");
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
   float A[] = { -0.753, -0.144 };
   float X[] = { 0.26, -0.696 };
   int incX = 1;
   float x_expected[] = { 0.26, -0.696 };
   cblas_ctrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrsv(case 945) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrsv(case 945) imag");
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
   float A[] = { -0.753, -0.144 };
   float X[] = { 0.26, -0.696 };
   int incX = 1;
   float x_expected[] = { -0.503626572748, 0.827991731108 };
   cblas_ctrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrsv(case 946) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrsv(case 946) imag");
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
   float A[] = { -0.753, -0.144 };
   float X[] = { 0.26, -0.696 };
   int incX = 1;
   float x_expected[] = { 0.26, -0.696 };
   cblas_ctrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrsv(case 947) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrsv(case 947) imag");
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
   float A[] = { -0.753, -0.144 };
   float X[] = { 0.26, -0.696 };
   int incX = 1;
   float x_expected[] = { -0.503626572748, 0.827991731108 };
   cblas_ctrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrsv(case 948) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrsv(case 948) imag");
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
   float A[] = { -0.753, -0.144 };
   float X[] = { 0.26, -0.696 };
   int incX = 1;
   float x_expected[] = { 0.26, -0.696 };
   cblas_ctrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrsv(case 949) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrsv(case 949) imag");
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
   float A[] = { -0.753, -0.144 };
   float X[] = { 0.26, -0.696 };
   int incX = 1;
   float x_expected[] = { -0.503626572748, 0.827991731108 };
   cblas_ctrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrsv(case 950) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrsv(case 950) imag");
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
   float A[] = { -0.753, -0.144 };
   float X[] = { 0.26, -0.696 };
   int incX = 1;
   float x_expected[] = { 0.26, -0.696 };
   cblas_ctrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctrsv(case 951) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctrsv(case 951) imag");
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
   double A[] = { 0.142, 0.549 };
   double X[] = { -0.805, -0.965 };
   int incX = 1;
   double x_expected[] = { -2.00300094849, 0.948221976894 };
   cblas_ztrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrsv(case 952) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrsv(case 952) imag");
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
   double A[] = { 0.142, 0.549 };
   double X[] = { -0.805, -0.965 };
   int incX = 1;
   double x_expected[] = { -0.805, -0.965 };
   cblas_ztrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrsv(case 953) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrsv(case 953) imag");
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
   double A[] = { 0.142, 0.549 };
   double X[] = { -0.805, -0.965 };
   int incX = 1;
   double x_expected[] = { -2.00300094849, 0.948221976894 };
   cblas_ztrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrsv(case 954) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrsv(case 954) imag");
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
   double A[] = { 0.142, 0.549 };
   double X[] = { -0.805, -0.965 };
   int incX = 1;
   double x_expected[] = { -0.805, -0.965 };
   cblas_ztrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrsv(case 955) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrsv(case 955) imag");
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
   double A[] = { 0.142, 0.549 };
   double X[] = { -0.805, -0.965 };
   int incX = 1;
   double x_expected[] = { -2.00300094849, 0.948221976894 };
   cblas_ztrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrsv(case 956) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrsv(case 956) imag");
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
   double A[] = { 0.142, 0.549 };
   double X[] = { -0.805, -0.965 };
   int incX = 1;
   double x_expected[] = { -0.805, -0.965 };
   cblas_ztrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrsv(case 957) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrsv(case 957) imag");
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
   double A[] = { 0.142, 0.549 };
   double X[] = { -0.805, -0.965 };
   int incX = 1;
   double x_expected[] = { -2.00300094849, 0.948221976894 };
   cblas_ztrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrsv(case 958) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrsv(case 958) imag");
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
   double A[] = { 0.142, 0.549 };
   double X[] = { -0.805, -0.965 };
   int incX = 1;
   double x_expected[] = { -0.805, -0.965 };
   cblas_ztrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrsv(case 959) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrsv(case 959) imag");
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
   double A[] = { 0.537, -0.059 };
   double X[] = { 0.751, 0.108 };
   int incX = 1;
   double x_expected[] = { 1.35999657358, 0.350539660785 };
   cblas_ztrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrsv(case 960) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrsv(case 960) imag");
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
   double A[] = { 0.537, -0.059 };
   double X[] = { 0.751, 0.108 };
   int incX = 1;
   double x_expected[] = { 0.751, 0.108 };
   cblas_ztrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrsv(case 961) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrsv(case 961) imag");
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
   double A[] = { 0.537, -0.059 };
   double X[] = { 0.751, 0.108 };
   int incX = 1;
   double x_expected[] = { 1.35999657358, 0.350539660785 };
   cblas_ztrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrsv(case 962) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrsv(case 962) imag");
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
   double A[] = { 0.537, -0.059 };
   double X[] = { 0.751, 0.108 };
   int incX = 1;
   double x_expected[] = { 0.751, 0.108 };
   cblas_ztrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrsv(case 963) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrsv(case 963) imag");
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
   double A[] = { 0.537, -0.059 };
   double X[] = { 0.751, 0.108 };
   int incX = 1;
   double x_expected[] = { 1.35999657358, 0.350539660785 };
   cblas_ztrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrsv(case 964) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrsv(case 964) imag");
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
   double A[] = { 0.537, -0.059 };
   double X[] = { 0.751, 0.108 };
   int incX = 1;
   double x_expected[] = { 0.751, 0.108 };
   cblas_ztrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrsv(case 965) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrsv(case 965) imag");
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
   double A[] = { 0.537, -0.059 };
   double X[] = { 0.751, 0.108 };
   int incX = 1;
   double x_expected[] = { 1.35999657358, 0.350539660785 };
   cblas_ztrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrsv(case 966) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrsv(case 966) imag");
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
   double A[] = { 0.537, -0.059 };
   double X[] = { 0.751, 0.108 };
   int incX = 1;
   double x_expected[] = { 0.751, 0.108 };
   cblas_ztrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrsv(case 967) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrsv(case 967) imag");
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
   double A[] = { -0.709, -0.047 };
   double X[] = { 0.672, 0.554 };
   int incX = 1;
   double x_expected[] = { -0.892095307889, -0.840519717166 };
   cblas_ztrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrsv(case 968) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrsv(case 968) imag");
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
   double A[] = { -0.709, -0.047 };
   double X[] = { 0.672, 0.554 };
   int incX = 1;
   double x_expected[] = { 0.672, 0.554 };
   cblas_ztrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrsv(case 969) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrsv(case 969) imag");
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
   double A[] = { -0.709, -0.047 };
   double X[] = { 0.672, 0.554 };
   int incX = 1;
   double x_expected[] = { -0.892095307889, -0.840519717166 };
   cblas_ztrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrsv(case 970) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrsv(case 970) imag");
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
   double A[] = { -0.709, -0.047 };
   double X[] = { 0.672, 0.554 };
   int incX = 1;
   double x_expected[] = { 0.672, 0.554 };
   cblas_ztrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrsv(case 971) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrsv(case 971) imag");
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
   double A[] = { -0.709, -0.047 };
   double X[] = { 0.672, 0.554 };
   int incX = 1;
   double x_expected[] = { -0.892095307889, -0.840519717166 };
   cblas_ztrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrsv(case 972) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrsv(case 972) imag");
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
   double A[] = { -0.709, -0.047 };
   double X[] = { 0.672, 0.554 };
   int incX = 1;
   double x_expected[] = { 0.672, 0.554 };
   cblas_ztrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrsv(case 973) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrsv(case 973) imag");
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
   double A[] = { -0.709, -0.047 };
   double X[] = { 0.672, 0.554 };
   int incX = 1;
   double x_expected[] = { -0.892095307889, -0.840519717166 };
   cblas_ztrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrsv(case 974) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrsv(case 974) imag");
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
   double A[] = { -0.709, -0.047 };
   double X[] = { 0.672, 0.554 };
   int incX = 1;
   double x_expected[] = { 0.672, 0.554 };
   cblas_ztrsv(order, uplo, trans, diag, N, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztrsv(case 975) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztrsv(case 975) imag");
     };
   };
  };


}
