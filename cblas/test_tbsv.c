#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

#include "tests.h"

void
test_tbsv (void) {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int trans = 111;
   int uplo = 121;
   int diag = 131;
   int N = 3;
   int K = 1;
   int lda = 3;
   float A[] = { -0.681, 0.209, 0.436, -0.369, 0.786, -0.84, 0.86, -0.233, 0.734 };
   float X[] = { -0.305, 0.61, -0.831 };
   int incX = -1;
   float x_expected[] = { -0.354651, -2.40855, 0.481076 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 1230)");
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
   float A[] = { -0.681, 0.209, 0.436, -0.369, 0.786, -0.84, 0.86, -0.233, 0.734 };
   float X[] = { -0.305, 0.61, -0.831 };
   int incX = -1;
   float x_expected[] = { -0.305, 0.84973, -1.00859 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 1231)");
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
   float A[] = { -0.681, 0.209, 0.436, -0.369, 0.786, -0.84, 0.86, -0.233, 0.734 };
   float X[] = { -0.305, 0.61, -0.831 };
   int incX = -1;
   float x_expected[] = { -2.71619, -1.09055, -3.97608 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 1232)");
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
   float A[] = { -0.681, 0.209, 0.436, -0.369, 0.786, -0.84, 0.86, -0.233, 0.734 };
   float X[] = { -0.305, 0.61, -0.831 };
   int incX = -1;
   float x_expected[] = { -0.56589, 0.303361, -0.831 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 1233)");
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
   float A[] = { -0.681, 0.209, 0.436, -0.369, 0.786, -0.84, 0.86, -0.233, 0.734 };
   float X[] = { -0.305, 0.61, -0.831 };
   int incX = -1;
   float x_expected[] = { 1.30901, -0.656172, -5.13458 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 1234)");
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
   float A[] = { -0.681, 0.209, 0.436, -0.369, 0.786, -0.84, 0.86, -0.233, 0.734 };
   float X[] = { -0.305, 0.61, -0.831 };
   int incX = -1;
   float x_expected[] = { -0.305, 0.8723, -0.509121 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 1235)");
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
   float A[] = { -0.681, 0.209, 0.436, -0.369, 0.786, -0.84, 0.86, -0.233, 0.734 };
   float X[] = { -0.305, 0.61, -0.831 };
   int incX = -1;
   float x_expected[] = { 0.524539, -0.961964, 1.22026 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 1236)");
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
   float A[] = { -0.681, 0.209, 0.436, -0.369, 0.786, -0.84, 0.86, -0.233, 0.734 };
   float X[] = { -0.305, 0.61, -0.831 };
   int incX = -1;
   float x_expected[] = { -0.920972, 0.783679, -0.831 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 1237)");
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
   float A[] = { 0.022, 0.795, -0.389, -0.205, -0.121, 0.323, 0.133, 0.679, 0.742 };
   float X[] = { 0.144, 0.635, 0.116 };
   int incX = -1;
   float x_expected[] = { 16.8676, 17.3503, 5.27273 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 1238)");
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
   float A[] = { 0.022, 0.795, -0.389, -0.205, -0.121, 0.323, 0.133, 0.679, 0.742 };
   float X[] = { 0.144, 0.635, 0.116 };
   int incX = -1;
   float x_expected[] = { 0.209676, 0.54278, 0.116 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 1239)");
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
   float A[] = { 0.022, 0.795, -0.389, -0.205, -0.121, 0.323, 0.133, 0.679, 0.742 };
   float X[] = { 0.144, 0.635, 0.116 };
   int incX = -1;
   float x_expected[] = { 0.212077, -5.01482, -1.14722 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 1240)");
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
   float A[] = { 0.022, 0.795, -0.389, -0.205, -0.121, 0.323, 0.133, 0.679, 0.742 };
   float X[] = { 0.144, 0.635, 0.116 };
   int incX = -1;
   float x_expected[] = { 0.144, 0.615848, 0.242249 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 1241)");
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
   float A[] = { 0.022, 0.795, -0.389, -0.205, -0.121, 0.323, 0.133, 0.679, 0.742 };
   float X[] = { 0.144, 0.635, 0.116 };
   int incX = -1;
   float x_expected[] = { 1.28844, -5.49514, 0.145912 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 1242)");
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
   float A[] = { 0.022, 0.795, -0.389, -0.205, -0.121, 0.323, 0.133, 0.679, 0.742 };
   float X[] = { 0.144, 0.635, 0.116 };
   int incX = -1;
   float x_expected[] = { 0.0563823, 0.65878, 0.116 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 1243)");
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
   float A[] = { 0.022, 0.795, -0.389, -0.205, -0.121, 0.323, 0.133, 0.679, 0.742 };
   float X[] = { 0.144, 0.635, 0.116 };
   int incX = -1;
   float x_expected[] = { 1.08271, -3.73662, 140.301 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 1244)");
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
   float A[] = { 0.022, 0.795, -0.389, -0.205, -0.121, 0.323, 0.133, 0.679, 0.742 };
   float X[] = { 0.144, 0.635, 0.116 };
   int incX = -1;
   float x_expected[] = { 0.144, 0.652424, -0.402677 };
   cblas_stbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "stbsv(case 1245)");
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
   double A[] = { 0.619, -0.443, 0.957, -0.633, -0.698, 0.783, -0.343, -0.603, 0.735 };
   double X[] = { 0.332, 0.588, 0.252 };
   int incX = -1;
   double x_expected[] = { -0.967930029155, 0.138412575592, 0.506166027443 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 1246)");
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
   double A[] = { 0.619, -0.443, 0.957, -0.633, -0.698, 0.783, -0.343, -0.603, 0.735 };
   double X[] = { 0.332, 0.588, 0.252 };
   int incX = -1;
   double x_expected[] = { 0.332, 0.819736, 0.615143048 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 1247)");
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
   double A[] = { 0.619, -0.443, 0.957, -0.633, -0.698, 0.783, -0.343, -0.603, 0.735 };
   double X[] = { 0.332, 0.588, 0.252 };
   int incX = -1;
   double x_expected[] = { -0.364842154056, -0.326531140246, -0.568848758465 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 1248)");
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
   double A[] = { 0.619, -0.443, 0.957, -0.633, -0.698, 0.783, -0.343, -0.603, 0.735 };
   double X[] = { 0.332, 0.588, 0.252 };
   int incX = -1;
   double x_expected[] = { 0.588397988, 0.747516, 0.252 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 1249)");
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
   double A[] = { 0.619, -0.443, 0.957, -0.633, -0.698, 0.783, -0.343, -0.603, 0.735 };
   double X[] = { 0.332, 0.588, 0.252 };
   int incX = -1;
   double x_expected[] = { -0.550580431177, -0.571849444278, 0.248263427151 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 1250)");
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
   double A[] = { 0.619, -0.443, 0.957, -0.633, -0.698, 0.783, -0.343, -0.603, 0.735 };
   double X[] = { 0.332, 0.588, 0.252 };
   int incX = -1;
   double x_expected[] = { 0.332, 0.701876, 0.696287508 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 1251)");
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
   double A[] = { 0.619, -0.443, 0.957, -0.633, -0.698, 0.783, -0.343, -0.603, 0.735 };
   double X[] = { 0.332, 0.588, 0.252 };
   int incX = -1;
   double x_expected[] = { 1.50217883761, -1.21382140588, 0.407108239095 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 1252)");
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
   double A[] = { 0.619, -0.443, 0.957, -0.633, -0.698, 0.783, -0.343, -0.603, 0.735 };
   double X[] = { 0.332, 0.588, 0.252 };
   int incX = -1;
   double x_expected[] = { 0.820345928, 0.699636, 0.252 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 1253)");
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
   double A[] = { -0.199, 0.303, -0.705, -0.013, -0.678, 0.547, 0.756, -0.177, -0.079 };
   double X[] = { 0.58, 0.558, -0.54 };
   int incX = -1;
   double x_expected[] = { 18.994209959, 20.323927329, 2.7135678392 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 1254)");
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
   double A[] = { -0.199, 0.303, -0.705, -0.013, -0.678, 0.547, 0.756, -0.177, -0.079 };
   double X[] = { 0.58, 0.558, -0.54 };
   int incX = -1;
   double x_expected[] = { 1.06925836, 0.72162, -0.54 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 1255)");
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
   double A[] = { -0.199, 0.303, -0.705, -0.013, -0.678, 0.547, 0.756, -0.177, -0.079 };
   double X[] = { 0.58, 0.558, -0.54 };
   int incX = -1;
   double x_expected[] = { -3.27683615819, -4.47682615869, -1.97425326753 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 1256)");
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
   double A[] = { -0.199, 0.303, -0.705, -0.013, -0.678, 0.547, 0.756, -0.177, -0.079 };
   double X[] = { 0.58, 0.558, -0.54 };
   int incX = -1;
   double x_expected[] = { 0.58, 0.11952, -0.53844624 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 1257)");
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
   double A[] = { -0.199, 0.303, -0.705, -0.013, -0.678, 0.547, 0.756, -0.177, -0.079 };
   double X[] = { 0.58, 0.558, -0.54 };
   int incX = -1;
   double x_expected[] = { -6.6461072986, -0.788837290809, -1.78217821782 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 1258)");
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
   double A[] = { -0.199, 0.303, -0.705, -0.013, -0.678, 0.547, 0.756, -0.177, -0.079 };
   double X[] = { 0.58, 0.558, -0.54 };
   int incX = -1;
   double x_expected[] = { 0.16345912, 0.55098, -0.54 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 1259)");
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
   double A[] = { -0.199, 0.303, -0.705, -0.013, -0.678, 0.547, 0.756, -0.177, -0.079 };
   double X[] = { 0.58, 0.558, -0.54 };
   int incX = -1;
   double x_expected[] = { 0.767195767196, -82.9352869353, -123.564783625 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 1260)");
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
   double A[] = { -0.199, 0.303, -0.705, -0.013, -0.678, 0.547, 0.756, -0.177, -0.079 };
   double X[] = { 0.58, 0.558, -0.54 };
   int incX = -1;
   double x_expected[] = { 0.58, 0.95124, -0.82822572 };
   cblas_dtbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "dtbsv(case 1261)");
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
   float A[] = { -0.975, -0.667, 0.813, -0.962, -0.961, 0.226, -0.503, 0.809, 0.81, -0.162, -0.027, -0.044, 0.212, 0.563, 0.446, -0.392, 0.798, -0.07 };
   float X[] = { 0.11, 0.787, -0.826, 0.809, -0.437, 0.592 };
   int incX = -1;
   float x_expected[] = { 1.28871, 0.289887, 1.76043, 1.27481, 1.56506, -2.35181 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1262) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1262) imag");
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
   float A[] = { -0.975, -0.667, 0.813, -0.962, -0.961, 0.226, -0.503, 0.809, 0.81, -0.162, -0.027, -0.044, 0.212, 0.563, 0.446, -0.392, 0.798, -0.07 };
   float X[] = { 0.11, 0.787, -0.826, 0.809, -0.437, 0.592 };
   int incX = -1;
   float x_expected[] = { 0.11, 0.787, -1.04259, 0.18935, 0.228474, -0.564917 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1263) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1263) imag");
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
   float A[] = { -0.975, -0.667, 0.813, -0.962, -0.961, 0.226, -0.503, 0.809, 0.81, -0.162, -0.027, -0.044, 0.212, 0.563, 0.446, -0.392, 0.798, -0.07 };
   float X[] = { 0.11, 0.787, -0.826, 0.809, -0.437, 0.592 };
   int incX = -1;
   float x_expected[] = { -0.0906249, 3.09442, -1.60036, 1.28475, -0.582941, 0.0383898 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1264) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1264) imag");
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
   float A[] = { -0.975, -0.667, 0.813, -0.962, -0.961, 0.226, -0.503, 0.809, 0.81, -0.162, -0.027, -0.044, 0.212, 0.563, 0.446, -0.392, 0.798, -0.07 };
   float X[] = { 0.11, 0.787, -0.826, 0.809, -0.437, 0.592 };
   int incX = -1;
   float x_expected[] = { 1.05233, 0.79657, -0.566883, 1.46031, -0.437, 0.592 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1265) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1265) imag");
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
   float A[] = { -0.975, -0.667, 0.813, -0.962, -0.961, 0.226, -0.503, 0.809, 0.81, -0.162, -0.027, -0.044, 0.212, 0.563, 0.446, -0.392, 0.798, -0.07 };
   float X[] = { 0.11, 0.787, -0.826, 0.809, -0.437, 0.592 };
   int incX = -1;
   float x_expected[] = { -0.735844, 1.11782, -0.28244, 1.16117, -0.66707, 0.938302 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1266) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1266) imag");
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
   float A[] = { -0.975, -0.667, 0.813, -0.962, -0.961, 0.226, -0.503, 0.809, 0.81, -0.162, -0.027, -0.044, 0.212, 0.563, 0.446, -0.392, 0.798, -0.07 };
   float X[] = { 0.11, 0.787, -0.826, 0.809, -0.437, 0.592 };
   int incX = -1;
   float x_expected[] = { 0.11, 0.787, -0.406239, 0.580226, -0.171935, 1.2125 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1267) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1267) imag");
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
   float A[] = { -0.975, -0.667, 0.813, -0.962, -0.961, 0.226, -0.503, 0.809, 0.81, -0.162, -0.027, -0.044, 0.212, 0.563, 0.446, -0.392, 0.798, -0.07 };
   float X[] = { 0.11, 0.787, -0.826, 0.809, -0.437, 0.592 };
   int incX = -1;
   float x_expected[] = { 1.70081, 2.20477, 1.32753, -0.522112, 0.0223652, -0.62248 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1268) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1268) imag");
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
   float A[] = { -0.975, -0.667, 0.813, -0.962, -0.961, 0.226, -0.503, 0.809, 0.81, -0.162, -0.027, -0.044, 0.212, 0.563, 0.446, -0.392, 0.798, -0.07 };
   float X[] = { 0.11, 0.787, -0.826, 0.809, -0.437, 0.592 };
   int incX = -1;
   float x_expected[] = { 0.967596, 0.693563, -1.04022, -0.09269, -0.437, 0.592 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1269) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1269) imag");
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
   float A[] = { 0.33, -0.236, 0.267, -0.139, 0.25, 0.509, 0.86, -0.089, -0.018, -0.847, 0.424, -0.573, 0.097, -0.663, 0.65, -0.811, 0.283, 0.032 };
   float X[] = { -0.509, 0.608, 0.021, -0.848, -0.134, 0.179 };
   int incX = -1;
   float x_expected[] = { -1.11985, 0.801655, 0.273814, -1.09438, -0.52531, 0.166748 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1270) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1270) imag");
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
   float A[] = { 0.33, -0.236, 0.267, -0.139, 0.25, 0.509, 0.86, -0.089, -0.018, -0.847, 0.424, -0.573, 0.097, -0.663, 0.65, -0.811, 0.283, 0.032 };
   float X[] = { -0.509, 0.608, 0.021, -0.848, -0.134, 0.179 };
   int incX = -1;
   float x_expected[] = { 0.266087, 0.618557, 0.031897, -0.914419, -0.134, 0.179 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1271) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1271) imag");
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
   float A[] = { 0.33, -0.236, 0.267, -0.139, 0.25, 0.509, 0.86, -0.089, -0.018, -0.847, 0.424, -0.573, 0.097, -0.663, 0.65, -0.811, 0.283, 0.032 };
   float X[] = { -0.509, 0.608, 0.021, -0.848, -0.134, 0.179 };
   int incX = -1;
   float x_expected[] = { -0.762749, -0.016292, 1.59299, 0.158751, -4.75603, -1.78591 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1272) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1272) imag");
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
   float A[] = { 0.33, -0.236, 0.267, -0.139, 0.25, 0.509, 0.86, -0.089, -0.018, -0.847, 0.424, -0.573, 0.097, -0.663, 0.65, -0.811, 0.283, 0.032 };
   float X[] = { -0.509, 0.608, 0.021, -0.848, -0.134, 0.179 };
   int incX = -1;
   float x_expected[] = { -0.509, 0.608, -0.332731, -1.24444, 0.262904, 1.21961 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1273) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1273) imag");
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
   float A[] = { 0.33, -0.236, 0.267, -0.139, 0.25, 0.509, 0.86, -0.089, -0.018, -0.847, 0.424, -0.573, 0.097, -0.663, 0.65, -0.811, 0.283, 0.032 };
   float X[] = { -0.509, 0.608, 0.021, -0.848, -0.134, 0.179 };
   int incX = -1;
   float x_expected[] = { -1.76046, 0.0455463, 1.38348, 0.700097, -0.669451, 0.321896 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1274) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1274) imag");
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
   float A[] = { 0.33, -0.236, 0.267, -0.139, 0.25, 0.509, 0.86, -0.089, -0.018, -0.847, 0.424, -0.573, 0.097, -0.663, 0.65, -0.811, 0.283, 0.032 };
   float X[] = { -0.509, 0.608, 0.021, -0.848, -0.134, 0.179 };
   int incX = -1;
   float x_expected[] = { 0.151523, 0.78611, 0.120309, -1.01387, -0.134, 0.179 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1275) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1275) imag");
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
   float A[] = { 0.33, -0.236, 0.267, -0.139, 0.25, 0.509, 0.86, -0.089, -0.018, -0.847, 0.424, -0.573, 0.097, -0.663, 0.65, -0.811, 0.283, 0.032 };
   float X[] = { -0.509, 0.608, 0.021, -0.848, -0.134, 0.179 };
   int incX = -1;
   float x_expected[] = { -1.00779, -0.620278, 0.81164, -1.90759, -1.32022, 1.48356 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1276) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1276) imag");
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
   float A[] = { 0.33, -0.236, 0.267, -0.139, 0.25, 0.509, 0.86, -0.089, -0.018, -0.847, 0.424, -0.573, 0.097, -0.663, 0.65, -0.811, 0.283, 0.032 };
   float X[] = { -0.509, 0.608, 0.021, -0.848, -0.134, 0.179 };
   int incX = -1;
   float x_expected[] = { -0.509, 0.608, -0.503138, -1.26818, 0.176615, 0.447668 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1277) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1277) imag");
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
   float A[] = { 0.041, -0.61, 0.099, -0.393, 0.357, -0.984, -0.576, -0.342, -0.903, -0.083, -0.157, -0.694, 0.768, 0.688, 0.203, -0.079, 0.298, -0.424 };
   float X[] = { -0.037, -0.599, 0.959, -0.499, 0.296, 0.034 };
   int incX = -1;
   float x_expected[] = { -0.613838, -1.13321, -1.34847, 0.0432903, 0.0879552, -0.479334 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1278) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1278) imag");
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
   float A[] = { 0.041, -0.61, 0.099, -0.393, 0.357, -0.984, -0.576, -0.342, -0.903, -0.083, -0.157, -0.694, 0.768, 0.688, 0.203, -0.079, 0.298, -0.424 };
   float X[] = { -0.037, -0.599, 0.959, -0.499, 0.296, 0.034 };
   int incX = -1;
   float x_expected[] = { 0.76323, -1.23595, 0.943058, -0.618694, 0.296, 0.034 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1279) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1279) imag");
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
   float A[] = { 0.041, -0.61, 0.099, -0.393, 0.357, -0.984, -0.576, -0.342, -0.903, -0.083, -0.157, -0.694, 0.768, 0.688, 0.203, -0.079, 0.298, -0.424 };
   float X[] = { -0.037, -0.599, 0.959, -0.499, 0.296, 0.034 };
   int incX = -1;
   float x_expected[] = { -1.15557, -2.50103, -3.85402, -1.04833, 0.414582, 5.91218 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1280) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1280) imag");
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
   float A[] = { 0.041, -0.61, 0.099, -0.393, 0.357, -0.984, -0.576, -0.342, -0.903, -0.083, -0.157, -0.694, 0.768, 0.688, 0.203, -0.079, 0.298, -0.424 };
   float X[] = { -0.037, -0.599, 0.959, -0.499, 0.296, 0.034 };
   int incX = -1;
   float x_expected[] = { -0.037, -0.599, 1.39953, -0.064424, 1.0801, -0.481747 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1281) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1281) imag");
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
   float A[] = { 0.041, -0.61, 0.099, -0.393, 0.357, -0.984, -0.576, -0.342, -0.903, -0.083, -0.157, -0.694, 0.768, 0.688, 0.203, -0.079, 0.298, -0.424 };
   float X[] = { -0.037, -0.599, 0.959, -0.499, 0.296, 0.034 };
   int incX = -1;
   float x_expected[] = { -3.0802, -9.09377, -1.05845, 0.99239, 0.259763, -0.687744 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1282) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1282) imag");
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
   float A[] = { 0.041, -0.61, 0.099, -0.393, 0.357, -0.984, -0.576, -0.342, -0.903, -0.083, -0.157, -0.694, 0.768, 0.688, 0.203, -0.079, 0.298, -0.424 };
   float X[] = { -0.037, -0.599, 0.959, -0.499, 0.296, 0.034 };
   int incX = -1;
   float x_expected[] = { -0.513897, 0.632031, 1.14112, -0.580648, 0.296, 0.034 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1283) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1283) imag");
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
   float A[] = { 0.041, -0.61, 0.099, -0.393, 0.357, -0.984, -0.576, -0.342, -0.903, -0.083, -0.157, -0.694, 0.768, 0.688, 0.203, -0.079, 0.298, -0.424 };
   float X[] = { -0.037, -0.599, 0.959, -0.499, 0.296, 0.034 };
   int incX = -1;
   float x_expected[] = { 0.360899, -0.456643, -2.31803, 0.257877, 1.56928, -0.922115 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1284) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1284) imag");
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
   float A[] = { 0.041, -0.61, 0.099, -0.393, 0.357, -0.984, -0.576, -0.342, -0.903, -0.083, -0.157, -0.694, 0.768, 0.688, 0.203, -0.079, 0.298, -0.424 };
   float X[] = { -0.037, -0.599, 0.959, -0.499, 0.296, 0.034 };
   int incX = -1;
   float x_expected[] = { -0.037, -0.599, 0.875872, -1.03683, -0.198184, -0.207572 };
   cblas_ctbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], flteps, "ctbsv(case 1285) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], flteps, "ctbsv(case 1285) imag");
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
   double A[] = { 0.474, 0.715, 0.061, 0.532, 0.004, -0.318, 0.37, -0.692, -0.166, 0.039, -0.946, 0.857, -0.922, -0.491, 0.012, -0.217, -0.674, -0.429 };
   double X[] = { -0.123, 0.122, 0.981, 0.321, 0.942, 0.98 };
   int incX = -1;
   double x_expected[] = { 0.0490338308139, -0.158433417494, 0.261604043488, 1.28058846321, 1.77633350191, -1.07039599422 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1286) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1286) imag");
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
   double A[] = { 0.474, 0.715, 0.061, 0.532, 0.004, -0.318, 0.37, -0.692, -0.166, 0.039, -0.946, 0.857, -0.922, -0.491, 0.012, -0.217, -0.674, -0.429 };
   double X[] = { -0.123, 0.122, 0.981, 0.321, 0.942, 0.98 };
   int incX = -1;
   double x_expected[] = { -0.123, 0.122, 0.96534, 0.346049, 1.067212328, 0.445330131 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1287) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1287) imag");
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
   double A[] = { 0.474, 0.715, 0.061, 0.532, 0.004, -0.318, 0.37, -0.692, -0.166, 0.039, -0.946, 0.857, -0.922, -0.491, 0.012, -0.217, -0.674, -0.429 };
   double X[] = { -0.123, 0.122, 0.981, 0.321, 0.942, 0.98 };
   int incX = -1;
   double x_expected[] = { 72.7437666278, 10.4206532927, -4.34946941374, -14.8012581742, 2.01859491883, -1.53922125931 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1288) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1288) imag");
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
   double A[] = { 0.474, 0.715, 0.061, 0.532, 0.004, -0.318, 0.37, -0.692, -0.166, 0.039, -0.946, 0.857, -0.922, -0.491, 0.012, -0.217, -0.674, -0.429 };
   double X[] = { -0.123, 0.122, 0.981, 0.321, 0.942, 0.98 };
   int incX = -1;
   double x_expected[] = { -0.464775024, 0.662224708, -0.0457, 0.610264, 0.942, 0.98 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1289) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1289) imag");
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
   double A[] = { 0.474, 0.715, 0.061, 0.532, 0.004, -0.318, 0.37, -0.692, -0.166, 0.039, -0.946, 0.857, -0.922, -0.491, 0.012, -0.217, -0.674, -0.429 };
   double X[] = { -0.123, 0.122, 0.981, 0.321, 0.942, 0.98 };
   int incX = -1;
   double x_expected[] = { -0.591747295323, -0.534096923761, -4.60251824353, 1.70172936273, -4.94687072873, -3.32536493524 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1290) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1290) imag");
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
   double A[] = { 0.474, 0.715, 0.061, 0.532, 0.004, -0.318, 0.37, -0.692, -0.166, 0.039, -0.946, 0.857, -0.922, -0.491, 0.012, -0.217, -0.674, -0.429 };
   double X[] = { -0.123, 0.122, 0.981, 0.321, 0.942, 0.98 };
   int incX = -1;
   double x_expected[] = { -0.123, 0.122, 0.807692, 0.373091, 0.384974988, 1.400879194 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1291) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1291) imag");
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
   double A[] = { 0.474, 0.715, 0.061, 0.532, 0.004, -0.318, 0.37, -0.692, -0.166, 0.039, -0.946, 0.857, -0.922, -0.491, 0.012, -0.217, -0.674, -0.429 };
   double X[] = { -0.123, 0.122, 0.981, 0.321, 0.942, 0.98 };
   int incX = -1;
   double x_expected[] = { -0.129998778267, -0.116630230861, 0.993340886904, 0.530739563688, 1.55891621291, -0.284019181928 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1292) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1292) imag");
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
   double A[] = { 0.474, 0.715, 0.061, 0.532, 0.004, -0.318, 0.37, -0.692, -0.166, 0.039, -0.946, 0.857, -0.922, -0.491, 0.012, -0.217, -0.674, -0.429 };
   double X[] = { -0.123, 0.122, 0.981, 0.321, 0.942, 0.98 };
   int incX = -1;
   double x_expected[] = { 0.107496032, 0.025821594, 1.444898, -0.239924, 0.942, 0.98 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1293) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1293) imag");
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
   double A[] = { -0.872, -0.841, 0.108, -0.744, 0.231, -0.513, -0.973, 0.087, 0.348, 0.196, 0.447, 0.307, 0.632, -0.949, 0.322, 0.277, 0.282, 0.831 };
   double X[] = { -0.373, 0.566, 0.92, 0.627, 0.293, -0.434 };
   int incX = -1;
   double x_expected[] = { -0.825842176606, 0.212941473892, -0.548817434511, -0.703261551538, 0.0746069436827, 0.425751789407 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1294) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1294) imag");
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
   double A[] = { -0.872, -0.841, 0.108, -0.744, 0.231, -0.513, -0.973, 0.087, 0.348, 0.196, 0.447, 0.307, 0.632, -0.949, 0.322, 0.277, 0.282, 0.831 };
   double X[] = { -0.373, 0.566, 0.92, 0.627, 0.293, -0.434 };
   int incX = -1;
   double x_expected[] = { -0.619710352, 0.018225936, 1.211252, 0.891864, 0.293, -0.434 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1295) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1295) imag");
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
   double A[] = { -0.872, -0.841, 0.108, -0.744, 0.231, -0.513, -0.973, 0.087, 0.348, 0.196, 0.447, 0.307, 0.632, -0.949, 0.322, 0.277, 0.282, 0.831 };
   double X[] = { -0.373, 0.566, 0.92, 0.627, 0.293, -0.434 };
   int incX = -1;
   double x_expected[] = { 0.203289119964, 1.58288482537, -1.7720160159, 0.479463518178, -0.511241930019, -1.79333888299 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1296) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1296) imag");
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
   double A[] = { -0.872, -0.841, 0.108, -0.744, 0.231, -0.513, -0.973, 0.087, 0.348, 0.196, 0.447, 0.307, 0.632, -0.949, 0.322, 0.277, 0.282, 0.831 };
   double X[] = { -0.373, 0.566, 0.92, 0.627, 0.293, -0.434 };
   int incX = -1;
   double x_expected[] = { -0.373, 0.566, 0.618602, -0.084689, 0.887531803, -0.570220771 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1297) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1297) imag");
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
   double A[] = { -0.872, -0.841, 0.108, -0.744, 0.231, -0.513, -0.973, 0.087, 0.348, 0.196, 0.447, 0.307, 0.632, -0.949, 0.322, 0.277, 0.282, 0.831 };
   double X[] = { -0.373, 0.566, 0.92, 0.627, 0.293, -0.434 };
   int incX = -1;
   double x_expected[] = { 1.72799012007, 13.4612400765, 4.46126528205, -0.0212528722047, 0.627282377919, 0.302760084926 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1298) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1298) imag");
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
   double A[] = { -0.872, -0.841, 0.108, -0.744, 0.231, -0.513, -0.973, 0.087, 0.348, 0.196, 0.447, 0.307, 0.632, -0.949, 0.322, 0.277, 0.282, 0.831 };
   double X[] = { -0.373, 0.566, 0.92, 0.627, 0.293, -0.434 };
   int incX = -1;
   double x_expected[] = { -1.280839615, 1.560525655, 1.167331, 0.179227, 0.293, -0.434 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1299) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1299) imag");
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
   double A[] = { -0.872, -0.841, 0.108, -0.744, 0.231, -0.513, -0.973, 0.087, 0.348, 0.196, 0.447, 0.307, 0.632, -0.949, 0.322, 0.277, 0.282, 0.831 };
   double X[] = { -0.373, 0.566, 0.92, 0.627, 0.293, -0.434 };
   int incX = -1;
   double x_expected[] = { -0.594503951847, 0.00287302167266, -1.08185265666, -0.859860374254, 0.0331027077244, 1.28233265933 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1300) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1300) imag");
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
   double A[] = { -0.872, -0.841, 0.108, -0.744, 0.231, -0.513, -0.973, 0.087, 0.348, 0.196, 0.447, 0.307, 0.632, -0.949, 0.322, 0.277, 0.282, 0.831 };
   double X[] = { -0.373, 0.566, 0.92, 0.627, 0.293, -0.434 };
   int incX = -1;
   double x_expected[] = { -0.373, 0.566, 1.16074, 0.50314, -0.20669608, 0.37525144 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1301) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1301) imag");
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
   double A[] = { 0.404, 0.667, 0.861, 0.22, 0.298, -0.858, -0.682, -0.969, 0.327, -0.86, 0.125, 0.606, -0.143, -0.865, -0.036, 0.23, -0.776, 0.079 };
   double X[] = { 0.028, -0.804, 0.582, -0.078, -0.126, 0.459 };
   int incX = -1;
   double x_expected[] = { 0.0654496252357, 0.224007771015, -0.752486084395, -0.554870892947, -0.587163401057, 0.166737652215 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1302) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1302) imag");
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
   double A[] = { 0.404, 0.667, 0.861, 0.22, 0.298, -0.858, -0.682, -0.969, 0.327, -0.86, 0.125, 0.606, -0.143, -0.865, -0.036, 0.23, -0.776, 0.079 };
   double X[] = { 0.028, -0.804, 0.582, -0.078, -0.126, 0.459 };
   int incX = -1;
   double x_expected[] = { -0.595558802, -1.147174647, 0.589506, -0.500919, -0.126, 0.459 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1303) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1303) imag");
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
   double A[] = { 0.404, 0.667, 0.861, 0.22, 0.298, -0.858, -0.682, -0.969, 0.327, -0.86, 0.125, 0.606, -0.143, -0.865, -0.036, 0.23, -0.776, 0.079 };
   double X[] = { 0.028, -0.804, 0.582, -0.078, -0.126, 0.459 };
   int incX = -1;
   double x_expected[] = { 3.39346077201, 0.652889512141, -2.33602680355, -2.7859245153, -5.04672104102, -0.334110541026 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1304) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1304) imag");
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
   double A[] = { 0.404, 0.667, 0.861, 0.22, 0.298, -0.858, -0.682, -0.969, 0.327, -0.86, 0.125, 0.606, -0.143, -0.865, -0.036, 0.23, -0.776, 0.079 };
   double X[] = { 0.028, -0.804, 0.582, -0.078, -0.126, 0.459 };
   int incX = -1;
   double x_expected[] = { 0.028, -0.804, -0.109456, -0.217192, -0.41110804, 0.41693792 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1305) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1305) imag");
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
   double A[] = { 0.404, 0.667, 0.861, 0.22, 0.298, -0.858, -0.682, -0.969, 0.327, -0.86, 0.125, 0.606, -0.143, -0.865, -0.036, 0.23, -0.776, 0.079 };
   double X[] = { 0.028, -0.804, 0.582, -0.078, -0.126, 0.459 };
   int incX = -1;
   double x_expected[] = { 7.16970224467, -0.772071373678, 0.833386981173, -0.673826630129, -0.26524050899, 0.465327628365 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1306) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1306) imag");
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
   double A[] = { 0.404, 0.667, 0.861, 0.22, 0.298, -0.858, -0.682, -0.969, 0.327, -0.86, 0.125, 0.606, -0.143, -0.865, -0.036, 0.23, -0.776, 0.079 };
   double X[] = { 0.028, -0.804, 0.582, -0.078, -0.126, 0.459 };
   int incX = -1;
   double x_expected[] = { 0.471459157, -1.566755859, 0.940839, 0.357132, -0.126, 0.459 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1307) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1307) imag");
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
   double A[] = { 0.404, 0.667, 0.861, 0.22, 0.298, -0.858, -0.682, -0.969, 0.327, -0.86, 0.125, 0.606, -0.143, -0.865, -0.036, 0.23, -0.776, 0.079 };
   double X[] = { 0.028, -0.804, 0.582, -0.078, -0.126, 0.459 };
   int incX = -1;
   double x_expected[] = { -0.909961830373, 0.118063054039, -0.0169425582229, -1.00055409731, -1.37205489923, 0.994032418785 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1308) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1308) imag");
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
   double A[] = { 0.404, 0.667, 0.861, 0.22, 0.298, -0.858, -0.682, -0.969, 0.327, -0.86, 0.125, 0.606, -0.143, -0.865, -0.036, 0.23, -0.776, 0.079 };
   double X[] = { 0.028, -0.804, 0.582, -0.078, -0.126, 0.459 };
   int incX = -1;
   double x_expected[] = { 0.028, -0.804, -0.118596, 0.160828, -0.059271004, 0.294435972 };
   cblas_ztbsv(order, uplo, trans, diag, N, K, A, lda, X, incX);
   {
     int i;
     for (i = 0; i < 3; i++) {
       gsl_test_rel(X[2*i], x_expected[2*i], dbleps, "ztbsv(case 1309) real");
       gsl_test_rel(X[2*i+1], x_expected[2*i+1], dbleps, "ztbsv(case 1309) imag");
     };
   };
  };


}
