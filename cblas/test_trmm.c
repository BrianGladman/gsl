#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_trmm () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = 1;
   float A[] = { 0.565, 0.967, -0.969, 0.184 };
   int lda = 2;
   float B[] = { 0.842, -0.918, -0.748, -0.859, -0.463, 0.292 };
   int ldb = 3;
   float B_expected[] = { -0.354923, -0.966391, -0.140256, -0.158056, -0.085192, 0.053728 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1670)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = 1;
   float A[] = { -0.748, 0.548, 0.245, 0.761 };
   int lda = 2;
   float B[] = { 0.349, -0.552, -0.682, -0.71, 0.475, -0.59 };
   int ldb = 3;
   float B_expected[] = { -0.04008, -0.2917, -1.00532, -0.71, 0.475, -0.59 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1671)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = 1;
   float A[] = { 0.788, 0.617, -0.998, -0.97 };
   int lda = 2;
   float B[] = { -0.4, 0.773, 0.074, -0.388, 0.825, -0.608 };
   int ldb = 3;
   float B_expected[] = { -0.3152, 0.609124, 0.058312, 0.77556, -1.5717, 0.515908 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1672)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = 1;
   float A[] = { 0.01, 0.387, -0.953, -0.374 };
   int lda = 2;
   float B[] = { 0.364, 0.09, 0.588, -0.263, 0.584, 0.463 };
   int ldb = 3;
   float B_expected[] = { 0.364, 0.09, 0.588, -0.609892, 0.49823, -0.097364 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1673)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = 1;
   float A[] = { -0.586, -0.426, 0.765, -0.239 };
   int lda = 2;
   float B[] = { -0.673, -0.724, 0.217, -0.672, -0.378, -0.005 };
   int ldb = 2;
   float B_expected[] = { -0.159482, 0.173036, -0.641242, 0.160608, 0.217683, 0.001195 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1674)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = 1;
   float A[] = { -0.668, 0.962, 0.515, 0.292 };
   int lda = 2;
   float B[] = { -0.145, -0.337, 0.718, -0.866, -0.454, -0.439 };
   int ldb = 2;
   float B_expected[] = { -0.318555, -0.337, 0.27201, -0.866, -0.680085, -0.439 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1675)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = 1;
   float A[] = { -0.125, -0.676, 0.181, 0.741 };
   int lda = 2;
   float B[] = { 0.354, -0.366, 0.455, 0.134, -0.564, -0.303 };
   int ldb = 2;
   float B_expected[] = { -0.04425, -0.51051, -0.056875, -0.208286, 0.0705, 0.156741 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1676)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = 1;
   float A[] = { -0.162, 0.542, -0.839, -0.935 };
   int lda = 2;
   float B[] = { 0.216, 0.766, -0.228, -0.097, 0.205, 0.875 };
   int ldb = 2;
   float B_expected[] = { 0.216, 0.883072, -0.228, -0.220576, 0.205, 0.98611 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1677)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = 1;
   float A[] = { -0.353, -0.854, -0.502, 0.591, -0.934, -0.729, 0.063, 0.352, 0.126 };
   int lda = 3;
   float B[] = { 0.2, -0.626, -0.694, -0.889, -0.251, -0.42 };
   int ldb = 3;
   float B_expected[] = { -0.0706, 0.413884, 0.26851, 0.313817, 0.99364, 0.576337 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1678)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = 1;
   float A[] = { -0.864, -0.046, -0.755, 0.12, 0.525, 0.917, 0.571, -0.098, -0.226 };
   int lda = 3;
   float B[] = { -0.905, -0.296, -0.927, -0.813, 0.624, -0.366 };
   int ldb = 3;
   float B_expected[] = { -0.905, -0.25437, -0.515157, -0.813, 0.661398, 0.820023 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1679)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = 1;
   float A[] = { -0.69, -0.927, -0.281, -0.918, -0.527, -0.652, -0.393, -0.954, 0.651 };
   int lda = 3;
   float B[] = { -0.587, 0.788, -0.629, -0.444, 0.515, 0.081 };
   int ldb = 3;
   float B_expected[] = { -0.071157, 0.18479, -0.409479, -0.198243, -0.348679, 0.052731 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1680)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = 1;
   float A[] = { -0.082, -0.077, 0.811, 0.852, 0.224, 0.443, -0.509, 0.171, 0.986 };
   int lda = 3;
   float B[] = { -0.982, 0.388, -0.493, -0.497, -0.605, 0.433 };
   int ldb = 3;
   float B_expected[] = { -0.400487, 0.303697, -0.493, -1.23286, -0.530957, 0.433 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1681)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = 1;
   float A[] = { 0.97, -0.666, 0.066, -0.176, 0.402, 0.286, -0.703, 0.962, 0.912 };
   int lda = 3;
   float B[] = { -0.644, -0.97, 0.814, -0.777, 0.812, 0.254 };
   int ldb = 2;
   float B_expected[] = { -0.62468, -0.9409, 0.440572, -0.141634, 1.97634, 0.166084 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1682)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = 1;
   float A[] = { 0.714, 0.468, 0.859, -0.547, 0.076, 0.542, 0.512, -0.987, -0.167 };
   int lda = 3;
   float B[] = { -0.238, -0.336, 0.402, 0.945, -0.242, -0.062 };
   int ldb = 2;
   float B_expected[] = { -0.238, -0.336, 0.532186, 1.12879, -0.76063, -1.16675 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1683)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = 1;
   float A[] = { -0.723, 0.041, 0.333, -0.682, 0.193, 0.581, 0.963, -0.757, 0.396 };
   int lda = 3;
   float B[] = { 0.047, -0.701, -0.25, -0.779, 0.435, 0.612 };
   int ldb = 2;
   float B_expected[] = { 0.100624, 0.67868, 0.204485, 0.205225, 0.17226, 0.242352 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1684)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = 1;
   float A[] = { -0.13, 0.511, -0.544, 0.938, -0.126, -0.873, 0.118, -0.75, 0.674 };
   int lda = 3;
   float B[] = { -0.927, -0.558, -0.289, -0.66, 0.83, 0.363 };
   int ldb = 2;
   float B_expected[] = { -1.5262, -1.09273, -1.01359, -0.976899, 0.83, 0.363 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1685)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = 0.1;
   float A[] = { -0.625, -0.123, -0.48, -0.088 };
   int lda = 2;
   float B[] = { 0.376, -0.46, -0.813, 0.419, 0.792, 0.226 };
   int ldb = 3;
   float B_expected[] = { -0.0235, 0.02875, 0.0508125, -0.008312, -0.0013116, 0.0080111 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1686)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = 0.1;
   float A[] = { -0.038, -0.105, -0.946, 0.474 };
   int lda = 2;
   float B[] = { -0.757, 0.974, -0.045, -0.809, 0.654, 0.611 };
   int ldb = 3;
   float B_expected[] = { -0.0757, 0.0974, -0.0045, -0.0729515, 0.055173, 0.0615725 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1687)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = 0.1;
   float A[] = { -0.328, 0.713, 0.781, 0.084 };
   int lda = 2;
   float B[] = { -0.097, 0.442, -0.563, 0.065, -0.18, 0.63 };
   int ldb = 3;
   float B_expected[] = { 0.0082581, -0.0285556, 0.0676694, 5.460000e-04, -0.001512, 0.005292 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1688)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = 0.1;
   float A[] = { 0.261, -0.659, -0.536, 0.694 };
   int lda = 2;
   float B[] = { -0.498, 0.692, 0.125, 0.706, -0.118, -0.907 };
   int ldb = 3;
   float B_expected[] = { -0.0876416, 0.0755248, 0.0611152, 0.0706, -0.0118, -0.0907 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1689)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = 0.1;
   float A[] = { -0.669, 0.416, 0.761, -0.359 };
   int lda = 2;
   float B[] = { -0.305, -0.675, -0.442, 0.566, 0.064, 0.962 };
   int ldb = 2;
   float B_expected[] = { 0.0204045, 0.001022, 0.0295698, -0.0539556, -0.0042816, -0.0296654 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1690)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = 0.1;
   float A[] = { 0.565, 0.386, 0.643, -0.028 };
   int lda = 2;
   float B[] = { 0.863, -0.241, 0.766, 0.656, -0.977, 0.274 };
   int ldb = 2;
   float B_expected[] = { 0.0863, 0.0313909, 0.0766, 0.114854, -0.0977, -0.0354211 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1691)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = 0.1;
   float A[] = { 0.116, 0.534, 0.043, 0.73 };
   int lda = 2;
   float B[] = { -0.758, -0.63, -0.043, 0.666, -0.088, 0.382 };
   int ldb = 2;
   float B_expected[] = { -0.0424348, -0.04599, 0.0350656, 0.048618, 0.019378, 0.027886 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1692)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = 0.1;
   float A[] = { 0.48, -0.63, -0.786, -0.437 };
   int lda = 2;
   float B[] = { 0.945, 0.528, -0.855, -0.587, 0.062, 0.372 };
   int ldb = 2;
   float B_expected[] = { 0.061236, 0.0528, -0.048519, -0.0587, -0.017236, 0.0372 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1693)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = 0.1;
   float A[] = { -0.822, -0.068, 0.119, -0.244, -0.05, 0.685, 0.752, -0.059, -0.935 };
   int lda = 3;
   float B[] = { -0.431, -0.753, -0.319, 0.164, 0.979, 0.885 };
   int ldb = 3;
   float B_expected[] = { 0.0367525, -0.0180865, 0.0298265, -0.0096065, 0.0557275, -0.0827475 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1694)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = 0.1;
   float A[] = { 0.97, -0.408, 0.174, -0.308, 0.997, -0.484, 0.322, -0.183, 0.849 };
   int lda = 3;
   float B[] = { -0.571, 0.696, -0.256, -0.178, 0.098, 0.004 };
   int ldb = 3;
   float B_expected[] = { -0.0899512, 0.0819904, -0.0256, -0.0217288, 0.0096064, 4.000000e-04 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1695)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = 0.1;
   float A[] = { -0.831, 0.73, 0.407, 0.721, 0.086, -0.294, 0.941, -0.656, -0.066 };
   int lda = 3;
   float B[] = { -0.051, -0.343, -0.98, 0.722, -0.372, 0.466 };
   int ldb = 3;
   float B_expected[] = { 0.0042381, -0.0066269, 0.0241697, -0.0599982, 0.048857, 0.0892678 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1696)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = 0.1;
   float A[] = { 0.472, 0.137, -0.341, 0.386, -0.578, 0.863, -0.415, -0.547, -0.023 };
   int lda = 3;
   float B[] = { 0.582, 0.141, -0.306, -0.047, -0.162, -0.784 };
   int ldb = 3;
   float B_expected[] = { 0.0582, 0.0365652, -0.0624657, -0.0047, -0.0180142, -0.0675881 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1697)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = 0.1;
   float A[] = { -0.775, 0.762, -0.038, -0.8, 0.626, -0.701, 0.639, 0.239, 0.34 };
   int lda = 3;
   float B[] = { 0.42, 0.917, 0.485, 0.844, -0.832, 0.179 };
   int ldb = 2;
   float B_expected[] = { -0.124515, -0.127149, 0.0104762, 0.0571125, -0.028288, 0.006086 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1698)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = 0.1;
   float A[] = { -0.675, 0.283, 0.785, -0, -0.592, -0.661, 0.149, -0.129, 0.149 };
   int lda = 3;
   float B[] = { 0.964, -0.575, -0.215, 0.953, 0.527, -0.418 };
   int ldb = 2;
   float B_expected[] = { 0.104252, -0.0637282, -0.0282983, 0.100692, 0.0527, -0.0418 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1699)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = 0.1;
   float A[] = { -0.225, -0.943, 0.839, 0.759, 0.752, 0.807, 0.288, -0.276, 0.434 };
   int lda = 3;
   float B[] = { -0.234, 0.275, 0.658, -0.423, -0.807, -0.683 };
   int ldb = 2;
   float B_expected[] = { 0.005265, -0.0061875, 0.0715478, -0.0577421, -0.0015558, -0.0407058 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1700)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = 0.1;
   float A[] = { -0.043, -0.983, 0.479, -0.136, 0.048, 0.745, -0.408, -0.731, -0.953 };
   int lda = 3;
   float B[] = { 0.917, 0.682, -0.32, 0.557, -0.302, 0.989 };
   int ldb = 2;
   float B_expected[] = { 0.0917, 0.0682, -0.122141, -0.0113406, -0.0101157, 0.173064 };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1701)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { -0.561, -0.114, -0.148, 0.488 };
   int lda = 2;
   double B[] = { 0.684, 0.38, 0.419, -0.361, 0.378, -0.423 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1702)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { -0.378, 0.607, 0.41, 0.418 };
   int lda = 2;
   double B[] = { 0.146, -0.688, -0.953, -0.983, 0.237, 0.128 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1703)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { -0.31, 0.277, -0.587, 0.885 };
   int lda = 2;
   double B[] = { -0.221, -0.831, -0.319, -0.547, -0.577, 0.295 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1704)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { -0.577, 0.861, -0.439, -0.916 };
   int lda = 2;
   double B[] = { -0.933, -0.582, 0.528, 0.268, -0.804, 0.62 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1705)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { -0.824, -0.119, -0.399, -0.653 };
   int lda = 2;
   double B[] = { 0.452, -0.168, 0.256, 0.554, 0.342, 0.318 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1706)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { -0.299, 0.837, -0.03, 0.552 };
   int lda = 2;
   double B[] = { -0.83, -0.82, -0.362, -0.252, -0.062, -0.942 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1707)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { -0.545, -0.107, 0.096, 0.183 };
   int lda = 2;
   double B[] = { -0.43, 0.841, 0.035, 0.7, 0.637, 0.095 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1708)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { 0.626, 0.123, -0.959, 0.971 };
   int lda = 2;
   double B[] = { 0.185, -0.218, -0.074, 0.49, 0.802, -0.454 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1709)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { -0.131, 0.048, 0.148, 0.834, -0.98, -0.009, -0.727, 0.241, 0.276 };
   int lda = 3;
   double B[] = { 0.75, -0.664, -0.136, -0.793, -0.742, 0.126 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1710)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { 0.431, -0.387, 0.427, 0.495, 0.282, 0.158, -0.335, 0.535, -0.978 };
   int lda = 3;
   double B[] = { 0.518, -0.489, 0.899, -0.375, 0.376, -0.831 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1711)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { -0.669, -0.976, -0.2, 0.661, -0.975, -0.965, -0.861, -0.779, -0.73 };
   int lda = 3;
   double B[] = { 0.31, 0.023, -0.853, 0.632, -0.174, 0.608 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1712)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { 0.153, -0.408, -0.127, -0.634, -0.384, -0.815, 0.051, -0.096, 0.476 };
   int lda = 3;
   double B[] = { 0.343, -0.665, -0.348, 0.748, 0.893, 0.91 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1713)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { -0.918, -0.19, 0.829, 0.942, 0.885, 0.087, 0.321, 0.67, -0.475 };
   int lda = 3;
   double B[] = { 0.377, 0.931, 0.291, -0.603, -0.617, 0.402 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1714)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { -0.598, -0.232, -0.64, 0.595, 0.642, -0.921, -0.679, -0.846, -0.921 };
   int lda = 3;
   double B[] = { 0.032, -0.036, -0.278, -0.83, 0.922, -0.701 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1715)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { 0.341, -0.858, -0.559, 0.499, -0.114, 0.57, 0.847, -0.612, 0.593 };
   int lda = 3;
   double B[] = { 0.672, 0.292, 0.752, 0.842, 0.625, 0.967 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1716)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { 0.958, 0.823, -0.181, 0.141, 0.932, 0.097, -0.636, 0.844, 0.205 };
   int lda = 3;
   double B[] = { 0.113, -0.658, 0.703, -0.023, -0.384, 0.439 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1717)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { 0.675, -0.468, -0.564, 0.71 };
   int lda = 2;
   double B[] = { -0.401, -0.823, 0.342, -0.384, 0.344, 0.18 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1718)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { 0.932, -0.388, 0.432, -0.167 };
   int lda = 2;
   double B[] = { -0.624, 0.023, 0.065, 0.678, 0.044, -0.472 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1719)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { -0.738, 0.649, -0.171, -0.462 };
   int lda = 2;
   double B[] = { -0.277, -0.519, -0.501, -0.024, -0.767, -0.591 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1720)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { -0.17, -0.184, -0.243, 0.907 };
   int lda = 2;
   double B[] = { 0.593, 0.131, -0.317, -0.254, -0.948, 0.002 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1721)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { 0.06, -0.838, -0.455, -0.715 };
   int lda = 2;
   double B[] = { -0.423, 0.665, -0.023, -0.872, -0.313, -0.698 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1722)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { -0.506, 0.792, 0.338, -0.155 };
   int lda = 2;
   double B[] = { -0.257, -0.19, 0.201, 0.685, 0.663, 0.302 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1723)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { 0.739, -0.996, 0.182, 0.626 };
   int lda = 2;
   double B[] = { 0.009, 0.485, -0.633, -0.08, -0.579, 0.223 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1724)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { 0.777, 0.723, 0.378, 0.98 };
   int lda = 2;
   double B[] = { 0.291, -0.267, -0.076, 0.103, -0.021, -0.866 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1725)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { -0.771, 0.469, 0.822, -0.619, 0.953, -0.706, 0.318, 0.559, -0.68 };
   int lda = 3;
   double B[] = { -0.32, 0.362, 0.719, -0.661, -0.504, 0.595 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1726)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { 0.073, -0.501, -0.561, -0.229, -0.533, -0.138, 0.924, -0.164, -0.023 };
   int lda = 3;
   double B[] = { -0.208, 0.49, 0.827, 0.641, -0.884, -0.624 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1727)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { 0.33, -0.649, -0.43, -0.266, 0.787, 0.449, 0.435, -0.774, -0.447 };
   int lda = 3;
   double B[] = { -0.687, -0.459, 0.189, 0.762, -0.039, 0.047 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1728)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { 0.981, 0.242, 0.581, 0.064, 0.792, -0.529, 0.461, 0.224, -0.419 };
   int lda = 3;
   double B[] = { 0.285, 0.274, -0.912, 0.601, 0.24, 0.06 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1729)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { -0.582, 0.269, -0.587, 0.68, -0.59, -0.936, 0.236, -0.728, -0.434 };
   int lda = 3;
   double B[] = { 0.113, 0.468, 0.943, 0.48, 0.215, -0.525 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1730)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { -0.344, -0.938, 0.556, -0.678, -0.612, -0.519, -0.578, -0.848, 0.699 };
   int lda = 3;
   double B[] = { 0.915, -0.118, 0.538, -0.186, -0.413, -0.216 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1731)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { -0.843, 0.54, -0.892, -0.296, 0.786, 0.136, 0.731, -0.418, -0.118 };
   int lda = 3;
   double B[] = { -0.775, 0.5, -0.399, -0.709, 0.779, 0.774 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1732)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0;
   double A[] = { -0.765, 0.233, 0.318, 0.547, -0.469, 0.023, -0.867, 0.687, -0.912 };
   int lda = 3;
   double B[] = { 0.019, -0.145, 0.472, 0.333, 0.527, -0.224 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1733)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.852, -0.409, 0.871, -0.854, -0.493, 0.444, 0.973, 0.027 };
   int lda = 2;
   float B[] = { -0.561, 0.132, 0.689, 0.653, -0.758, -0.109, -0.596, 0.395, -0.561, 0.378, 0.21, 0.51 };
   int ldb = 3;
   float B_expected[] = { -0.0970014, 0.0350174, 0.0029825, -0.048577, -0.066776, 0.121969, -0.0368243, -0.0590573, -0.0352647, -0.0556059, -0.05019, 0.019056 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1734) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1734) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { 0.349, 0, -0.462, 0.91, -0.693, 0.587, -0.617, 0.112 };
   int lda = 2;
   float B[] = { 0.842, -0.473, 0.825, 0.866, 0.986, 0.686, 0.346, 0.299, -0.659, 0.009, 0.007, -0.478 };
   int ldb = 3;
   float B_expected[] = { 0.0296278, 0.0410058, -0.0262152, 0.112127, -0.0913206, 0.141775, -0.0299, 0.0346, -0.0009, -0.0659, 0.0478, 0.0007 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1735) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1735) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { 0.661, -0.823, 0.28, 0.171, 0.267, 0.66, 0.844, 0.472 };
   int lda = 2;
   float B[] = { -0.256, -0.518, -0.933, 0.066, -0.513, -0.286, 0.109, 0.372, -0.183, 0.482, 0.362, -0.436 };
   int ldb = 3;
   float B_expected[] = { 0.013171, -0.059553, -0.0811485, -0.0562395, -0.0233153, -0.0574471, -0.005815, 0.018994, 0.0277726, -0.0674627, 0.0612062, 0.0563109 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1736) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1736) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { 0.623, 0.314, -0.594, 0.717, 0.566, 0.001, -0.411, -0.387 };
   int lda = 2;
   float B[] = { -0.083, 0.937, -0.814, 0.9, -0.042, 0.678, -0.928, 0.228, 0.965, -0.16, 0.006, -0.281 };
   int ldb = 3;
   float B_expected[] = { -0.0937, -0.0083, -0.09, -0.0814, -0.0678, -0.0042, -0.0758259, -0.0975915, -0.0348586, 0.0503376, -0.0102706, -0.001845 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1737) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1737) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.247, -0.582, 0.651, -0.534, -0.491, 0.346, 0.936, -0.227 };
   int lda = 2;
   float B[] = { -0.002, -0.02, 0.162, -0.62, 0.632, -0.07, 0.352, 0.042, 0.574, 0.272, -0.139, 0.012 };
   int ldb = 2;
   float B_expected[] = { -0.0366576, 0.0123832, 0.0617094, 0.0010892, 0.0249364, -0.0384208, 0.0040592, 0.0339006, 0.0455238, 0.0080623, -0.0042785, -0.012738 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1738) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1738) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.152, 0.395, -0.077, -0.191, -0.757, 0.858, -0.494, -0.734 };
   int lda = 2;
   float B[] = { -0.166, -0.413, -0.373, 0.915, -0.824, -0.066, -0.114, -0.921, 0.862, 0.312, 0.221, 0.699 };
   int ldb = 2;
   float B_expected[] = { 0.142569, -0.0668709, -0.0915, -0.0373, -0.0533385, 0.0052516, 0.0921, -0.0114, 0.0027525, 0.0094961, -0.0699, 0.0221 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1739) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1739) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.426, 0.817, -0.993, -0.882, 0.615, 0.627, -0.238, -0.903 };
   int lda = 2;
   float B[] = { 0.895, 0.849, 0.811, 0.402, 0.074, -0.493, -0.548, -0.82, 0.323, 0.301, 0.612, -0.092 };
   int ldb = 2;
   float B_expected[] = { -0.0369541, -0.10749, 0.246046, 0.0030071, -0.0270476, 0.0371257, -0.111428, -0.111834, -0.0135665, -0.0383515, 0.111452, -0.0283989 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1740) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1740) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.451, -0.754, -0.673, 0.433, -0.712, -0.033, -0.588, 0.116 };
   int lda = 2;
   float B[] = { 0.787, -0.377, -0.854, -0.464, 0.118, 0.231, 0.362, -0.457, -0.076, 0.373, -0.286, -0.468 };
   int ldb = 2;
   float B_expected[] = { 0.0377, 0.0787, -0.0130492, -0.122041, -0.0231, 0.0118, 0.0561369, 0.0182563, -0.0373, -0.0076, 0.0751937, -0.0396361 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1741) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1741) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.454, 0.494, 0.424, -0.907, 0.339, -0.141, 0.169, 0.364, -0.607, 0.955, -0.156, 0.962, -0.254, 0.079, 0.209, 0.946, 0.93, 0.677 };
   int lda = 3;
   float B[] = { -0.99, -0.484, 0.915, -0.383, 0.228, 0.797, 0.597, 0.765, -0.629, 0.002, -0.89, 0.077 };
   int ldb = 3;
   float B_expected[] = { 0.0269324, 0.0688556, -0.179902, -0.104839, -0.181106, -0.0505677, 0.0052392, -0.0648948, 0.0819028, 0.132688, 0.0961172, -0.0473381 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1742) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1742) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.008, -0.654, 0.174, 0.448, 0.388, -0.108, -0.479, -0.708, -0.035, 0.816, 0.487, 0.22, -0.482, 0.57, -0.317, 0.203, -0.547, -0.415 };
   int lda = 3;
   float B[] = { 0.651, 0.187, 0.591, -0.007, 0.171, -0.923, -0.029, -0.685, -0.049, 0.135, 0.578, 0.979 };
   int ldb = 3;
   float B_expected[] = { -0.0187, 0.0651, -0.0317186, 0.0620498, 0.0794141, 0.0733141, 0.0685, -0.0029, -0.0002818, 0.0252834, -0.0771317, 0.0439205 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1743) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1743) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { 0.952, 0.29, 0.944, 0.294, -0.762, -0.7, -0.949, 0.167, 0.307, 0.904, -0.428, -0.411, 0.496, 0.004, -0.611, -0.09, -0.846, 0.081 };
   int lda = 3;
   float B[] = { 0.782, -0.035, -0.441, -0.791, -0.09, -0.56, -0.438, -0.691, 0.88, 0.545, -0.55, 0.595 };
   int ldb = 3;
   float B_expected[] = { -0.0592352, 0.126282, 0.0291241, 0.0584267, -0.046647, 0.01215, 0.0862177, -0.14179, -0.064879, 0.016708, 0.054792, 0.0417105 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1744) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1744) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.519, 0.708, -0.934, -0.219, 0.376, -0.967, 0.322, -0.355, 0.972, -0.156, -0.735, 0.928, 0.084, -0.267, -0.152, 0.434, 0.267, 0.983 };
   int lda = 3;
   float B[] = { -0.54, 0.149, 0.574, 0.742, 0.704, 0.459, -0.9, 0.04, 0.538, -0.858, 0.467, 0.686 };
   int ldb = 3;
   float B_expected[] = { -0.0034742, 0.0089927, -0.0977768, 0.0267786, -0.0459, 0.0704, 0.0494331, -0.0808964, 0.0759594, 0.0169292, -0.0686, 0.0467 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1745) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1745) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.541, 0.67, 0.014, 0.446, 0.086, -0.525, 0.033, -0.932, 0.977, 0.321, -0.651, 0.027, 0.409, 0.328, 0.359, -0.615, 0.419, -0.25 };
   int lda = 3;
   float B[] = { -0.156, 0.666, -0.231, 0.691, 0.935, -0.481, -0.142, -0.117, 0.529, 0.526, 0.266, 0.417 };
   int ldb = 2;
   float B_expected[] = { 0.0464826, -0.0361824, 0.0528601, -0.0337999, 0.0002432, 0.168346, -0.0078204, 0.0535212, 0.0438334, 0.0110749, -0.0360401, -0.0228356 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1746) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1746) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.459, -0.349, -0.335, 0.008, 0.866, 0.978, -0.869, -0.361, -0.711, 0.712, 0.207, 0.305, 0.766, -0.262, 0.012, -0.333, 0.617, 0.91 };
   int lda = 3;
   float B[] = { -0.138, -0.256, -0.319, -0.771, 0.674, -0.565, -0.779, -0.516, -0.017, -0.097, -0.555, 0.308 };
   int ldb = 2;
   float B_expected[] = { 0.0256, -0.0138, 0.0771, -0.0319, 0.0292718, 0.0701506, -0.0269158, -0.078012, 0.0488162, -0.0369837, -0.0054207, -0.118253 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1747) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1747) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.825, -0.785, -0.605, -0.508, 0.763, -0.578, -0.167, -0.233, 0.011, -0.853, 0.24, 0.192, 0.293, -0.72, -0.348, 0.023, -0.145, -0.493 };
   int lda = 3;
   float B[] = { 0.305, -0.255, 0.882, 0.883, 0.088, -0.473, 0.135, -0.063, -0.671, 0.473, 0.874, 0.548 };
   int ldb = 2;
   float B_expected[] = { -0.0961148, -0.0983903, 0.153836, 0.0835432, 0.0095579, -0.0654357, -0.018348, 0.005229, -0.0262218, 0.0330484, 0.0510342, 0.0143434 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1748) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1748) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { 0.63, 0.353, 0.445, 0.845, 0.273, -0.135, 0.03, 0.936, 0.141, 0.638, -0.399, 0.343, -0.037, -0.335, -0.089, 0.081, 0.987, -0.256 };
   int lda = 3;
   float B[] = { -0.567, 0.803, 0.168, 0.744, -0.328, 0.835, -0.852, 0.702, 0.21, -0.618, 0.666, -0.303 };
   int ldb = 2;
   float B_expected[] = { -0.0700351, -0.144464, -0.0163821, -0.0663417, -0.115361, -0.0199816, -0.105134, -0.10138, 0.0618, 0.021, 0.0303, 0.0666 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1749) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1749) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.741, 0.904, -0.599, 0.753, -0.297, 0.38, -0.056, -0.715 };
   int lda = 2;
   float B[] = { 0.646, -0.447, -0.147, 0.314, -0.713, 0.187, -0.589, 0.287, -0.809, -0.293, 0.418, 0.778 };
   int ldb = 3;
   float B_expected[] = { -0.0915211, -0.0074598, 0.0365562, -0.0174929, 0.0783119, 0.0359285, -0.115925, 0.0187826, -0.0296066, -0.031258, 0.099134, 0.0819138 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1750) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1750) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.645, 0.756, 0.709, -0.657, -0.023, -0.714, 0.03, 0.239 };
   int lda = 2;
   float B[] = { -0.16, 0.254, -0.68, 0.183, -0.402, -0.259, 0.104, -0.09, 0.944, 0.729, -0.378, -0.792 };
   int ldb = 3;
   float B_expected[] = { -0.0254, -0.016, -0.0183, -0.068, 0.0259, -0.0402, -0.0195206, 0.0157438, -0.130551, 0.0582111, 0.0711517, -0.0833181 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1751) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1751) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { 0.25, -0.038, 0.377, -0.209, 0.166, -0.073, -0.24, 0.938 };
   int lda = 2;
   float B[] = { 0.26, 0.696, -0.183, 0.668, -0.08, -0.938, -0.837, -0.509, 0.781, -0.063, -0.953, 0.227 };
   int ldb = 3;
   float B_expected[] = { -0.0140727, -0.0084651, -0.0106483, 0.0104681, 0.0124209, -0.0197271, 0.0662946, 0.0678322, -0.0747698, -0.0128346, 0.0948394, 0.0015794 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1752) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1752) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.668, 0.804, 0.608, -0.682, -0.513, 0.521, 0.878, -0.664 };
   int lda = 2;
   float B[] = { -0.871, 0.699, 0.561, 0.823, -0.787, 0.055, -0.686, 0.361, -0.662, -0.192, -0.301, -0.167 };
   int ldb = 3;
   float B_expected[] = { -0.0156401, -0.0707163, -0.0576594, 0.100064, 0.001615, -0.054558, -0.0361, -0.0686, 0.0192, -0.0662, 0.0167, -0.0301 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1753) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1753) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { 0.091, 0.189, -0.221, 0.749, 0.354, -0.397, 0.105, -0.944 };
   int lda = 2;
   float B[] = { 0.731, -0.446, 0.983, 0.793, 0.533, 0.386, -0.781, -0.063, 0.875, -0.128, -0.179, -0.079 };
   int ldb = 2;
   float B_expected[] = { -0.0097573, 0.0150815, 0.129278, 0.0933519, -0.0135863, -0.0024451, -0.0655692, 0.0200447, -0.0153727, 0.0103817, 0.0232006, 0.0165563 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1754) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1754) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { 0.676, 0.644, 0.03, 0.456, 0.002, -0.909, 0.984, 0.771 };
   int lda = 2;
   float B[] = { 0.65, 0.005, -0.883, -0.154, -0.137, -0.137, 0.531, -0.49, 0.052, 0.273, -0.602, 0.655 };
   int ldb = 2;
   float B_expected[] = { -0.0005, 0.065, 0.074484, -0.0877155, 0.0137, -0.0137, 0.0365741, 0.0406193, -0.0273, 0.0052, -0.0608278, -0.0353739 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1755) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1755) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { 0.832, -0.559, 0.188, -0.488, -0.051, -0.057, 0.909, 0.006 };
   int lda = 2;
   float B[] = { -0.408, 0.303, 0.03, 0.529, -0.584, -0.976, 0.443, -0.762, 0.43, 0.812, -0.075, 0.06 };
   int ldb = 2;
   float B_expected[] = { -0.056498, 0.0093713, -0.0481041, 0.0024096, 0.0845016, -0.132004, 0.069, 0.0407259, -0.0483094, 0.0826848, -0.005409, -0.0068535 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1756) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1756) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.15, -0.297, 0.821, -0.576, -0.572, 0.924, 0.106, -0.131 };
   int lda = 2;
   float B[] = { -0.271, 0.793, -0.232, -0.967, -0.466, 0.37, -0.745, -0.156, -0.091, -0.877, 0.595, 0.448 };
   int ldb = 2;
   float B_expected[] = { -0.0132725, -0.101846, 0.0967, -0.0232, -0.0671044, -0.11675, 0.0156, -0.0745, 0.0851912, 0.0655543, -0.0448, 0.0595 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1757) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1757) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { 0.002, 0, 0.626, -0.148, 0.874, 0.229, -0.227, -0.55, -0.895, 0.586, 0.934, 0.618, 0.958, -0.543, 0.49, 0.671, -0.871, 0.227 };
   int lda = 3;
   float B[] = { -0.415, 0.156, -0.539, -0.247, -0.725, 0.932, 0.565, 0.454, -0.118, 0.693, -0.968, -0.601 };
   int ldb = 3;
   float B_expected[] = { -0.0574005, -0.122188, -0.0327649, -0.0625979, 0.0976347, 0.0419911, 0.0294756, -0.0678577, 0.184894, -0.0833182, -0.0303735, 0.0979555 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1758) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1758) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { 0.89, 0.309, -0.786, 0.999, 0.511, 0.599, 0.385, -0.615, 0.527, -0.328, -0.078, -0.666, 0.004, -0.69, -0.281, -0.438, 0.456, 0.524 };
   int lda = 3;
   float B[] = { -0.648, -0.189, -0.295, 0.477, 0.509, 0.685, 0.875, 0.277, -0.34, -0.632, -0.453, -0.798 };
   int ldb = 3;
   float B_expected[] = { 0.0203701, -0.104287, -0.0084576, 0.0121508, -0.0685, 0.0509, 0.0245033, 0.202013, 0.0268058, -0.0836134, 0.0798, -0.0453 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1759) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1759) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { 0.772, 0.686, 0.693, 0.803, -0.328, -0.627, -0.869, -0.656, -0.055, -0.366, -0.981, -0.151, 0.147, -0.368, -0.824, -0.454, -0.445, -0.794 };
   int lda = 3;
   float B[] = { -0.268, -0.521, -0.685, -0.618, 0.508, 0.525, -0.492, -0.502, -0.997, 0.28, 0.63, 0.664 };
   int ldb = 3;
   float B_expected[] = { 0.058606, 0.015051, -0.0913257, -0.0297397, -0.0205282, 0.0243534, 0.0725056, -0.0035452, -0.110849, 0.0255551, 0.046652, 0.0938454 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1760) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1760) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.317, -0.822, 0.732, 0.383, 0.457, 0.443, 0.529, -0.949, -0.927, -0.65, -0.471, -0.624, -0.731, 0.107, -0.142, 0.623, 0.159, -0.419 };
   int lda = 3;
   float B[] = { 0.292, -0.665, -0.93, 0.517, 0.123, -0.181, 0.325, 0.954, -0.988, -0.128, 0.637, -0.997 };
   int ldb = 3;
   float B_expected[] = { 0.0665, 0.0292, 0.0111893, -0.140662, 0.0316445, -0.0209328, -0.0954, 0.0325, -0.0068241, 0.0089271, 0.225695, 0.0517387 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1761) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1761) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.809, 0.393, -0.015, -0.273, -0.956, 0.49, 0.365, -0.386, 0.941, 0.992, 0.297, 0.761, 0.425, -0.605, 0.672, 0.725, -0.077, -0.628 };
   int lda = 3;
   float B[] = { 0.21, 0.153, 0.218, -0.129, 0.736, -0.006, 0.502, -0.165, 0.242, 0.915, 0.67, 0.07 };
   int ldb = 2;
   float B_expected[] = { 0.0085068, 0.069273, 0.0439562, 0.0320975, -0.15148, 0.0197777, -0.0875509, 0.103555, 0.0222431, 0.0555986, 0.042615, -0.000763 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1762) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1762) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.187, -0.508, -0.987, -0.861, 0.519, 0.752, -0.117, 0.972, 0.068, -0.752, 0.344, 0.074, -0.343, 0, -0.876, 0.857, -0.148, -0.933 };
   int lda = 3;
   float B[] = { 0.827, 0.958, 0.395, 0.878, 0.88, -0.896, -0.771, -0.355, -0.979, 0.329, -0.166, -0.644 };
   int ldb = 2;
   float B_expected[] = { -0.180535, 0.193075, -0.0391015, 0.0887205, 0.202321, 0.145565, -0.0066882, -0.0073676, -0.0329, -0.0979, 0.0644, -0.0166 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1763) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1763) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { -0.622, 0.022, -0.966, 0.704, 0.43, -0.451, -0.221, 0.969, 0.977, 0.021, -0.725, -0.382, 0.779, 0.957, 0.25, 0.832, 0.029, -0.903 };
   int lda = 3;
   float B[] = { 0.315, -0.297, -0.864, 0.519, -0.601, -0.119, 0.028, 0.072, -0.171, 0.648, 0.159, -0.623 };
   int ldb = 2;
   float B_expected[] = { -0.0191664, -0.0189396, 0.0341826, 0.052599, -0.0379778, -0.067988, 0.103868, 0.0495092, -0.0219287, 0.0971955, -0.0388294, -0.0688205 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1764) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1764) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0.1};
   float A[] = { 0.106, 0.87, 0.21, 0.463, -0.496, -0.981, -0.354, -0.604, -0.149, -0.384, -0.958, -0.502, -0.579, 0.736, -0.322, 0.028, 0.193, 0.14 };
   int lda = 3;
   float B[] = { -0.812, 0.518, 0.085, -0.447, -0.443, 0.928, -0.972, 0.889, 0.605, -0.258, -0.025, 0.98 };
   int ldb = 2;
   float B_expected[] = { -0.0518, -0.0812, 0.0447, 0.0085, -0.0660824, -0.0853354, -0.0834485, -0.0747189, 0.0384994, 0.240616, -0.0754609, 0.0871787 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1765) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1765) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0};
   float A[] = { -0.553, 0.204, -0.793, -0.558, 0.741, 0.26, 0.945, -0.757 };
   int lda = 2;
   float B[] = { -0.515, 0.532, -0.321, 0.326, -0.81, -0.924, 0.474, 0.985, -0.03, 0.406, 0.923, -0.956 };
   int ldb = 3;
   float B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1766) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1766) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0};
   float A[] = { -0.41, -0.804, 0.988, -0.715, -0.281, -0.89, 0.389, -0.408 };
   int lda = 2;
   float B[] = { 0.917, 0.541, -0.108, -0.965, 0.524, 0.04, -0.736, -0.643, -0.202, 0.86, 0.346, -0.017 };
   int ldb = 3;
   float B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1767) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1767) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0};
   float A[] = { 0.153, -0.812, -0.742, -0.18, 0.473, 0.023, -0.433, 0.559 };
   int lda = 2;
   float B[] = { 0.078, -0.691, -0.717, -0.637, -0.016, 0.375, -0.902, -0.343, 0.155, 0.563, 0.419, 0.451 };
   int ldb = 3;
   float B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1768) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1768) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0};
   float A[] = { -0.288, 0.241, 0.593, -0.597, -0.469, 0.735, 0.193, -0.104 };
   int lda = 2;
   float B[] = { -0.835, 0.037, -0.762, 0.782, -0.874, -0.867, -0.81, -0.577, 0.352, 0.827, 0.237, -0.861 };
   int ldb = 3;
   float B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1769) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1769) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0};
   float A[] = { 0.441, -0.217, 0.679, 0.106, -0.76, -0.258, -0.956, -0.858 };
   int lda = 2;
   float B[] = { -0.802, 0.163, 0.293, 0.54, 0.228, 0.071, 0.942, 0.345, 0.591, 0.654, 0.382, -0.892 };
   int ldb = 2;
   float B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1770) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1770) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0};
   float A[] = { 0.916, 0.909, 0.834, 0.38, 0.391, -0.412, -0.714, -0.456 };
   int lda = 2;
   float B[] = { -0.151, 0.818, 0.717, -0.812, -0.649, -0.107, -0.454, 0.785, 0.86, 0.992, -0.244, -0.242 };
   int ldb = 2;
   float B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1771) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1771) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0};
   float A[] = { -0.992, 0.284, -0.01, 0.182, 0.527, -0.348, -0.509, 0.839 };
   int lda = 2;
   float B[] = { 0.504, -0.782, -0.88, 0.079, 0.216, 0.525, 0.198, 0.851, -0.102, -0.046, 0.079, -0.045 };
   int ldb = 2;
   float B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1772) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1772) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0};
   float A[] = { 0.985, 0.068, -0.095, -0.575, -0.607, 0.893, 0.085, 0.145 };
   int lda = 2;
   float B[] = { -0.149, 0.592, 0.588, -0.62, -0.409, -0.344, 0.263, 0.759, -0.026, -0.609, 0.507, -0.084 };
   int ldb = 2;
   float B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1773) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1773) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0};
   float A[] = { 0.36, 0.508, -0.771, -0.442, -0.671, -0.691, -0.771, 0.113, 0.282, 0.312, 0.564, -0.568, -0.743, 0.912, -0.395, 0.503, -0.167, -0.581 };
   int lda = 3;
   float B[] = { -0.018, 0.574, -0.144, -0.758, 0.53, 0.623, -0.771, -0.733, 0.932, -0.192, 0.997, 0.773 };
   int ldb = 3;
   float B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1774) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1774) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0};
   float A[] = { -0.627, 0.511, -0.246, -0.091, 0.66, -0.983, 0.99, 0.057, -0.259, 0.18, 0.606, 0.058, -0.238, 0.717, 0.358, -0.851, -0.71, -0.683 };
   int lda = 3;
   float B[] = { -0.907, 0.956, 0.56, -0.057, 0.054, -0.77, 0.868, -0.843, 0.645, -0.554, -0.958, 0.988 };
   int ldb = 3;
   float B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1775) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1775) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0};
   float A[] = { -0.882, 0.431, -0.868, -0.098, -0.006, -0.639, 0.757, -0.009, -0.821, 0.45, 0.347, 0.801, 0.314, 0.936, -0.725, 0.956, 0.536, 0.771 };
   int lda = 3;
   float B[] = { 0.38, -0.435, 0.977, 0.296, -0.624, -0.53, 0.73, -0.837, 0.105, 0.189, 0.362, -0.664 };
   int ldb = 3;
   float B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1776) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1776) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0};
   float A[] = { 0.595, -0.775, 0.75, 0.16, -0.572, 0.658, 0.216, 0.557, -0.279, 0.095, -0.495, 0.503, 0.071, -0.03, -0.116, 0.78, -0.104, 0.073 };
   int lda = 3;
   float B[] = { 0.948, 0.749, -0.854, 0.972, 0.704, 0.187, 0.347, 0.303, -0.865, 0.123, -0.041, 0.152 };
   int ldb = 3;
   float B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1777) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1777) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0};
   float A[] = { -0.617, -0.331, -0.074, 0.719, -0.469, -0.852, 0.25, -0.175, -0.719, -0.613, -0.321, 0.973, -0.337, -0.35, 0.607, -0.553, 0.688, 0.463 };
   int lda = 3;
   float B[] = { 0.568, -0.471, -0.947, -0.205, 0.835, -0.859, 0.27, -0.599, 0.171, -0.514, 0.939, 0.176 };
   int ldb = 2;
   float B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1778) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1778) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0};
   float A[] = { 0.99, -0.857, 0.728, -0.31, -0.506, -0.393, 0.97, 0.282, 0.375, -0.286, -0.496, -0.057, 0.186, -0.34, 0.608, -0.52, 0.921, -0.875 };
   int lda = 3;
   float B[] = { -0.929, 0.885, 0.864, -0.548, 0.393, 0.391, 0.033, 0.186, 0.949, -0.435, 0.986, -0.995 };
   int ldb = 2;
   float B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1779) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1779) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0};
   float A[] = { -0.101, -0.92, 0.969, -0.017, -0.016, -0.024, -0.11, 0.219, -0.287, -0.937, 0.619, 0.166, -0.068, 0.753, 0.374, 0.076, 0.79, -0.64 };
   int lda = 3;
   float B[] = { 0.255, 0.564, -0.478, -0.818, -0.043, 0.224, -0.268, 0.253, 0.021, 0.654, 0.98, -0.774 };
   int ldb = 2;
   float B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1780) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1780) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0, 0};
   float A[] = { -0.068, -0.603, -0.055, 0.14, 0.664, 0.987, 0.861, -0.691, -0.897, -0.778, 0.516, -0.073, -0.156, -0.42, 0.57, 0.628, 0.116, 0.344 };
   int lda = 3;
   float B[] = { 0.922, 0.39, -0.724, 0.421, 0.418, 0.92, -0.222, 0.835, 0.417, -0.392, 0.012, -0.346 };
   int ldb = 2;
   float B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1781) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1781) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { 0.904, 0.243, 0.206, 0.68, -0.946, 0.946, -0.675, 0.729 };
   int lda = 2;
   double B[] = { 0.427, 0.116, 0.916, -0.384, -0.372, -0.754, 0.148, 0.089, -0.924, 0.974, -0.307, -0.55 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1782) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1782) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { -0.898, 0.709, 0.719, -0.207, -0.841, -0.017, 0.202, -0.385 };
   int lda = 2;
   double B[] = { 0.308, 0.507, -0.838, 0.594, -0.811, 0.152, 0.118, -0.024, -0.632, 0.992, -0.942, 0.901 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1783) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1783) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { -0.849, 0.455, -0.273, -0.668, 0.196, -0.985, -0.39, 0.564 };
   int lda = 2;
   double B[] = { -0.874, 0.188, -0.039, 0.692, 0.33, 0.119, 0.012, 0.425, 0.787, -0.918, 0.739, -0.871 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1784) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1784) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { -0.325, 0.28, 0.902, -0.603, 0.091, -0.92, 0.209, -0.009 };
   int lda = 2;
   double B[] = { -0.202, -0.53, -0.88, -0.688, -0.215, 0.837, 0.917, 0.755, 0.477, 0.892, -0.524, -0.741 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1785) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1785) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { -0.756, 0.874, 0.56, 0.157, -0.831, -0.991, -0.531, 0.813 };
   int lda = 2;
   double B[] = { 0.271, 0.783, -0.861, 0.635, -0.088, 0.434, 0.256, -0.34, -0.724, -0.277, -0.604, 0.986 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1786) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1786) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { -0.371, -0.609, -0.812, -0.818, 0.45, -0.41, -0.704, -0.917 };
   int lda = 2;
   double B[] = { -0.268, 0.929, 0.82, 0.253, -0.883, 0.497, -0.265, 0.623, 0.131, -0.946, -0.365, 0.333 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1787) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1787) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { -0.265, 0.8, -0.676, -0.592, 0.78, -0.838, -0.651, 0.115 };
   int lda = 2;
   double B[] = { 0.942, 0.692, -0.516, 0.378, 0.028, 0.265, 0.289, -0.721, -0.25, -0.952, 0.463, -0.34 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1788) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1788) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { -0.852, -0.478, 0.16, 0.824, 0.073, 0.962, 0.509, -0.58 };
   int lda = 2;
   double B[] = { -0.789, 0.015, -0.779, -0.565, 0.048, -0.095, -0.272, 0.405, 0.272, 0.082, -0.693, -0.365 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1789) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1789) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { 0.251, 0.28, -0.092, 0.724, 0.928, -0.309, -0.222, -0.791, 0.113, -0.528, 0.148, 0.421, -0.833, 0.371, 0.354, 0.616, 0.313, 0.323 };
   int lda = 3;
   double B[] = { -0.769, -0.059, -0.068, 0.945, 0.938, -0.358, -0.17, 0.751, -0.248, -0.321, -0.818, 0.183 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1790) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1790) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { -0.707, -0.802, 0.13, -0.19, -0.564, -0.74, 0.118, -0.194, -0.124, -0.421, 0.665, 0.308, 0.505, -0.278, 0.588, 0.957, -0.727, 0.976 };
   int lda = 3;
   double B[] = { 0.153, -0.09, -0.4, 0.669, 0.689, -0.238, -0.259, 0.891, 0.993, 0.996, -0.829, -0.736 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1791) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1791) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { 0.83, 0.316, -0.099, 0.824, 0.767, 0.662, 0.244, 0.872, 0.35, 0.969, -0.084, 0.907, -0.752, -0.675, 0.129, -0.649, -0.539, 0.969 };
   int lda = 3;
   double B[] = { -0.145, 0.254, -0.497, -0.713, -0.742, 0.183, 0.272, -0.858, -0.606, -0.605, -0.807, 0.686 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1792) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1792) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { -0.091, 0.658, -0.834, -0.171, -0.126, -0.268, 0.879, -0.431, 0.678, -0.749, 0.136, -0.757, -0.578, 0.456, 0.978, -0.315, 0.333, 0.327 };
   int lda = 3;
   double B[] = { 0.963, -0.859, 0.599, 0.856, -0.924, 0.382, -0.531, 0.567, -0.454, 0.018, 0.97, 0.578 };
   int ldb = 3;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1793) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1793) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { -0.849, -0.819, 0.673, 0.574, -0.869, -0.969, -0.338, -0.097, -0.601, 0.903, 0.634, 0.313, 0.228, -0.028, 0.419, -0.762, 0.21, -0.532 };
   int lda = 3;
   double B[] = { -0.283, 0.999, -0.356, -0.459, 0.508, -0.132, -0.804, 0.173, 0.779, -0.427, 0.019, 0.347 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1794) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1794) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { -0.117, -0.663, -0.95, -0.273, -0.497, -0.037, 0.084, -0.831, 0.023, -0.241, 0.063, -0.023, -0.498, -0.137, -0.77, 0.457, -0.021, -0.69 };
   int lda = 3;
   double B[] = { 0.308, -0.004, 0.013, 0.354, 0.077, -0.944, -0.877, 0.741, -0.807, -0.3, 0.891, -0.056 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1795) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1795) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { -0.964, -0.653, 0.379, 0.994, -0.378, -0.409, 0.24, -0.333, 0.558, -0.099, -0.402, -0.812, 0.421, 0.823, -0.771, 0.998, 0.697, 0.253 };
   int lda = 3;
   double B[] = { 0.34, 0.479, 0.539, -0.133, 0.876, -0.347, 0.706, -0.623, 0.399, 0.903, -0.7, -0.088 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1796) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1796) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { -0.104, 0.643, -0.253, -0.988, -0.051, -0.805, 0.451, -0.421, -0.177, -0.534, -0.714, -0.581, -0.177, -0.582, -0.57, 0.259, -0.66, -0.864 };
   int lda = 3;
   double B[] = { 0.636, -0.365, -0.107, -0.279, 0.425, 0.976, 0.657, 0.294, 0.827, 0.187, 0.353, 0.31 };
   int ldb = 2;
   double B_expected[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1797) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1797) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.273, 0.812, 0.295, -0.415, -0.227, 0.901, 0.623, 0.786 };
   int lda = 2;
   double B[] = { -0.539, -0.551, -0.969, 0.09, -0.581, -0.594, -0.833, 0.457, -0.284, 0.434, -0.459, -0.662 };
   int ldb = 3;
   double B_expected[] = { -0.0312704, 0.2064538, 0.1775109, 0.1949157, -0.0337211, 0.2225517, 0.410638, -0.033917, 0.182384, -0.219409, 0.1257905, 0.1938415 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1798) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1798) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.323, 0.02, 0.718, 0.152, 0.665, 0.289, 0.317, 0.705 };
   int lda = 2;
   double B[] = { 0.448, -0.75, 0.851, 0.172, -0.244, 0.398, 0.602, 0.31, -0.017, 0.181, -0.119, 0.402 };
   int ldb = 3;
   double B_expected[] = { -0.0594, 0.2698, -0.2725, 0.0335, 0.0334, -0.1438, -0.2952588, 0.1518876, -0.213747, -0.073367, 0.0413388, -0.2306716 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1799) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1799) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.578, 0.018, -0.093, 0.964, 0.414, -0.729, 0.696, 0.874 };
   int lda = 2;
   double B[] = { -0.735, 0.788, -0.942, -0.71, -0.254, 0.265, 0.304, 0.218, 0.247, -0.172, 0.419, 0.448 };
   int ldb = 3;
   double B_expected[] = { -0.1486214, 0.2495598, -0.1744531, 0.0107667, -0.1648579, 0.1475263, -0.048058, -0.123122, -0.1062886, 0.0033742, -0.037823, -0.213397 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1800) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1800) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.358, -0.773, -0.065, 0.532, -0.319, 0.455, 0.578, 0.493 };
   int lda = 2;
   double B[] = { 0.744, -0.958, 0.162, 0.555, -0.131, 0.971, -0.467, 0.175, -0.794, 0.191, 0.361, 0.882 };
   int ldb = 3;
   double B_expected[] = { -0.1213734, 0.4492278, -0.1117944, -0.0070022, 0.108851, -0.320916, 0.1226, -0.0992, 0.2191, -0.1367, -0.1965, -0.2285 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1801) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1801) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.354, -0.504, -0.177, 0.186, -0.762, -0.506, 0.758, -0.994 };
   int lda = 2;
   double B[] = { -0.944, 0.562, 0.142, 0.742, 0.632, -0.627, -0.101, 0.476, 0.476, 0.675, 0.912, -0.33 };
   int ldb = 2;
   double B_expected[] = { -0.21291, -0.021306, -0.601736, 0.043676, 0.1715778, -0.0250026, 0.0587596, -0.2259812, -0.0036234, 0.1608258, 0.0885532, 0.6077736 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1802) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1802) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.001, 0.015, 0.942, 0.497, -0.104, 0.803, 0.679, 0.026 };
   int lda = 2;
   double B[] = { 0.889, -0.216, -0.912, -0.263, -0.329, 0.681, 0.332, -0.5, -0.484, 0.741, -0.728, -0.912 };
   int ldb = 2;
   double B_expected[] = { -0.2451, 0.1537, 0.2019693, -0.2251001, 0.0306, -0.2372, 0.1376892, 0.2324406, 0.0711, -0.2707, 0.5195777, 0.2860461 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1803) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1803) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.563, 0.394, -0.902, -0.27, 0.461, 0.939, -0.597, 0.803 };
   int lda = 2;
   double B[] = { 0.535, -0.111, 0.379, -0.036, 0.803, -0.341, 0.667, 0.001, 0.775, 0.714, 0.908, -0.508 };
   int ldb = 2;
   double B_expected[] = { 0.1623722, -0.1219324, 0.0266236, -0.1174842, 0.2429924, -0.1901218, 0.0662002, -0.2004014, 0.4905027, -0.2023089, -0.0629944, -0.3231352 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1804) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1804) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.159, 0.032, 0.785, 0.049, -0.128, 0.132, -0.735, -0.235 };
   int lda = 2;
   double B[] = { -0.331, -0.257, -0.725, 0.689, -0.793, 0.398, 0.127, -0.098, -0.498, -0.307, -0.019, 0.517 };
   int ldb = 2;
   double B_expected[] = { 0.2553318, -0.1678906, 0.1486, -0.2792, 0.1738216, -0.1670382, -0.0283, 0.0421, 0.151683, -0.083199, -0.046, -0.157 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1805) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1805) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.416, -0.424, -0.088, 0.614, -0.371, 0.983, -0.737, -0.647, 0.321, -0.518, 0.058, -0.533, 0.153, 0.283, 0.342, 0.993, -0.071, 0.225 };
   int lda = 3;
   double B[] = { -0.09, -0.844, -0.707, 0.903, 0.632, -0.294, -0.558, 0.74, -0.99, -0.855, -0.189, 0.543 };
   int ldb = 3;
   double B_expected[] = { 0.1668304, -0.2576208, -0.0664464, -0.0785782, -0.0226908, -0.0467944, -0.1091876, 0.3667652, 0.1076073, -0.1594011, 0.0407346, 0.0134478 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1806) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1806) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.67, -0.423, -0.165, 0.157, -0.43, 0.674, -0.35, 0.434, 0.972, -0.116, -0.029, 0.316, 0.914, 0.321, 0.132, 0.034, -0.907, -0.401 };
   int lda = 3;
   double B[] = { -0.396, 0.71, -0.588, 0.709, -0.024, -0.704, -0.988, 0.656, 0.665, -0.085, -0.778, 0.264 };
   int ldb = 3;
   double B_expected[] = { -0.1010812, -0.2287206, 0.0372688, -0.2530336, 0.0776, 0.2088, 0.264679, -0.133739, -0.147391, 0.161965, 0.207, -0.157 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1807) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1807) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.756, -0.149, -0.706, -0.162, -0.145, 0.67, 0.416, -0.27, -0.916, 0.995, -0.863, -0.25, -0.079, 0.248, -0.191, -0.195, 0.981, 0.834 };
   int lda = 3;
   double B[] = { 0.329, 0.921, -0.018, -0.02, 0.095, -0.892, -0.105, -0.799, -0.583, 0.564, -0.436, 0.965 };
   int ldb = 3;
   double B_expected[] = { -0.1805114, -0.1555812, -0.1560482, -0.0462226, -0.0967127, 0.2921239, 0.1183692, 0.1566766, 0.2260429, 0.3915667, 0.1788155, -0.2682995 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1808) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1808) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.552, -0.668, -0.013, 0.088, -0.766, 0.977, 0.088, -0.06, -0.311, 0.872, -0.328, -0.01, 0.659, -0.327, -0.276, 0.553, -0.734, -0.079 };
   int lda = 3;
   double B[] = { -0.87, 0.728, 0.997, -0.36, -0.046, -0.505, 0.082, -0.787, 0.414, 0.965, -0.048, -0.591 };
   int ldb = 3;
   double B_expected[] = { 0.1882, -0.3054, -0.2648624, 0.1695328, 0.0462155, -0.3187195, 0.0541, 0.2443, -0.2012812, -0.2298476, 0.3871505, 0.2622315 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1809) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1809) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.349, -0.072, 0.545, 0.212, -0.306, -0.009, 0.757, -0.925, 0.159, 0.308, 0.476, 0.1, 0.725, -0.757, -0.245, 0.571, 0.515, 0.993 };
   int lda = 3;
   double B[] = { 0.865, 0.501, 0.165, -0.63, -0.513, 0.351, -0.521, -0.062, 0.54, -0.634, -0.719, 0.216 };
   int ldb = 2;
   double B_expected[] = { -0.054193, 0.023274, 0.1487731, -0.3509657, -0.0481592, -0.1044386, 0.0666567, 0.1890461, -0.2932696, 0.0278532, 0.2357046, 0.1223408 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1810) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1810) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.941, -0.496, 0.492, 0.356, 0.353, 0.346, -0.519, -0.86, -0.677, -0.154, 0.313, 0.228, -0.56, -0.451, -0.78, 0.174, -0.663, 0.22 };
   int lda = 3;
   double B[] = { 0.162, -0.345, 0.188, 0.578, -0.675, 0.775, -0.018, 0.198, -0.222, -0.52, 0.672, -0.438 };
   int ldb = 2;
   double B_expected[] = { -0.3430472, 0.0394834, 0.0185782, -0.1505014, 0.0092108, -0.3837276, 0.0741276, -0.2435652, 0.1186, 0.1338, -0.1578, 0.1986 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1811) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1811) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.592, 0.708, 0.442, 0.212, 0.815, -0.638, 0.55, -0.512, -0.487, 0.181, 0.708, -0.126, 0.408, -0.51, 0.175, 0.114, -0.919, -0.268 };
   int lda = 3;
   double B[] = { 0.858, -0.004, 0.59, -0.395, -0.943, 0.824, 0.01, 0.455, -0.775, 0.062, -0.644, 0.03 };
   int ldb = 2;
   double B_expected[] = { -0.21374, -0.130452, -0.20707, 0.00773, -0.16787, 0.186571, -0.05026, 0.106515, -0.2887485, -0.0045065, -0.2446935, 0.1590455 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1812) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1812) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.988, -0.915, 0.963, 0.103, 0.921, 0.555, 0.846, 0.148, -0.43, 0.336, -0.371, 0.381, -0.487, 0.717, 0.881, -0.777, 0.774, -0.962 };
   int lda = 3;
   double B[] = { -0.805, 0.605, 0.481, 0.163, -0.057, -0.017, -0.886, 0.809, 0.875, 0.905, 0.095, 0.894 };
   int ldb = 2;
   double B_expected[] = { 0.181, -0.262, -0.1606, -0.0008, 0.220089, -0.234263, 0.0303246, -0.3486122, -0.0476352, -0.3174616, -0.2077412, -0.1552106 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1813) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1813) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.513, -0.385, -0.524, 0.726, 0.823, 0.839, -0.355, -0.881 };
   int lda = 2;
   double B[] = { -0.707, 0.016, 0.481, 0.935, 0.052, 0.719, 0.277, 0.169, 0.894, 0.352, -0.216, -0.741 };
   int ldb = 3;
   double B_expected[] = { -0.078919, 0.119774, 0.2114654, 0.0276682, 0.12593, 0.074299, -0.109352, -0.193196, 0.077864, 0.032876, -0.3330992, 0.2249494 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1814) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1814) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.136, -0.37, 0.669, -0.731, -0.4, 0.638, 0.833, -0.29 };
   int lda = 2;
   double B[] = { -0.861, -0.278, 0.941, 0.822, 0.88, 0.501, 0.911, -0.502, 0.573, -0.498, -0.517, -0.518 };
   int ldb = 3;
   double B_expected[] = { 0.2861, -0.0027, -0.3645, -0.1525, -0.3141, -0.0623, -0.0297254, 0.4490328, -0.254473, -0.161772, 0.0423084, -0.1675858 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1815) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1815) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.641, -0.058, 0.246, 0.884, -0.686, 0.123, -0.869, 0.891 };
   int lda = 2;
   double B[] = { 0.107, -0.333, 0.556, 0.124, 0.206, 0.049, -0.573, -0.9, -0.417, -0.734, -0.719, 0.76 };
   int ldb = 3;
   double B_expected[] = { -0.1591469, -0.1071617, -0.2301499, -0.1454657, -0.1758188, 0.1884616, -0.0380754, -0.4181892, -0.013453, -0.33198, -0.3886102, 0.1361404 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1816) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1816) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.083, 0.441, 0.995, 0.338, -0.988, -0.828, -0.254, -0.036 };
   int lda = 2;
   double B[] = { -0.792, 0.552, 0.033, -0.178, -0.225, 0.553, 0.348, 0.229, -0.151, -0.594, 0.711, -0.335 };
   int ldb = 3;
   double B_expected[] = { 0.3362416, -0.3167112, -0.2305904, -0.0177512, 0.0477576, -0.5068152, -0.1273, -0.0339, 0.1047, 0.1631, -0.1798, 0.1716 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1817) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1817) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.105, 0.584, -0.33, -0.182, -0.096, -0.257, 0.327, -0.123 };
   int lda = 2;
   double B[] = { -0.249, -0.274, -0.197, -0.899, 0.85, -0.318, 0.596, -0.237, 0.179, 0.046, -0.859, -0.459 };
   int ldb = 2;
   double B_expected[] = { 0.0441837, -0.0536099, -0.0065547, 0.1208159, 0.0819176, 0.1492908, -0.0917294, -0.0510192, -0.0037271, 0.0344777, 0.0974489, 0.0389047 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1818) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1818) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.972, 0.794, -0.968, -0.406, -0.2, -0.512, 0.436, 0.161 };
   int lda = 2;
   double B[] = { 0.817, -0.17, -0.613, -0.565, -0.494, 0.129, -0.593, -0.516, -0.695, -0.42, 0.848, 0.122 };
   int ldb = 2;
   double B_expected[] = { -0.2281, 0.1327, 0.2180776, -0.0351272, 0.1353, -0.0881, 0.2475472, 0.1823936, 0.2505, 0.0565, -0.345628, 0.165156 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1819) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1819) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.373, -0.316, -0.052, 0.025, -0.878, 0.612, 0.486, 0.953 };
   int lda = 2;
   double B[] = { -0.626, 0.408, 0.536, 0.66, -0.666, -0.127, 0.622, 0.036, -0.761, 0.773, -0.137, 0.074 };
   int ldb = 2;
   double B_expected[] = { 0.1214746, -0.0093742, -0.247838, 0.145962, 0.0994439, 0.0586017, -0.043453, 0.206241, 0.1510011, -0.0661437, -0.0178345, -0.0495635 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1820) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1820) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.621, -0.252, -0.942, 0.073, 0.416, -0.724, -0.972, 0.028 };
   int lda = 2;
   double B[] = { -0.006, 0.427, 0.292, -0.212, -0.319, -0.08, -0.401, 0.465, -0.493, -0.529, 0.003, -0.19 };
   int ldb = 2;
   double B_expected[] = { 0.0284232, -0.2112704, -0.0664, 0.0928, 0.0210696, 0.1558958, 0.0738, -0.1796, 0.1879327, 0.0541021, 0.0181, 0.0573 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1821) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1821) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.415, 0.215, 0.507, 0.094, 0.697, 0.633, 0.206, -0.383, -0.974, 0.734, -0.533, -0.15, -0.982, -0.232, -0.297, 0.501, -0.092, 0.663 };
   int lda = 3;
   double B[] = { 0.812, 0.323, 0.294, -0.423, -0.85, 0.043, -0.338, -0.568, 0.976, -0.375, 0.913, -0.119 };
   int ldb = 3;
   double B_expected[] = { 0.2153111, -0.0775367, 0.0404927, -0.0287599, -0.0879721, -0.1572073, -0.2481947, 0.2941819, 0.5234716, -0.1242382, 0.108305, 0.162022 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1822) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1822) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.827, -0.754, 0.719, 0.88, -0.942, -0.152, 0.051, 0.033, -0.603, -0.557, 0.668, 0.024, 0.082, 0.458, 0.733, 0.669, 0.722, -0.661 };
   int lda = 3;
   double B[] = { -0.523, 0.365, -0.811, -0.632, -0.06, 0.151, -0.962, -0.71, -0.543, 0.8, -0.264, 0.994 };
   int ldb = 3;
   double B_expected[] = { 0.4413193, -0.3047431, 0.307206, 0.074162, 0.0029, -0.0513, 0.2285887, 0.1349491, 0.061616, -0.510648, -0.0202, -0.3246 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1823) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1823) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.958, 0.948, -0.161, -0.34, -0.184, 0.43, -0.045, -0.465, -0.278, 0.461, 0.584, 0.003, -0.794, -0.778, -0.65, -0.91, 0.24, -0.944 };
   int lda = 3;
   double B[] = { 0.279, 0.041, -0.033, 0.332, 0.788, 0.611, -0.644, -0.133, 0.247, 0.06, 0.125, -0.407 };
   int ldb = 3;
   double B_expected[] = { -0.0693236, 0.0981792, -0.0442625, -0.0021815, 0.1936084, -0.3409328, 0.174601, -0.219233, 0.0274565, 0.1321885, -0.2252264, 0.1381888 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1824) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1824) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.983, -0.795, -0.115, -0.542, 0.837, 0.518, -0.164, 0.776, -0.453, -0.28, 0.135, -0.377, -0.199, -0.965, 0.784, -0.39, -0.499, 0.257 };
   int lda = 3;
   double B[] = { -0.712, 0.364, -0.28, 0.05, 0.314, 0.748, -0.719, 0.619, 0.474, -0.906, -0.859, 0.943 };
   int ldb = 3;
   double B_expected[] = { 0.1772, -0.1804, -0.0900512, -0.1509216, 0.0485292, 0.0109956, 0.1538, -0.2576, -0.2767208, 0.2420976, 0.2164354, 0.0610082 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1825) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1825) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.105, 0.503, -0.17, 0.2, -0.861, -0.279, -0.231, 0.058, 0.699, 0.437, 0.578, 0.462, 0.473, -0.793, -0.34, -0.162, -0.128, -0.844 };
   int lda = 3;
   double B[] = { -0.802, 0.292, -0.155, -0.916, -0.099, -0.082, 0.057, 0.215, 0.94, 0.911, -0.714, 0.41 };
   int ldb = 2;
   double B_expected[] = { -0.1044001, -0.5102243, 0.3865174, 0.0189802, 0.1888166, -0.0057672, -0.0800722, 0.0699214, 0.199086, -0.291946, 0.141904, 0.171064 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1826) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1826) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.468, 0.378, -0.498, 0.251, 0.777, -0.543, -0.913, 0.095, 0.779, -0.933, 0.068, -0.669, 0.715, 0.03, 0.012, 0.392, -0.785, -0.056 };
   int lda = 3;
   double B[] = { 0.143, -0.242, -0.379, -0.831, -0.46, -0.663, -0.735, -0.098, -0.861, -0.894, 0.772, -0.059 };
   int ldb = 2;
   double B_expected[] = { 0.0633681, 0.0476643, -0.1761819, 0.3044093, 0.2798556, 0.0187868, 0.2647924, 0.0455132, 0.3477, 0.1821, -0.2257, 0.0949 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1827) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1827) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.764, 0.908, 0.899, 0.119, -0.447, 0.279, 0.338, 0.73, -0.74, -0.366, -0.572, 0.583, 0.75, 0.519, 0.603, 0.831, 0.697, 0.822 };
   int lda = 3;
   double B[] = { 0.399, 0.572, -0.489, 0.964, -0.167, -0.104, 0.75, -0.199, 0.777, 0.503, -0.025, -0.386 };
   int ldb = 2;
   double B_expected[] = { 0.015568, 0.261244, -0.345424, 0.212636, -0.2247824, -0.0859342, 0.1074596, -0.4846822, -0.2415227, 0.2465939, 0.2042976, 0.2206978 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1828) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1828) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.432, 0.063, 0.065, -0.546, 0.099, 0.892, 0.48, -0.085, 0.746, -0.541, -0.739, -0.207, 0.695, 0.765, 0.197, -0.86, 0.621, -0.653 };
   int lda = 3;
   double B[] = { 0.182, 0.731, 0.571, 0.01, -0.357, -0.612, 0.581, 0.756, -0.911, -0.225, 0.438, 0.546 };
   int ldb = 2;
   double B_expected[] = { -0.1277, -0.2011, -0.1723, 0.0541, 0.2698001, 0.0651043, -0.2906381, -0.2592593, -0.0512125, -0.0040605, 0.0647965, 0.1119875 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1829) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1829) imag");
     };
   };
  };


}
