#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_scal () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int N = 1;
   float alpha = 0;
   float X[] = { 0.651 };
   int incX = -1;
   float expected[] = { 0.651 };
   cblas_sscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], expected[i], flteps, "sscal(case 112)");
     }
   };
  };


  {
   int N = 1;
   float alpha = 0.1;
   float X[] = { 0.651 };
   int incX = -1;
   float expected[] = { 0.651 };
   cblas_sscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], expected[i], flteps, "sscal(case 113)");
     }
   };
  };


  {
   int N = 1;
   float alpha = 1;
   float X[] = { 0.651 };
   int incX = -1;
   float expected[] = { 0.651 };
   cblas_sscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], expected[i], flteps, "sscal(case 114)");
     }
   };
  };


  {
   int N = 1;
   double alpha = 0;
   double X[] = { 0.686 };
   int incX = -1;
   double expected[] = { 0.686 };
   cblas_dscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], expected[i], dbleps, "dscal(case 115)");
     }
   };
  };


  {
   int N = 1;
   double alpha = 0.1;
   double X[] = { 0.686 };
   int incX = -1;
   double expected[] = { 0.686 };
   cblas_dscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], expected[i], dbleps, "dscal(case 116)");
     }
   };
  };


  {
   int N = 1;
   double alpha = 1;
   double X[] = { 0.686 };
   int incX = -1;
   double expected[] = { 0.686 };
   cblas_dscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], expected[i], dbleps, "dscal(case 117)");
     }
   };
  };


  {
   int N = 1;
   float alpha[2] = {0, 0};
   float X[] = { 0.986, -0.775 };
   int incX = -1;
   float expected[] = { 0.986, -0.775 };
   cblas_cscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], expected[2*i], flteps, "cscal(case 118) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], flteps, "cscal(case 118) imag");
     };
   };
  };


  {
   int N = 1;
   float alpha[2] = {0.1, 0};
   float X[] = { 0.986, -0.775 };
   int incX = -1;
   float expected[] = { 0.986, -0.775 };
   cblas_cscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], expected[2*i], flteps, "cscal(case 119) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], flteps, "cscal(case 119) imag");
     };
   };
  };


  {
   int N = 1;
   float alpha[2] = {1, 0};
   float X[] = { 0.986, -0.775 };
   int incX = -1;
   float expected[] = { 0.986, -0.775 };
   cblas_cscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], expected[2*i], flteps, "cscal(case 120) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], flteps, "cscal(case 120) imag");
     };
   };
  };


  {
   int N = 1;
   float alpha[2] = {0, 0.1};
   float X[] = { 0.986, -0.775 };
   int incX = -1;
   float expected[] = { 0.986, -0.775 };
   cblas_cscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], expected[2*i], flteps, "cscal(case 121) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], flteps, "cscal(case 121) imag");
     };
   };
  };


  {
   int N = 1;
   float alpha[2] = {0.1, 0.2};
   float X[] = { 0.986, -0.775 };
   int incX = -1;
   float expected[] = { 0.986, -0.775 };
   cblas_cscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], expected[2*i], flteps, "cscal(case 122) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], flteps, "cscal(case 122) imag");
     };
   };
  };


  {
   int N = 1;
   float alpha[2] = {1, 0.3};
   float X[] = { 0.986, -0.775 };
   int incX = -1;
   float expected[] = { 0.986, -0.775 };
   cblas_cscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], expected[2*i], flteps, "cscal(case 123) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], flteps, "cscal(case 123) imag");
     };
   };
  };


  {
   int N = 1;
   double alpha[2] = {0, 0};
   double X[] = { 0.454, -0.478 };
   int incX = -1;
   double expected[] = { 0.454, -0.478 };
   cblas_zscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], expected[2*i], dbleps, "zscal(case 124) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], dbleps, "zscal(case 124) imag");
     };
   };
  };


  {
   int N = 1;
   double alpha[2] = {0.1, 0};
   double X[] = { 0.454, -0.478 };
   int incX = -1;
   double expected[] = { 0.454, -0.478 };
   cblas_zscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], expected[2*i], dbleps, "zscal(case 125) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], dbleps, "zscal(case 125) imag");
     };
   };
  };


  {
   int N = 1;
   double alpha[2] = {1, 0};
   double X[] = { 0.454, -0.478 };
   int incX = -1;
   double expected[] = { 0.454, -0.478 };
   cblas_zscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], expected[2*i], dbleps, "zscal(case 126) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], dbleps, "zscal(case 126) imag");
     };
   };
  };


  {
   int N = 1;
   double alpha[2] = {0, 0.1};
   double X[] = { 0.454, -0.478 };
   int incX = -1;
   double expected[] = { 0.454, -0.478 };
   cblas_zscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], expected[2*i], dbleps, "zscal(case 127) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], dbleps, "zscal(case 127) imag");
     };
   };
  };


  {
   int N = 1;
   double alpha[2] = {0.1, 0.2};
   double X[] = { 0.454, -0.478 };
   int incX = -1;
   double expected[] = { 0.454, -0.478 };
   cblas_zscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], expected[2*i], dbleps, "zscal(case 128) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], dbleps, "zscal(case 128) imag");
     };
   };
  };


  {
   int N = 1;
   double alpha[2] = {1, 0.3};
   double X[] = { 0.454, -0.478 };
   int incX = -1;
   double expected[] = { 0.454, -0.478 };
   cblas_zscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[2*i], expected[2*i], dbleps, "zscal(case 129) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], dbleps, "zscal(case 129) imag");
     };
   };
  };


  {
   int N = 2;
   float alpha = 0;
   float X[] = { 0.389, -0.236 };
   int incX = 1;
   float expected[] = { 0, -0 };
   cblas_sscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], expected[i], flteps, "sscal(case 130)");
     }
   };
  };


  {
   int N = 2;
   float alpha = 0.1;
   float X[] = { 0.389, -0.236 };
   int incX = 1;
   float expected[] = { 0.0389, -0.0236 };
   cblas_sscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], expected[i], flteps, "sscal(case 131)");
     }
   };
  };


  {
   int N = 2;
   float alpha = 1;
   float X[] = { 0.389, -0.236 };
   int incX = 1;
   float expected[] = { 0.389, -0.236 };
   cblas_sscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], expected[i], flteps, "sscal(case 132)");
     }
   };
  };


  {
   int N = 2;
   double alpha = 0;
   double X[] = { -0.429, -0.183 };
   int incX = 1;
   double expected[] = { -0, -0 };
   cblas_dscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], expected[i], dbleps, "dscal(case 133)");
     }
   };
  };


  {
   int N = 2;
   double alpha = 0.1;
   double X[] = { -0.429, -0.183 };
   int incX = 1;
   double expected[] = { -0.0429, -0.0183 };
   cblas_dscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], expected[i], dbleps, "dscal(case 134)");
     }
   };
  };


  {
   int N = 2;
   double alpha = 1;
   double X[] = { -0.429, -0.183 };
   int incX = 1;
   double expected[] = { -0.429, -0.183 };
   cblas_dscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], expected[i], dbleps, "dscal(case 135)");
     }
   };
  };


  {
   int N = 2;
   float alpha[2] = {0, 0};
   float X[] = { -0.603, 0.239, 0.339, -0.58 };
   int incX = 1;
   float expected[] = { -0, 0, 0, 0 };
   cblas_cscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], expected[2*i], flteps, "cscal(case 136) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], flteps, "cscal(case 136) imag");
     };
   };
  };


  {
   int N = 2;
   float alpha[2] = {0.1, 0};
   float X[] = { -0.603, 0.239, 0.339, -0.58 };
   int incX = 1;
   float expected[] = { -0.0603, 0.0239, 0.0339, -0.058 };
   cblas_cscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], expected[2*i], flteps, "cscal(case 137) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], flteps, "cscal(case 137) imag");
     };
   };
  };


  {
   int N = 2;
   float alpha[2] = {1, 0};
   float X[] = { -0.603, 0.239, 0.339, -0.58 };
   int incX = 1;
   float expected[] = { -0.603, 0.239, 0.339, -0.58 };
   cblas_cscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], expected[2*i], flteps, "cscal(case 138) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], flteps, "cscal(case 138) imag");
     };
   };
  };


  {
   int N = 2;
   float alpha[2] = {0, 0.1};
   float X[] = { -0.603, 0.239, 0.339, -0.58 };
   int incX = 1;
   float expected[] = { -0.0239, -0.0603, 0.058, 0.0339 };
   cblas_cscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], expected[2*i], flteps, "cscal(case 139) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], flteps, "cscal(case 139) imag");
     };
   };
  };


  {
   int N = 2;
   float alpha[2] = {0.1, 0.2};
   float X[] = { -0.603, 0.239, 0.339, -0.58 };
   int incX = 1;
   float expected[] = { -0.1081, -0.0967, 0.1499, 0.0098 };
   cblas_cscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], expected[2*i], flteps, "cscal(case 140) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], flteps, "cscal(case 140) imag");
     };
   };
  };


  {
   int N = 2;
   float alpha[2] = {1, 0.3};
   float X[] = { -0.603, 0.239, 0.339, -0.58 };
   int incX = 1;
   float expected[] = { -0.6747, 0.0581, 0.513, -0.4783 };
   cblas_cscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], expected[2*i], flteps, "cscal(case 141) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], flteps, "cscal(case 141) imag");
     };
   };
  };


  {
   int N = 2;
   double alpha[2] = {0, 0};
   double X[] = { -0.956, 0.613, 0.443, 0.503 };
   int incX = 1;
   double expected[] = { -0, 0, 0, 0 };
   cblas_zscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], expected[2*i], dbleps, "zscal(case 142) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], dbleps, "zscal(case 142) imag");
     };
   };
  };


  {
   int N = 2;
   double alpha[2] = {0.1, 0};
   double X[] = { -0.956, 0.613, 0.443, 0.503 };
   int incX = 1;
   double expected[] = { -0.0956, 0.0613, 0.0443, 0.0503 };
   cblas_zscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], expected[2*i], dbleps, "zscal(case 143) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], dbleps, "zscal(case 143) imag");
     };
   };
  };


  {
   int N = 2;
   double alpha[2] = {1, 0};
   double X[] = { -0.956, 0.613, 0.443, 0.503 };
   int incX = 1;
   double expected[] = { -0.956, 0.613, 0.443, 0.503 };
   cblas_zscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], expected[2*i], dbleps, "zscal(case 144) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], dbleps, "zscal(case 144) imag");
     };
   };
  };


  {
   int N = 2;
   double alpha[2] = {0, 0.1};
   double X[] = { -0.956, 0.613, 0.443, 0.503 };
   int incX = 1;
   double expected[] = { -0.0613, -0.0956, -0.0503, 0.0443 };
   cblas_zscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], expected[2*i], dbleps, "zscal(case 145) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], dbleps, "zscal(case 145) imag");
     };
   };
  };


  {
   int N = 2;
   double alpha[2] = {0.1, 0.2};
   double X[] = { -0.956, 0.613, 0.443, 0.503 };
   int incX = 1;
   double expected[] = { -0.2182, -0.1299, -0.0563, 0.1389 };
   cblas_zscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], expected[2*i], dbleps, "zscal(case 146) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], dbleps, "zscal(case 146) imag");
     };
   };
  };


  {
   int N = 2;
   double alpha[2] = {1, 0.3};
   double X[] = { -0.956, 0.613, 0.443, 0.503 };
   int incX = 1;
   double expected[] = { -1.1399, 0.3262, 0.2921, 0.6359 };
   cblas_zscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], expected[2*i], dbleps, "zscal(case 147) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], dbleps, "zscal(case 147) imag");
     };
   };
  };


  {
   int N = 2;
   float alpha = 0;
   float X[] = { 0.629, -0.419 };
   int incX = -1;
   float expected[] = { 0.629, -0.419 };
   cblas_sscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], expected[i], flteps, "sscal(case 148)");
     }
   };
  };


  {
   int N = 2;
   float alpha = 0.1;
   float X[] = { 0.629, -0.419 };
   int incX = -1;
   float expected[] = { 0.629, -0.419 };
   cblas_sscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], expected[i], flteps, "sscal(case 149)");
     }
   };
  };


  {
   int N = 2;
   float alpha = 1;
   float X[] = { 0.629, -0.419 };
   int incX = -1;
   float expected[] = { 0.629, -0.419 };
   cblas_sscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], expected[i], flteps, "sscal(case 150)");
     }
   };
  };


  {
   int N = 2;
   double alpha = 0;
   double X[] = { 0.398, -0.656 };
   int incX = -1;
   double expected[] = { 0.398, -0.656 };
   cblas_dscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], expected[i], dbleps, "dscal(case 151)");
     }
   };
  };


  {
   int N = 2;
   double alpha = 0.1;
   double X[] = { 0.398, -0.656 };
   int incX = -1;
   double expected[] = { 0.398, -0.656 };
   cblas_dscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], expected[i], dbleps, "dscal(case 152)");
     }
   };
  };


  {
   int N = 2;
   double alpha = 1;
   double X[] = { 0.398, -0.656 };
   int incX = -1;
   double expected[] = { 0.398, -0.656 };
   cblas_dscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[i], expected[i], dbleps, "dscal(case 153)");
     }
   };
  };


  {
   int N = 2;
   float alpha[2] = {0, 0};
   float X[] = { 0.736, 0.331, -0.318, 0.622 };
   int incX = -1;
   float expected[] = { 0.736, 0.331, -0.318, 0.622 };
   cblas_cscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], expected[2*i], flteps, "cscal(case 154) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], flteps, "cscal(case 154) imag");
     };
   };
  };


  {
   int N = 2;
   float alpha[2] = {0.1, 0};
   float X[] = { 0.736, 0.331, -0.318, 0.622 };
   int incX = -1;
   float expected[] = { 0.736, 0.331, -0.318, 0.622 };
   cblas_cscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], expected[2*i], flteps, "cscal(case 155) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], flteps, "cscal(case 155) imag");
     };
   };
  };


  {
   int N = 2;
   float alpha[2] = {1, 0};
   float X[] = { 0.736, 0.331, -0.318, 0.622 };
   int incX = -1;
   float expected[] = { 0.736, 0.331, -0.318, 0.622 };
   cblas_cscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], expected[2*i], flteps, "cscal(case 156) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], flteps, "cscal(case 156) imag");
     };
   };
  };


  {
   int N = 2;
   float alpha[2] = {0, 0.1};
   float X[] = { 0.736, 0.331, -0.318, 0.622 };
   int incX = -1;
   float expected[] = { 0.736, 0.331, -0.318, 0.622 };
   cblas_cscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], expected[2*i], flteps, "cscal(case 157) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], flteps, "cscal(case 157) imag");
     };
   };
  };


  {
   int N = 2;
   float alpha[2] = {0.1, 0.2};
   float X[] = { 0.736, 0.331, -0.318, 0.622 };
   int incX = -1;
   float expected[] = { 0.736, 0.331, -0.318, 0.622 };
   cblas_cscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], expected[2*i], flteps, "cscal(case 158) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], flteps, "cscal(case 158) imag");
     };
   };
  };


  {
   int N = 2;
   float alpha[2] = {1, 0.3};
   float X[] = { 0.736, 0.331, -0.318, 0.622 };
   int incX = -1;
   float expected[] = { 0.736, 0.331, -0.318, 0.622 };
   cblas_cscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], expected[2*i], flteps, "cscal(case 159) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], flteps, "cscal(case 159) imag");
     };
   };
  };


  {
   int N = 2;
   double alpha[2] = {0, 0};
   double X[] = { 0.521, -0.811, 0.556, -0.147 };
   int incX = -1;
   double expected[] = { 0.521, -0.811, 0.556, -0.147 };
   cblas_zscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], expected[2*i], dbleps, "zscal(case 160) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], dbleps, "zscal(case 160) imag");
     };
   };
  };


  {
   int N = 2;
   double alpha[2] = {0.1, 0};
   double X[] = { 0.521, -0.811, 0.556, -0.147 };
   int incX = -1;
   double expected[] = { 0.521, -0.811, 0.556, -0.147 };
   cblas_zscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], expected[2*i], dbleps, "zscal(case 161) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], dbleps, "zscal(case 161) imag");
     };
   };
  };


  {
   int N = 2;
   double alpha[2] = {1, 0};
   double X[] = { 0.521, -0.811, 0.556, -0.147 };
   int incX = -1;
   double expected[] = { 0.521, -0.811, 0.556, -0.147 };
   cblas_zscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], expected[2*i], dbleps, "zscal(case 162) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], dbleps, "zscal(case 162) imag");
     };
   };
  };


  {
   int N = 2;
   double alpha[2] = {0, 0.1};
   double X[] = { 0.521, -0.811, 0.556, -0.147 };
   int incX = -1;
   double expected[] = { 0.521, -0.811, 0.556, -0.147 };
   cblas_zscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], expected[2*i], dbleps, "zscal(case 163) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], dbleps, "zscal(case 163) imag");
     };
   };
  };


  {
   int N = 2;
   double alpha[2] = {0.1, 0.2};
   double X[] = { 0.521, -0.811, 0.556, -0.147 };
   int incX = -1;
   double expected[] = { 0.521, -0.811, 0.556, -0.147 };
   cblas_zscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], expected[2*i], dbleps, "zscal(case 164) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], dbleps, "zscal(case 164) imag");
     };
   };
  };


  {
   int N = 2;
   double alpha[2] = {1, 0.3};
   double X[] = { 0.521, -0.811, 0.556, -0.147 };
   int incX = -1;
   double expected[] = { 0.521, -0.811, 0.556, -0.147 };
   cblas_zscal(N, alpha, X, incX);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(X[2*i], expected[2*i], dbleps, "zscal(case 165) real");
       gsl_test_rel(X[2*i+1], expected[2*i+1], dbleps, "zscal(case 165) imag");
     };
   };
  };


}
