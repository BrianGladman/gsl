#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_rotm () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int N = 1;
   float h[] = { -1, -5.75980058783, 3.321447378776e+03, 0.011381907396, -0.00100289944225 };
   float X[] = { -0.927 };
   int incX = 1;
   float Y[] = { -0.554 };
   int incY = -1;
   float x_expected[] = { 5.33302956822 };
   float y_expected[] = { -3.078981164519e+03 };
   cblas_srotm(N, X, incX, Y, incY, h);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "srotm(case 480)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "srotm(case 481)");
     }
   };
  };


  {
   int N = 1;
   float h[] = { 0.000000000000e+00, 3.538824295371e+04, 5.709241129910e-04, -22.6895696448, 38.7214117842 };
   float X[] = { -0.927 };
   int incX = 1;
   float Y[] = { -0.554 };
   int incY = -1;
   float x_expected[] = { 11.6430215832 };
   float y_expected[] = { -0.554529246653 };
   cblas_srotm(N, X, incX, Y, incY, h);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "srotm(case 482)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "srotm(case 483)");
     }
   };
  };


  {
   int N = 1;
   float h[] = { 1, 1.129608743775e-04, 275.291846192, 57.6202549144, 485.94817593 };
   float X[] = { -0.927 };
   int incX = 1;
   float Y[] = { -0.554 };
   int incY = -1;
   float x_expected[] = { -0.554104714731 };
   float y_expected[] = { -268.288289465 };
   cblas_srotm(N, X, incX, Y, incY, h);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "srotm(case 484)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "srotm(case 485)");
     }
   };
  };


  {
   int N = 1;
   float h[] = { -2, -9.259319704766e+04, -3.685997388537e+04, 7.25806759797, 0.00658117449935 };
   float X[] = { -0.927 };
   int incX = 1;
   float Y[] = { -0.554 };
   int incY = -1;
   float x_expected[] = { -0.927 };
   float y_expected[] = { -0.554 };
   cblas_srotm(N, X, incX, Y, incY, h);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "srotm(case 486)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "srotm(case 487)");
     }
   };
  };


  {
   int N = 1;
   float h[] = { -1, -1.2422321269, -3.612766930622e-05, 4.491993173898e-04, -0.00194090979621 };
   float X[] = { -0.927 };
   int incX = 1;
   float Y[] = { -0.554 };
   int incY = -1;
   float x_expected[] = { 1.15130032522 };
   float y_expected[] = { 0.00110875437655 };
   cblas_srotm(N, X, incX, Y, incY, h);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "srotm(case 488)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "srotm(case 489)");
     }
   };
  };


  {
   int N = 1;
   float h[] = { 0.000000000000e+00, -0.00502428748094, -0.00269096326389, 0.0540388390637, 5.243917751912e+03 };
   float X[] = { -0.927 };
   int incX = 1;
   float Y[] = { -0.554 };
   int incY = -1;
   float x_expected[] = { -0.956937516841 };
   float y_expected[] = { -0.551505477054 };
   cblas_srotm(N, X, incX, Y, incY, h);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "srotm(case 490)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "srotm(case 491)");
     }
   };
  };


  {
   int N = 1;
   float h[] = { 1, 7.222773083206e+03, 1.094492130176e-05, 0.0102788444513, -9.40292098981 };
   float X[] = { -0.927 };
   int incX = 1;
   float Y[] = { -0.554 };
   int incY = -1;
   float x_expected[] = { -6.696064648132e+03 };
   float y_expected[] = { 6.13621822835 };
   cblas_srotm(N, X, incX, Y, incY, h);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "srotm(case 492)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "srotm(case 493)");
     }
   };
  };


  {
   int N = 1;
   float h[] = { -2, 3.04579892527, 4.387338021520e-04, 1.739032938502e+04, 104.97320748 };
   float X[] = { -0.927 };
   int incX = 1;
   float Y[] = { -0.554 };
   int incY = -1;
   float x_expected[] = { -0.927 };
   float y_expected[] = { -0.554 };
   cblas_srotm(N, X, incX, Y, incY, h);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "srotm(case 494)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "srotm(case 495)");
     }
   };
  };


  {
   int N = 1;
   float h[] = { -1, 1.727658381567e+03, 34.6979983774, -0.939619472089, -3.149966367331e+04 };
   float X[] = { -0.927 };
   int incX = 1;
   float Y[] = { -0.554 };
   int incY = -1;
   float x_expected[] = { -1.601018770525e+03 };
   float y_expected[] = { 1.741864863052e+04 };
   cblas_srotm(N, X, incX, Y, incY, h);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "srotm(case 496)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "srotm(case 497)");
     }
   };
  };


  {
   int N = 1;
   float h[] = { 0.000000000000e+00, -5.99683738995, 6.958113482311e+03, -8.232107101893e+04, 4.586037806461e+03 };
   float X[] = { -0.927 };
   int incX = 1;
   float Y[] = { -0.554 };
   int incY = -1;
   float x_expected[] = { 4.560494634449e+04 };
   float y_expected[] = { -6.450725198102e+03 };
   cblas_srotm(N, X, incX, Y, incY, h);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "srotm(case 498)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "srotm(case 499)");
     }
   };
  };


  {
   int N = 1;
   double h[] = { -1, -1.02995859704, 2.082161760569e+04, 0.147405205715, -37.8310845634 };
   double X[] = { 0.898 };
   int incX = 1;
   double Y[] = { 0.699 };
   int incY = -1;
   double x_expected[] = { -0.821866581348 };
   double y_expected[] = { 1.867136868180e+04 };
   cblas_drotm(N, X, incX, Y, incY, h);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "drotm(case 500)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "drotm(case 501)");
     }
   };
  };


  {
   int N = 1;
   double h[] = { 0.000000000000e+00, -0.0274034457938, 1.28224974662, 12.9274829091, 76.2938135407 };
   double X[] = { 0.898 };
   int incX = 1;
   double Y[] = { 0.699 };
   int incY = -1;
   double x_expected[] = { 9.93431055348 };
   double y_expected[] = { 1.85046027246 };
   cblas_drotm(N, X, incX, Y, incY, h);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "drotm(case 502)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "drotm(case 503)");
     }
   };
  };


  {
   int N = 1;
   double h[] = { 1, -5.2442689774, 0.164179586168, -0.0471943419602, -2.36823965913 };
   double X[] = { 0.898 };
   int incX = 1;
   double Y[] = { 0.699 };
   int incY = -1;
   double x_expected[] = { -4.01035354171 };
   double y_expected[] = { -2.55339952173 };
   cblas_drotm(N, X, incX, Y, incY, h);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "drotm(case 504)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "drotm(case 505)");
     }
   };
  };


  {
   int N = 1;
   double h[] = { -2, -0.0133921537363, -35.5540447903, 2.65138334384, -1.435758005093e+04 };
   double X[] = { 0.898 };
   int incX = 1;
   double Y[] = { 0.699 };
   int incY = -1;
   double x_expected[] = { 0.898 };
   double y_expected[] = { 0.699 };
   cblas_drotm(N, X, incX, Y, incY, h);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "drotm(case 506)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "drotm(case 507)");
     }
   };
  };


  {
   int N = 1;
   double h[] = { -1, 7.922514040997e-04, 3.425634659628e+04, -497.787532255, -0.0261220542513 };
   double X[] = { 0.898 };
   int incX = 1;
   double Y[] = { 0.699 };
   int incY = -1;
   double x_expected[] = { -347.952773604 };
   double y_expected[] = { 3.076218098414e+04 };
   cblas_drotm(N, X, incX, Y, incY, h);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "drotm(case 508)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "drotm(case 509)");
     }
   };
  };


  {
   int N = 1;
   double h[] = { 0.000000000000e+00, 0.0272191223261, -0.0743939716696, -367.115537691, -47.9630083232 };
   double X[] = { 0.898 };
   int incX = 1;
   double Y[] = { 0.699 };
   int incY = -1;
   double x_expected[] = { -255.715760846 };
   double y_expected[] = { 0.632194213441 };
   cblas_drotm(N, X, incX, Y, incY, h);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "drotm(case 510)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "drotm(case 511)");
     }
   };
  };


  {
   int N = 1;
   double h[] = { 1, -6.645488820841e+03, 5.60330418541, 0.063811792489, -0.13931093654 };
   double X[] = { 0.898 };
   int incX = 1;
   double Y[] = { 0.699 };
   int incY = -1;
   double x_expected[] = { -5.966949961115e+03 };
   double y_expected[] = { -0.995378344642 };
   cblas_drotm(N, X, incX, Y, incY, h);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "drotm(case 512)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "drotm(case 513)");
     }
   };
  };


  {
   int N = 1;
   double h[] = { -2, -298.760433562, -1.565597465631e-04, -0.0094402167, 1.849385201220e+04 };
   double X[] = { 0.898 };
   int incX = 1;
   double Y[] = { 0.699 };
   int incY = -1;
   double x_expected[] = { 0.898 };
   double y_expected[] = { 0.699 };
   cblas_drotm(N, X, incX, Y, incY, h);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "drotm(case 514)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "drotm(case 515)");
     }
   };
  };


  {
   int N = 1;
   double h[] = { -1, -1.11056328843, -8.791376276919e+04, 0.0723081690309, 1.816755222220e+03 };
   double X[] = { 0.898 };
   int incX = 1;
   double Y[] = { 0.699 };
   int incY = -1;
   double x_expected[] = { -0.946742422857 };
   double y_expected[] = { -7.767664706640e+04 };
   cblas_drotm(N, X, incX, Y, incY, h);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "drotm(case 516)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "drotm(case 517)");
     }
   };
  };


  {
   int N = 1;
   double h[] = { 0.000000000000e+00, -8.579562121919e+04, -88.5040351412, -187.285164165, -0.00712906621648 };
   double X[] = { 0.898 };
   int incX = 1;
   double Y[] = { 0.699 };
   int incY = -1;
   double x_expected[] = { -130.014329752 };
   double y_expected[] = { -78.7776235568 };
   cblas_drotm(N, X, incX, Y, incY, h);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "drotm(case 518)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "drotm(case 519)");
     }
   };
  };


}
