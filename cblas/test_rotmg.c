#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

void
test_rotmg () {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   float d1 = 0.00145765239243;
   float d2 = 0.000592069951704;
   float b1 = -0.585181586191;
   float b2 = 0.000438369659147;
   float h[] = { -999, -999.1, -999.2, -999.3, -999.4 };
   float d1_expected = 0.00145765206017;
   float d2_expected = 0.000592069816748;
   float b1_expected = -0.585181719577;
   float h0_expected = 0;
   float h11_expected = -999.1;
   float h21_expected = 0.000749117315875;
   float h12_expected = -0.000304276832621;
   float h22_expected = -999.4;
   cblas_srotmg(&d1, &d2, &b1, b2, h);
   gsl_test_rel(d1, d1_expected, flteps, "srotmg(case 464)");
   gsl_test_rel(d2, d2_expected, flteps, "srotmg(case 465)");
   gsl_test_rel(b1, b1_expected, flteps, "srotmg(case 466)");
   gsl_test_rel(h[0], h0_expected, flteps, "srotmg(case 467)");
   gsl_test_rel(h[1], h11_expected, flteps, "srotmg(case 468)");
   gsl_test_rel(h[2], h21_expected, flteps, "srotmg(case 469)");
   gsl_test_rel(h[3], h12_expected, flteps, "srotmg(case 470)");
   gsl_test_rel(h[4], h22_expected, flteps, "srotmg(case 471)");
  };


  {
   double d1 = 7.01947923253;
   double d2 = -0.00116963299111;
   double b1 = -0.0267926523512;
   double b2 = -0.0256880295441;
   double h[] = { -999, -999.1, -999.2, -999.3, -999.4 };
   double d1_expected = 7.02055457378;
   double d2_expected = -0.00116981217172;
   double b1_expected = -0.0267885485096;
   double h0_expected = 0;
   double h11_expected = -999.1;
   double h21_expected = -0.958771427603;
   double h12_expected = -0.000159756964229;
   double h22_expected = -999.4;
   cblas_drotmg(&d1, &d2, &b1, b2, h);
   gsl_test_rel(d1, d1_expected, dbleps, "drotmg(case 472)");
   gsl_test_rel(d2, d2_expected, dbleps, "drotmg(case 473)");
   gsl_test_rel(b1, b1_expected, dbleps, "drotmg(case 474)");
   gsl_test_rel(h[0], h0_expected, dbleps, "drotmg(case 475)");
   gsl_test_rel(h[1], h11_expected, dbleps, "drotmg(case 476)");
   gsl_test_rel(h[2], h21_expected, dbleps, "drotmg(case 477)");
   gsl_test_rel(h[3], h12_expected, dbleps, "drotmg(case 478)");
   gsl_test_rel(h[4], h22_expected, dbleps, "drotmg(case 479)");
  };


}
