#include <gsl_test.h>
#include <gsl_math.h>
#include "test_cases.h"
#include "gsl_blas_raw.h"

const double eps   = 1.0e-15;
const double eps_f = 1.0e-07;


double frac_diff(double a, double b)
{
  return fabs(a-b)/(fabs(a) + fabs(b));
}


/* Level 1 tests */

int test_L1(void)
{
  int status = 0;
  int s;

  float  x_f;
  double x_d;

  float c[2];
  double z[2];


  /* sdsdot */

  x_f = gsl_blas_raw_sdsdot(4, 5.0, vector_4_f, 1, vector_4_f, 1);
  s = ( frac_diff(x_f, 14.0 + 5.0) > eps );
  gsl_test(s, "gsl_blas_raw_sdsdot A");
  status += s;

  x_f = gsl_blas_raw_sdsdot(2, 5.0, vector_4_f, 3, vector_4_f, 1);
  s = ( frac_diff(x_f, 1.0 + 5.0) > eps );
  gsl_test(s, "gsl_blas_raw_sdsdot B");
  status += s;


  /* dsdot */

  x_d = gsl_blas_raw_dsdot(4, vector_4_f, 1, vector_4_f, 1);
  s = ( frac_diff(x_d, 14.0) > eps );
  gsl_test(s, "gsl_blas_raw_dsdot A");
  status += s;

  x_d = gsl_blas_raw_dsdot(2, vector_4_f, 3, vector_4_f, 1);
  s = ( frac_diff(x_d, 1.0) > eps );
  gsl_test(s, "gsl_blas_raw_dsdot B");
  status += s;


  /* sdot */

  x_f = gsl_blas_raw_sdot(4, vector_4_f, 1, vector_4_f, 1);
  s = ( frac_diff(x_f, 14.0) > eps );
  gsl_test(s, "gsl_blas_raw_sdot A");
  status += s;

  x_f = gsl_blas_raw_sdot(2, vector_4_f, 3, vector_4_f, 1);
  s = ( frac_diff(x_f, 1.0) > eps );
  gsl_test(s, "gsl_blas_raw_sdot B");
  status += s;


  /* ddot */

  x_d = gsl_blas_raw_ddot(4, vector_4_d, 1, vector_4_d, 1);
  s = ( frac_diff(x_d, 14.0) > eps );
  gsl_test(s, "gsl_blas_raw_ddot A");
  status += s;

  x_d = gsl_blas_raw_ddot(2, vector_4_d, 3, vector_4_d, 1);
  s = ( frac_diff(x_d, 1.0) > eps );
  gsl_test(s, "gsl_blas_raw_ddot B");
  status += s;


  /* cdotu */

  gsl_blas_raw_cdotu(4, vector_4_c, 1, vector_4_c, 1, c);
  s = ( frac_diff(c[0], -22.0) > eps || frac_diff(c[1], 24.0) > eps );
  gsl_test(s, "gsl_blas_raw_cdotu A");
  status += s;

  gsl_blas_raw_cdotu(2, vector_4_c, 3, vector_4_c, 1, c);
  s = ( frac_diff(c[0], -5.0) > eps || frac_diff(c[1], -6.0) > eps );
  gsl_test(s, "gsl_blas_raw_cdotu B");
  status += s;


  /* cdotc */

  gsl_blas_raw_cdotc(4, vector_4_c, 1, vector_4_c, 1, c);
  s = ( frac_diff(c[0], 50.0) > eps || fabs(c[1]) > eps );
  gsl_test(s, "gsl_blas_raw_cdotc A");
  status += s;

  gsl_blas_raw_cdotc(2, vector_4_c, 3, vector_4_c, 1, c);
  s = ( frac_diff(c[0], 7.0) > eps || frac_diff(c[1], 8.0) > eps );
  gsl_test(s, "gsl_blas_raw_cdotc B");
  status += s;


  /* zdotu */

  gsl_blas_raw_zdotu(4, vector_4_z, 1, vector_4_z, 1, z);
  s = ( frac_diff(z[0], -22.0) > eps || frac_diff(z[1], 24.0) > eps );
  gsl_test(s, "gsl_blas_raw_zdotu A");
  status += s;

  gsl_blas_raw_zdotu(2, vector_4_z, 3, vector_4_z, 1, z);
  s = ( frac_diff(z[0], -5.0) > eps || frac_diff(z[1], -6.0) > eps );
  gsl_test(s, "gsl_blas_raw_zdotu B");
  status += s;


  /* zdotc */

  gsl_blas_raw_zdotc(4, vector_4_z, 1, vector_4_z, 1, z);
  s = ( frac_diff(z[0], 50.0) > eps || fabs(z[1]) > eps );
  gsl_test(s, "gsl_blas_raw_zdotc A");
  status += s;

  gsl_blas_raw_zdotc(2, vector_4_z, 3, vector_4_z, 1, z);
  s = ( frac_diff(z[0], 7.0) > eps || frac_diff(z[1], 8.0) > eps );
  gsl_test(s, "gsl_blas_raw_zdotc B");
  status += s;


  /* snrm2 */

  x_f = gsl_blas_raw_snrm2(4, vector_4_f, 1);
  s = ( frac_diff(x_f, sqrt(14.0)) > eps_f );
  gsl_test(s, "gsl_blas_raw_snrm2 A");
  status += s;

  x_f = gsl_blas_raw_snrm2(2, vector_4_f, 3);
  s = ( frac_diff(x_f, sqrt(13.0)) > eps_f );
  gsl_test(s, "gsl_blas_raw_snrm2 B");
  status += s;


  /* dnrm2 */

  x_d = gsl_blas_raw_dnrm2(4, vector_4_d, 1);
  s = ( frac_diff(x_d, sqrt(14.0)) > eps );
  gsl_test(s, "gsl_blas_raw_dnrm2 A");
  status += s;

  x_d = gsl_blas_raw_dnrm2(2, vector_4_d, 3);
  s = ( frac_diff(x_d, sqrt(13.0)) > eps );
  gsl_test(s, "gsl_blas_raw_dnrm2 B");
  status += s;

  /* scnrm2 */

  x_f = gsl_blas_raw_scnrm2(4, vector_4_c, 1);
  s = ( frac_diff(x_f, sqrt(50.0)) > eps_f );
  gsl_test(s, "gsl_blas_raw_scnrm2 A");
  status += s;

  x_f = gsl_blas_raw_scnrm2(2, vector_4_c, 3);
  s = ( frac_diff(x_f, sqrt(39.0)) > eps_f );
  gsl_test(s, "gsl_blas_raw_scnrm2 B");
  status += s;


  /* dznrm2 */

  x_d = gsl_blas_raw_dznrm2(4, vector_4_z, 1);
  s = ( frac_diff(x_d, sqrt(50.0)) > eps );
  gsl_test(s, "gsl_blas_raw_zdnrm2 A");
  status += s;

  x_d = gsl_blas_raw_dznrm2(2, vector_4_z, 3);
  s = ( frac_diff(x_d, sqrt(39.0)) > eps );
  gsl_test(s, "gsl_blas_raw_zdnrm2 B");
  status += s;  


  /* sasum */

  x_f = gsl_blas_raw_sasum(4, vector_4_f, 1);
  s = ( frac_diff(x_f, 6.0) > eps_f );
  gsl_test(s, "gsl_blas_raw_sasum A");
  status += s;

  x_f = gsl_blas_raw_sasum(2, vector_4_f, 3);
  s = ( frac_diff(x_f, 5.0) > eps_f );
  gsl_test(s, "gsl_blas_raw_sasum B");
  status += s;


  /* dasum */

  x_d = gsl_blas_raw_dasum(4, vector_4_d, 1);
  s = ( frac_diff(x_d, 6.0) > eps );
  gsl_test(s, "gsl_blas_raw_dasum A");
  status += s;

  x_d = gsl_blas_raw_dasum(2, vector_4_d, 3);
  s = ( frac_diff(x_d, 5.0) > eps );
  gsl_test(s, "gsl_blas_raw_dasum B");


  /* scasum */

  x_f = gsl_blas_raw_scasum(4, vector_4_c, 1);
  s = ( frac_diff(x_f, 16.0) > eps_f );
  gsl_test(s, "gsl_blas_raw_scasum A");
  status += s;

  x_f = gsl_blas_raw_scasum(2, vector_4_c, 3);
  s = ( frac_diff(x_f, 11.0) > eps_f );
  gsl_test(s, "gsl_blas_raw_scasum B");
  status += s;


  /* dzasum */

  x_d = gsl_blas_raw_dzasum(4, vector_4_z, 1);
  s = ( frac_diff(x_d, 16.0) > eps );
  gsl_test(s, "gsl_blas_raw_dzasum A");
  status += s;

  x_d = gsl_blas_raw_dzasum(2, vector_4_z, 3);
  s = ( frac_diff(x_d, 11.0) > eps );
  gsl_test(s, "gsl_blas_raw_dzasum B");
  status += s;



  return status;
}



int main()
{
  int status = 0;

  status += test_L1();

  return status;
}
