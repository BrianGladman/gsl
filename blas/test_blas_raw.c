/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */

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


/**** level 1 tests ****/

int test_L1(void)
{
  int status = 0;
  int s;

  float  x_f;
  double x_d;
  CBLAS_INDEX bi;

  float c[2];
  double z[2];

  float  tmp_a_f[256];
  float  tmp_b_f[256];
  double tmp_a_d[256];
  double tmp_b_d[256];

  float  alpha_c[2] = { 1.0, 1.0 };
  double alpha_z[2] = { 1.0, 1.0 };


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
  s = ( frac_diff(x_f, 14.0) > eps_f );
  gsl_test(s, "gsl_blas_raw_sdot A");
  status += s;

  x_f = gsl_blas_raw_sdot(2, vector_4_f, 3, vector_4_f, 1);
  s = ( frac_diff(x_f, 1.0) > eps_f );
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
  s = ( frac_diff(c[0], -22.0) > eps_f || frac_diff(c[1], 24.0) > eps_f );
  gsl_test(s, "gsl_blas_raw_cdotu A");
  status += s;

  gsl_blas_raw_cdotu(2, vector_4_c, 3, vector_4_c, 1, c);
  s = ( frac_diff(c[0], -5.0) > eps_f || frac_diff(c[1], -6.0) > eps_f );
  gsl_test(s, "gsl_blas_raw_cdotu B");
  status += s;


  /* cdotc */

  gsl_blas_raw_cdotc(4, vector_4_c, 1, vector_4_c, 1, c);
  s = ( frac_diff(c[0], 50.0) > eps_f || fabs(c[1]) > eps_f );
  gsl_test(s, "gsl_blas_raw_cdotc A");
  status += s;

  gsl_blas_raw_cdotc(2, vector_4_c, 3, vector_4_c, 1, c);
  s = ( frac_diff(c[0], 7.0) > eps_f || frac_diff(c[1], 8.0) > eps_f );
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


  /* isamax */

  bi = gsl_blas_raw_isamax(4, vector_4_f, 1);
  s = ( bi != 3 );
  gsl_test(s, "gsl_blas_raw_isamax A");
  status += s;

  bi = gsl_blas_raw_isamax(2, vector_4_f, 2);
  s = ( bi != 0 );
  gsl_test(s, "gsl_blas_raw_isamax B");
  status += s;


  /* idamax */

  bi = gsl_blas_raw_idamax(4, vector_4_d, 1);
  s = ( bi != 3 );
  gsl_test(s, "gsl_blas_raw_idamax A");
  status += s;

  bi = gsl_blas_raw_idamax(2, vector_4_d, 2);
  s = ( bi != 0 );
  gsl_test(s, "gsl_blas_raw_idamax B");
  status += s;


  /* icamax */

  bi = gsl_blas_raw_icamax(4, vector_4_c, 1);
  s = ( bi != 3 );
  gsl_test(s, "gsl_blas_raw_icamax A");
  status += s;

  bi = gsl_blas_raw_icamax(2, vector_4_c, 2);
  s = ( bi != 0 );
  gsl_test(s, "gsl_blas_raw_icamax B");
  status += s;


  /* izamax */

  bi = gsl_blas_raw_izamax(4, vector_4_z, 1);
  s = ( bi != 3 );
  gsl_test(s, "gsl_blas_raw_izamax A");
  status += s;

  bi = gsl_blas_raw_izamax(2, vector_4_z, 2);
  s = ( bi != 0 );
  gsl_test(s, "gsl_blas_raw_izamax B");
  status += s;


  /* sswap */
  memcpy(tmp_a_f, vector_4_f, 4*sizeof(float));
  memcpy(tmp_b_f, vector_4_zero_f, 4*sizeof(float));
  gsl_blas_raw_sswap(4, tmp_a_f, 1, tmp_b_f, 1);
  s = ( tmp_a_f[0] != 0.0 || tmp_b_f[0] != vector_4_f[0] );
  gsl_test(s, "gsl_blas_raw_sswap");
  status += s;


  /* dswap */
  memcpy(tmp_a_d, vector_4_d, 4*sizeof(double));
  memcpy(tmp_b_d, vector_4_zero_d, 4*sizeof(double));
  gsl_blas_raw_dswap(4, tmp_a_d, 1, tmp_b_d, 1);
  s = ( tmp_a_d[0] != 0.0 || tmp_b_d[0] != vector_4_d[0] );
  gsl_test(s, "gsl_blas_raw_dswap");
  status += s;


  /* cswap */
  memcpy(tmp_a_f, vector_4_c, 2*4*sizeof(float));
  memcpy(tmp_b_f, vector_4_zero_c, 2*4*sizeof(float));
  gsl_blas_raw_cswap(4, tmp_a_f, 1, tmp_b_f, 1);
  s = (   tmp_a_f[0] != 0.0 || tmp_a_f[1] != 0.0 || tmp_a_f[5] != 0.0
       || tmp_b_f[4] != vector_4_c[4]
       || tmp_b_f[5] != vector_4_c[5]
       );
  gsl_test(s, "gsl_blas_raw_cswap");
  status += s;


  /* zswap */
  memcpy(tmp_a_d, vector_4_z, 2*4*sizeof(double));
  memcpy(tmp_b_d, vector_4_zero_c, 2*4*sizeof(double));
  gsl_blas_raw_zswap(4, tmp_a_d, 1, tmp_b_d, 1);
  s = (   tmp_a_d[0] != 0.0 || tmp_a_d[1] != 0.0 || tmp_a_d[5] != 0.0
       || tmp_b_d[4] != vector_4_z[4]
       || tmp_b_d[5] != vector_4_z[5]
       );
  gsl_test(s, "gsl_blas_raw_zswap");
  status += s;


  /* scopy */
  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_a_f, 2);
  s = ( tmp_a_f[0] != vector_4_f[0] || tmp_a_f[2] != vector_4_f[1] );
  status += s;
  gsl_test(s, "gsl_blas_raw_scopy A");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_c, 2, tmp_a_f, 2);
  s = ( tmp_a_f[0] != vector_4_f[0] || tmp_a_f[2] != vector_4_f[1] );
  status += s;
  gsl_test(s, "gsl_blas_raw_scopy B");
  status += s;


  /* dcopy */
  gsl_blas_raw_dcopy(4, vector_4_d, 1, tmp_a_d, 2);
  s = ( tmp_a_d[0] != vector_4_d[0] || tmp_a_d[2] != vector_4_d[1] );
  status += s;
  gsl_test(s, "gsl_blas_raw_dcopy");
  status += s;


  /* ccopy */
  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_a_f, 2);
  s = (   tmp_a_f[0] != vector_4_c[0] || tmp_a_f[1] != vector_4_c[1]
       || tmp_a_f[4] != vector_4_c[2] || tmp_a_f[5] != vector_4_c[3]
       );
  status += s;
  gsl_test(s, "gsl_blas_raw_ccopy");
  status += s;


  /* zcopy */
  gsl_blas_raw_zcopy(4, vector_4_z, 1, tmp_a_d, 2);
  s = (   tmp_a_d[0] != vector_4_z[0] || tmp_a_d[1] != vector_4_z[1]
       || tmp_a_d[4] != vector_4_z[2] || tmp_a_d[5] != vector_4_z[3]
       );
  status += s;
  gsl_test(s, "gsl_blas_raw_zcopy");
  status += s;


  /* saxpy */

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_a_f, 1);
  gsl_blas_raw_saxpy(4, 2.0, vector_4_f, 1, tmp_a_f, 1);
  s = (   tmp_a_f[0] != 3.0 * vector_4_f[0]
       || tmp_a_f[1] != 3.0 * vector_4_f[1]
       || tmp_a_f[2] != 3.0 * vector_4_f[2]
       );
  gsl_test(s, "gsl_blas_raw_saxpy A");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_a_f, 1);
  gsl_blas_raw_saxpy(4, 2.0, vector_4_c, 2, tmp_a_f, 1);
  s = (   tmp_a_f[0] != 3.0 * vector_4_f[0]
       || tmp_a_f[1] != 3.0 * vector_4_f[1]
       || tmp_a_f[2] != 3.0 * vector_4_f[2]
       );
  gsl_test(s, "gsl_blas_raw_saxpy B");
  status += s;


  /* daxpy */

  gsl_blas_raw_dcopy(4, vector_4_d, 1, tmp_a_d, 1);
  gsl_blas_raw_daxpy(4, 2.0, vector_4_d, 1, tmp_a_d, 1);
  s = (   tmp_a_d[0] != 3.0 * vector_4_d[0]
       || tmp_a_d[1] != 3.0 * vector_4_d[1]
       || tmp_a_d[2] != 3.0 * vector_4_d[2]
       );
  gsl_test(s, "gsl_blas_raw_daxpy");
  status += s;


  /* caxpy */

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_a_f, 1);
  gsl_blas_raw_caxpy(4, alpha_c, vector_4_c, 1, tmp_a_f, 1);
  s = (   tmp_a_f[0] != -5.0
       || tmp_a_f[1] !=  0.0
       || tmp_a_f[2] != -3.0
       || tmp_a_f[3] !=  1.0
       );
  gsl_test(s, "gsl_blas_raw_caxpy");
  status += s;


  /* zaxpy */

  gsl_blas_raw_zcopy(4, vector_4_z, 1, tmp_a_d, 1);
  gsl_blas_raw_zaxpy(4, alpha_z, vector_4_z, 1, tmp_a_d, 1);
  s = (   tmp_a_d[0] != -5.0
       || tmp_a_d[1] !=  0.0
       || tmp_a_d[2] != -3.0
       || tmp_a_d[3] !=  1.0
       );
  gsl_test(s, "gsl_blas_raw_zaxpy");
  status += s;


  /* rotg */

  tmp_a_f[0] = 1.0;
  tmp_a_f[1] = 2.0;
  tmp_a_f[2] = 3.0;
  tmp_a_f[3] = 4.0;
  gsl_blas_raw_srotg (&(tmp_a_f[0]), &(tmp_a_f[1]), &(tmp_a_f[2]), &(tmp_a_f[3]));
  s = ( frac_diff(tmp_a_f[0], 2.23606801) > eps_f ||
        frac_diff(tmp_a_f[1], 2.23606801) > eps_f ||
	frac_diff(tmp_a_f[2], 0.44721359) > eps_f ||
	frac_diff(tmp_a_f[3], 0.89442718) > eps_f
       );
  gsl_test(s, "gsl_blas_raw_srotg");
  status += s;

  tmp_a_d[0] = 1.0;
  tmp_a_d[1] = 2.0;
  tmp_a_d[2] = 3.0;
  tmp_a_d[3] = 4.0;
  gsl_blas_raw_drotg (&(tmp_a_d[0]), &(tmp_a_d[1]), &(tmp_a_d[2]), &(tmp_a_d[3]));
  s = ( frac_diff(tmp_a_d[0], 2.2360679774997894) > eps ||
        frac_diff(tmp_a_d[1], 2.2360679774997894) > eps ||
	frac_diff(tmp_a_d[2], 0.44721359549995798) > eps ||
	frac_diff(tmp_a_d[3], 0.89442719099991597) > eps
       );
  gsl_test(s, "gsl_blas_raw_drotg");
  status += s;


  /* rot */

  tmp_a_f[0] = 1.0;
  tmp_a_f[1] = 2.0;
  tmp_a_f[2] = 3.0;
  tmp_a_f[3] = 4.0;
  tmp_b_f[0] = 5.0;
  tmp_b_f[1] = 6.0;
  tmp_b_f[2] = 7.0;
  tmp_b_f[3] = 8.0;
  gsl_blas_raw_srot (4, tmp_a_f, 1, tmp_b_f, 1, 0.5, sqrt(1.0 - 0.25));
  s = ( frac_diff(tmp_a_f[0], 4.83012676) > eps_f ||
        frac_diff(tmp_a_f[1], 6.19615221) > eps_f ||
	frac_diff(tmp_a_f[2], 7.56217766) > eps_f ||
	frac_diff(tmp_a_f[3], 8.92820358) > eps_f ||
	frac_diff(tmp_b_f[0], 1.63397455) > eps_f ||
        frac_diff(tmp_b_f[1], 1.26794922) > eps_f ||
	frac_diff(tmp_b_f[2], 0.901923835) > eps_f ||
	frac_diff(tmp_b_f[3], 0.535898447) > eps_f
        );
  gsl_test(s, "gsl_blas_raw_srot A");
  status += s;

  tmp_a_f[0] = 1.0;
  tmp_a_f[1] = 2.0;
  tmp_a_f[2] = 3.0;
  tmp_a_f[3] = 4.0;
  tmp_b_f[0] = 5.0;
  tmp_b_f[1] = 6.0;
  tmp_b_f[2] = 7.0;
  tmp_b_f[3] = 8.0;
  gsl_blas_raw_srot (2, tmp_a_f, 2, tmp_b_f, 2, 0.5, sqrt(1.0 - 0.25));
  s = ( frac_diff(tmp_a_f[0], 4.83012676) > eps_f ||
	frac_diff(tmp_a_f[2], 7.56217766) > eps_f ||
	frac_diff(tmp_b_f[0], 1.63397455) > eps_f ||
	frac_diff(tmp_b_f[2], 0.901923835) > eps_f
	);
  gsl_test(s, "gsl_blas_raw_srot B");
  status += s;

  tmp_a_d[0] = 1.0;
  tmp_a_d[1] = 2.0;
  tmp_a_d[2] = 3.0;
  tmp_a_d[3] = 4.0;
  tmp_b_d[0] = 5.0;
  tmp_b_d[1] = 6.0;
  tmp_b_d[2] = 7.0;
  tmp_b_d[3] = 8.0;
  gsl_blas_raw_drot (4, tmp_a_d, 1, tmp_b_d, 1, 0.5, sqrt(1.0 - 0.25));
  s = ( frac_diff(tmp_a_d[0], 4.8301270189221928) > eps ||
        frac_diff(tmp_a_d[1], 6.1961524227066320) > eps ||
	frac_diff(tmp_a_d[2], 7.5621778264910704) > eps ||
	frac_diff(tmp_a_d[3], 8.9282032302755088) > eps ||
	frac_diff(tmp_b_d[0], 1.6339745962155614) > eps ||
        frac_diff(tmp_b_d[1], 1.2679491924311228) > eps ||
	frac_diff(tmp_b_d[2], 0.90192378864668421) > eps ||
	frac_diff(tmp_b_d[3], 0.53589838486224561) > eps
        );
  gsl_test(s, "gsl_blas_raw_drot");
  status += s;


/* FIXME: rotmg, rotm */


  /* sscal */

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_a_f, 2);
  gsl_blas_raw_sscal(4, 2.0, tmp_a_f, 2);
  s = (   tmp_a_f[0] != 2.0 * vector_4_f[0]
       || tmp_a_f[2] != 2.0 * vector_4_f[1]
       || tmp_a_f[4] != 2.0 * vector_4_f[2]
       );
  gsl_test(s, "gsl_blas_raw_sscal");
  status += s;


  /* dscal */

  gsl_blas_raw_dcopy(4, vector_4_d, 1, tmp_a_d, 2);
  gsl_blas_raw_dscal(4, 2.0, tmp_a_d, 2);
  s = (   tmp_a_d[0] != 2.0 * vector_4_d[0]
       || tmp_a_d[2] != 2.0 * vector_4_d[1]
       || tmp_a_d[4] != 2.0 * vector_4_d[2]
       );
  gsl_test(s, "gsl_blas_raw_dscal");
  status += s;


  /* cscal */

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_a_f, 2);
  gsl_blas_raw_cscal(4, alpha_c, tmp_a_f, 2);
  s = (   tmp_a_f[0] != -3.0
       || tmp_a_f[1] != -1.0
       || tmp_a_f[4] != -2.0
       || tmp_a_f[5] !=  0.0
       );
  gsl_test(s, "gsl_blas_raw_cscal");
  status += s;


  /* zscal */

  gsl_blas_raw_zcopy(4, vector_4_z, 1, tmp_a_d, 2);
  gsl_blas_raw_zscal(4, alpha_z, tmp_a_d, 2);
  s = (   tmp_a_d[0] != -3.0
       || tmp_a_d[1] != -1.0
       || tmp_a_d[4] != -2.0
       || tmp_a_d[5] !=  0.0
       );
  gsl_test(s, "gsl_blas_raw_zscal");
  status += s;


  /* csscal */

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_a_f, 2);
  gsl_blas_raw_csscal(4, 2.0, tmp_a_f, 2);
  s = (   tmp_a_f[0] != -4.0
       || tmp_a_f[1] !=  2.0
       || tmp_a_f[4] != -2.0
       || tmp_a_f[5] !=  2.0
       );
  gsl_test(s, "gsl_blas_raw_csscal");
  status += s;


  /* zdscal */

  gsl_blas_raw_zcopy(4, vector_4_z, 1, tmp_a_d, 2);
  gsl_blas_raw_zdscal(4, 2.0, tmp_a_d, 2);
  s = (   tmp_a_d[0] != -4.0
       || tmp_a_d[1] !=  2.0
       || tmp_a_d[4] != -2.0
       || tmp_a_d[5] !=  2.0
       );
  gsl_test(s, "gsl_blas_raw_zdscal");
  status += s;


  return status;
}


/**** level 2 tests ****/

int test_gemv(void)
{
  int status = 0;
  int s;

  float  tmp_f[32];
  double tmp_d[32];

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_sgemv(CblasNoTrans, 4, 4, 2.0, matrix_gen_4_f, 4, vector_4_f, 1, 3.0, tmp_f, 1);
  s = ( tmp_f[0] != -4.0 || tmp_f[1] != -5.0 || tmp_f[2] != -1.0 || tmp_f[3] != 7.0 );
  gsl_test(s, "gsl_blas_raw_sgemv A");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_sgemv(CblasTrans, 4, 4, 2.0, matrix_gen_4_f, 4, vector_4_f, 1, 3.0, tmp_f, 1);
  s = ( tmp_f[0] != -5.0 || tmp_f[1] != -15.0 || tmp_f[2] != 1.0 || tmp_f[3] != 3.0 );
  gsl_test(s, "gsl_blas_raw_sgemv B");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 2);
  gsl_blas_raw_sgemv(CblasNoTrans, 4, 4, 2.0, matrix_gen_4_f, 4, vector_4_f, 1, 3.0, tmp_f, 2);
  s = ( tmp_f[0] != -4.0 || tmp_f[2] != -5.0 || tmp_f[4] != -1.0 || tmp_f[6] != 7.0 );
  gsl_test(s, "gsl_blas_raw_sgemv C");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 2);
  gsl_blas_raw_sgemv(CblasNoTrans, 4, 4, 2.0, matrix_gen_4_f, 4, vector_4_c, 2, 3.0, tmp_f, 2);
  s = ( tmp_f[0] != -4.0 || tmp_f[2] != -5.0 || tmp_f[4] != -1.0 || tmp_f[6] != 7.0 );
  gsl_test(s, "gsl_blas_raw_sgemv D");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 2);
  gsl_blas_raw_sgemv(CblasTrans, 4, 4, 2.0, matrix_gen_4_f, 4, vector_4_f, 1, 3.0, tmp_f, 2);
  s = ( tmp_f[0] != -5.0 || tmp_f[2] != -15.0 || tmp_f[4] != 1.0 || tmp_f[6] != 3.0 );
  gsl_test(s, "gsl_blas_raw_sgemv E");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 2);
  gsl_blas_raw_sgemv(CblasTrans, 4, 4, 2.0, matrix_gen_4_f, 4, vector_4_c, 2, 3.0, tmp_f, 2);
  s = ( tmp_f[0] != -5.0 || tmp_f[2] != -15.0 || tmp_f[4] != 1.0 || tmp_f[6] != 3.0 );
  gsl_test(s, "gsl_blas_raw_sgemv F");
  status += s;


  gsl_blas_raw_dcopy(4, vector_4_d, 1, tmp_d, 1);
  gsl_blas_raw_dgemv(CblasNoTrans, 4, 4, 2.0, matrix_gen_4_d, 4, vector_4_d, 1, 3.0, tmp_d, 1);
  s = ( tmp_d[0] != -4.0 || tmp_d[1] != -5.0 || tmp_d[2] != -1.0 || tmp_d[3] != 7.0 );
  gsl_test(s, "gsl_blas_raw_dgemv A");
  status += s;

  gsl_blas_raw_dcopy(4, vector_4_d, 1, tmp_d, 1);
  gsl_blas_raw_dgemv(CblasNoTrans, 4, 4, 2.0, matrix_gen_4_d, 4, vector_4_z, 2, 3.0, tmp_d, 1);
  s = ( tmp_d[0] != -4.0 || tmp_d[1] != -5.0 || tmp_d[2] != -1.0 || tmp_d[3] != 7.0 );
  gsl_test(s, "gsl_blas_raw_dgemv B");
  status += s;


  return status;
}

int test_gbmv(void)
{
  int status = 0;
  int s;

  float  tmp_f[32];
  double tmp_d[32];


  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_sgbmv(CblasNoTrans, 4, 4, 1, 2, 2.0, matrix_gen_4_f, 4, vector_4_f, 1, 3.0, tmp_f, 1);
  s = ( tmp_f[0] != -10.0 || tmp_f[1] != -5.0 || tmp_f[2] != 7.0 || tmp_f[3] != 9.0 );
  gsl_test(s, "gsl_blas_raw_sgbmv A");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_sgbmv(CblasTrans, 4, 4, 2, 1, 2.0, matrix_gen_4_f, 4, vector_4_f, 1, 3.0, tmp_f, 1);
  s = ( tmp_f[0] != -11.0 || tmp_f[1] != -9.0 || tmp_f[2] != 1.0 || tmp_f[3] != 7.0 );
  gsl_test(s, "gsl_blas_raw_sgbmv B");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 2);
  gsl_blas_raw_sgbmv(CblasNoTrans, 4, 4, 1, 2, 2.0, matrix_gen_4_f, 4, vector_4_f, 1, 3.0, tmp_f, 2);
  s = ( tmp_f[0] != -10.0 || tmp_f[2] != -5.0 || tmp_f[4] != 7.0 || tmp_f[6] != 9.0 );
  gsl_test(s, "gsl_blas_raw_sgbmv C");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 2);
  gsl_blas_raw_sgbmv(CblasNoTrans, 4, 4, 1, 2, 2.0, matrix_gen_4_f, 4, vector_4_c, 2, 3.0, tmp_f, 2);
  s = ( tmp_f[0] != -10.0 || tmp_f[2] != -5.0 || tmp_f[4] != 7.0 || tmp_f[6] != 9.0 );
  gsl_test(s, "gsl_blas_raw_sgbmv D");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 2);
  gsl_blas_raw_sgbmv(CblasTrans, 4, 4, 2, 1, 2.0, matrix_gen_4_f, 4, vector_4_f, 1, 3.0, tmp_f, 2);
  s = ( tmp_f[0] != -11.0 || tmp_f[2] != -9.0 || tmp_f[4] != 1.0 || tmp_f[6] != 7.0 );
  gsl_test(s, "gsl_blas_raw_sgbmv E");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 2);
  gsl_blas_raw_sgbmv(CblasTrans, 4, 4, 2, 1, 2.0, matrix_gen_4_f, 4, vector_4_c, 2, 3.0, tmp_f, 2);
  s = ( tmp_f[0] != -11.0 || tmp_f[2] != -9.0 || tmp_f[4] != 1.0 || tmp_f[6] != 7.0 );
  gsl_test(s, "gsl_blas_raw_sgbmv E");
  status += s;


  gsl_blas_raw_dcopy(4, vector_4_d, 1, tmp_d, 1);
  gsl_blas_raw_dgbmv(CblasNoTrans, 4, 4, 1, 2, 2.0, matrix_gen_4_d, 4, vector_4_d, 1, 3.0, tmp_d, 1);
  s = ( tmp_d[0] != -10.0 || tmp_d[1] != -5.0 || tmp_d[2] != 7.0 || tmp_d[3] != 9.0 );
  gsl_test(s, "gsl_blas_raw_dgbmv A");
  status += s;

  gsl_blas_raw_dcopy(4, vector_4_d, 1, tmp_d, 1);
  gsl_blas_raw_dgbmv(CblasNoTrans, 4, 4, 1, 2, 2.0, matrix_gen_4_d, 4, vector_4_z, 2, 3.0, tmp_d, 1);
  s = ( tmp_d[0] != -10.0 || tmp_d[1] != -5.0 || tmp_d[2] != 7.0 || tmp_d[3] != 9.0 );
  gsl_test(s, "gsl_blas_raw_dgbmv B");
  status += s;


  return status;
}

int test_trmv(void)
{
  int status = 0;
  int s;

  float  tmp_f[32];
  double tmp_d[32];

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_strmv(CblasUpper, CblasNoTrans, CblasNonUnit, 4, matrix_gen_4_f, 4, tmp_f, 1);
  s = ( tmp_f[0] != 1.0 || tmp_f[1] != 0.0 || tmp_f[2] != 1.5 || tmp_f[3] != 0.0 );
  gsl_test(s, "gsl_blas_raw_strmv A");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_strmv(CblasUpper, CblasNoTrans, CblasUnit, 4, matrix_gen_4_f, 4, tmp_f, 1);
  s = ( tmp_f[0] != 1.0 || tmp_f[1] != 2.0 || tmp_f[2] != 1.5 || tmp_f[3] != 3.0 );
  gsl_test(s, "gsl_blas_raw_strmv B");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_strmv(CblasUpper, CblasTrans, CblasNonUnit, 4, matrix_gen_4_f, 4, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] != -3.0 || tmp_f[2] != -1.0 || tmp_f[3] != -3.0 );
  gsl_test(s, "gsl_blas_raw_strmv C");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_strmv(CblasUpper, CblasTrans, CblasUnit, 4, matrix_gen_4_f, 4, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] != -1.0 || tmp_f[2] != -1.0 || tmp_f[3] != 0.0 );
  gsl_test(s, "gsl_blas_raw_strmv D");
  status += s;

  gsl_blas_raw_scopy(8, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_strmv(CblasUpper, CblasTrans, CblasUnit, 4, matrix_gen_4_f, 4, tmp_f, 2);
  s = ( tmp_f[0] != -2.0 || tmp_f[2] != -1.0 || tmp_f[4] != -1.0 || tmp_f[6] != 0.0 );
  gsl_test(s, "gsl_blas_raw_strmv E");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_strmv(CblasLower, CblasNoTrans, CblasNonUnit, 4, matrix_gen_4_f, 4, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] != -4.0 || tmp_f[2] != -2.0 || tmp_f[3] != -1.0 );
  gsl_test(s, "gsl_blas_raw_strmv F");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_strmv(CblasLower, CblasTrans, CblasNonUnit, 4, matrix_gen_4_f, 4, tmp_f, 1);
  s = ( tmp_f[0] != 0.5 || tmp_f[1] != -6.0 || tmp_f[2] != 1.5 || tmp_f[3] != 0.0 );
  gsl_test(s, "gsl_blas_raw_strmv G");
  status += s;

  gsl_blas_raw_scopy(8, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_strmv(CblasLower, CblasTrans, CblasNonUnit, 4, matrix_gen_4_f, 4, tmp_f, 2);
  s = ( tmp_f[0] != 0.5 || tmp_f[2] != -6.0 || tmp_f[4] != 1.5 || tmp_f[6] != 0.0 );
  gsl_test(s, "gsl_blas_raw_strmv H");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_strmv(CblasLower, CblasTrans, CblasUnit, 4, matrix_gen_4_f, 4, tmp_f, 1);
  s = ( tmp_f[0] != 0.5 || tmp_f[1] != -4.0 || tmp_f[2] != 1.5 || tmp_f[3] != 3.0 );
  gsl_test(s, "gsl_blas_raw_strmv I");
  status += s;


  gsl_blas_raw_dcopy(4, vector_4_d, 1, tmp_d, 1);
  gsl_blas_raw_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit, 4, matrix_gen_4_d, 4, tmp_d, 1);
  s = ( tmp_d[0] != 1.0 || tmp_d[1] != 0.0 || tmp_d[2] != 1.5 || tmp_d[3] != 0.0 );
  gsl_test(s, "gsl_blas_raw_dtrmv A");
  status += s;

  gsl_blas_raw_dcopy(8, vector_4_z, 1, tmp_d, 1);
  gsl_blas_raw_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit, 4, matrix_gen_4_d, 4, tmp_d, 2);
  s = ( tmp_d[0] != 1.0 || tmp_d[2] != 0.0 || tmp_d[4] != 1.5 || tmp_d[6] != 0.0 );
  gsl_test(s, "gsl_blas_raw_dtrmv B");
  status += s;


  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctrmv(CblasUpper, CblasNoTrans, CblasNonUnit, 4, matrix_gen_4_c, 4, tmp_f, 1);
  s = ( tmp_f[0] != -7.0 || tmp_f[1] != 3.0 || tmp_f[2] != -17.0 || tmp_f[3] != 23.0 ||
        tmp_f[4] != -27.0/2.0 || tmp_f[5] != 17.0/2.0 || tmp_f[6] != -5.0 || tmp_f[7] != 3.0
       );
  gsl_test(s, "gsl_blas_raw_ctrmv A");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctrmv(CblasUpper, CblasNoTrans, CblasUnit, 4, matrix_gen_4_c, 4, tmp_f, 1);
  s = ( tmp_f[0] != -7.0 || tmp_f[1] != 3.0 || tmp_f[2] != -14.0 || tmp_f[3] != 22.0 ||
        tmp_f[4] != -27.0/2.0 || tmp_f[5] != 29.0/2.0 || tmp_f[6] != 3.0 || tmp_f[7] != 5.0
       );
  gsl_test(s, "gsl_blas_raw_ctrmv B");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctrmv(CblasUpper, CblasTrans, CblasNonUnit, 4, matrix_gen_4_c, 4, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] != 1.0 || tmp_f[2] != -4.0 || tmp_f[3] != 2.0 ||
        tmp_f[4] != -3.0 || tmp_f[5] != -3.0 || tmp_f[6] != -19.0 || tmp_f[7] != 1.5
       );
  gsl_test(s, "gsl_blas_raw_ctrmv C");
  status += s;


  gsl_blas_raw_zcopy(4, vector_4_z, 1, tmp_d, 1);
  gsl_blas_raw_ztrmv(CblasUpper, CblasNoTrans, CblasNonUnit, 4, matrix_gen_4_z, 4, tmp_d, 1);
  s = ( tmp_d[0] != -7.0 || tmp_d[1] != 3.0 || tmp_d[2] != -17.0 || tmp_d[3] != 23.0 ||
        tmp_d[4] != -27.0/2.0 || tmp_d[5] != 17.0/2.0 || tmp_d[6] != -5.0 || tmp_d[7] != 3.0
       );
  gsl_test(s, "gsl_blas_raw_ztrmv A");
  status += s;


  return status;
}

int test_tbmv(void)
{
  int status = 0;
  int s;

  float  tmp_f[32];
  double tmp_d[32];

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_stbmv(CblasUpper, CblasNoTrans, CblasNonUnit, 4, 2, matrix_gen_4_f, 4, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] != 0.0 || tmp_f[2] != 1.5 || tmp_f[3] != 0.0 );
  gsl_test(s, "gsl_blas_raw_stbmv A");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_stbmv(CblasUpper, CblasNoTrans, CblasUnit, 4, 2, matrix_gen_4_f, 4, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] != 2.0 || tmp_f[2] != 1.5 || tmp_f[3] != 3.0 );
  gsl_test(s, "gsl_blas_raw_stbmv B");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_stbmv(CblasUpper, CblasTrans, CblasNonUnit, 4, 2, matrix_gen_4_f, 4, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] != -3.0 || tmp_f[2] != -1.0 || tmp_f[3] != -1.0 );
  gsl_test(s, "gsl_blas_raw_stbmv C");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_stbmv(CblasUpper, CblasTrans, CblasUnit, 4, 2, matrix_gen_4_f, 4, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] != -1.0 || tmp_f[2] != -1.0 || tmp_f[3] != 2.0 );
  gsl_test(s, "gsl_blas_raw_stbmv D");
  status += s;

  gsl_blas_raw_scopy(8, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_stbmv(CblasUpper, CblasTrans, CblasUnit, 4, 2, matrix_gen_4_f, 4, tmp_f, 2);
  s = ( tmp_f[0] != -2.0 || tmp_f[2] != -1.0 || tmp_f[4] != -1.0 || tmp_f[6] != 2.0 );
  gsl_test(s, "gsl_blas_raw_stbmv E");
  status += s;


  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctbmv(CblasUpper, CblasNoTrans, CblasNonUnit, 4, 2, matrix_gen_4_c, 4, tmp_f, 1);
  s = ( tmp_f[0] != -5.0 || tmp_f[1] != -5.0 || tmp_f[2] != -17.0 || tmp_f[3] != 23.0 ||
        tmp_f[4] != -27.0/2.0 || tmp_f[5] != 17.0/2.0 || tmp_f[6] != -5.0 || tmp_f[7] != 3.0
       );
  gsl_test(s, "gsl_blas_raw_ctbmv A");
  status += s;


  gsl_blas_raw_zcopy(4, vector_4_z, 1, tmp_d, 1);
  gsl_blas_raw_ztbmv(CblasUpper, CblasNoTrans, CblasNonUnit, 4, 2, matrix_gen_4_z, 4, tmp_d, 1);
  s = ( tmp_d[0] != -5.0 || tmp_d[1] != -5.0 || tmp_d[2] != -17.0 || tmp_d[3] != 23.0 ||
        tmp_d[4] != -27.0/2.0 || tmp_d[5] != 17.0/2.0 || tmp_d[6] != -5.0 || tmp_d[7] != 3.0
       );
  gsl_test(s, "gsl_blas_raw_ztbmv A");
  status += s;


  return status;
}

int test_tpmv(void)
{
  int status = 0;
  int s;

  float  tmp_f[32];
  double tmp_d[32];

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_stpmv(CblasUpper, CblasNoTrans, CblasNonUnit, 4, matrix_gen_4_f, tmp_f, 1);
  s = ( tmp_f[0] != 1.0 || tmp_f[1] != 29.0/2.0 || tmp_f[2] != 6.0 || tmp_f[3] != -6.0 );
  gsl_test(s, "gsl_blas_raw_stpmv A");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_stpmv(CblasUpper, CblasNoTrans, CblasUnit, 4, matrix_gen_4_f, tmp_f, 1);
  s = ( tmp_f[0] != 1.0 || tmp_f[1] != 14.0 || tmp_f[2] != 6.0 || tmp_f[3] != 3.0 );
  gsl_test(s, "gsl_blas_raw_stpmv B");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_stpmv(CblasLower, CblasNoTrans, CblasNonUnit, 4, matrix_gen_4_f, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] != 2.0 || tmp_f[2] != -5.0/2.0 || tmp_f[3] != -17.0 );
  gsl_test(s, "gsl_blas_raw_stpmv C");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_stpmv(CblasLower, CblasNoTrans, CblasUnit, 4, matrix_gen_4_f, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] != -1.0 || tmp_f[2] != -5.0/2.0 || tmp_f[3] != -8.0 );
  gsl_test(s, "gsl_blas_raw_stpmv D");
  status += s;

  gsl_blas_raw_scopy(8, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_stpmv(CblasLower, CblasNoTrans, CblasUnit, 4, matrix_gen_4_f, tmp_f, 2);
  s = ( tmp_f[0] != -2.0 || tmp_f[2] != -1.0 || tmp_f[4] != -5.0/2.0 || tmp_f[6] != -8.0 );
  gsl_test(s, "gsl_blas_raw_stpmv E");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_stpmv(CblasUpper, CblasTrans, CblasNonUnit, 4, matrix_gen_4_f, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] != -0.5 || tmp_f[2] != 1.0 || tmp_f[3] != -13.0 );
  gsl_test(s, "gsl_blas_raw_stpmv F");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_stpmv(CblasUpper, CblasTrans, CblasUnit, 4, matrix_gen_4_f, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] != -1.0 || tmp_f[2] != 1.0 || tmp_f[3] != -4.0 );
  gsl_test(s, "gsl_blas_raw_stpmv G");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_stpmv(CblasLower, CblasTrans, CblasNonUnit, 4, matrix_gen_4_f, tmp_f, 1);
  s = ( tmp_f[0] != 13.0 || tmp_f[1] != 5.0 || tmp_f[2] != 6.0 || tmp_f[3] != -6.0 );
  gsl_test(s, "gsl_blas_raw_stpmv H");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_stpmv(CblasLower, CblasTrans, CblasUnit, 4, matrix_gen_4_f, tmp_f, 1);
  s = ( tmp_f[0] != 13.0 || tmp_f[1] != 2.0 || tmp_f[2] != 6.0 || tmp_f[3] != 3.0 );
  gsl_test(s, "gsl_blas_raw_stpmv I");
  status += s;

  gsl_blas_raw_scopy(8, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_stpmv(CblasLower, CblasTrans, CblasUnit, 4, matrix_gen_4_f, tmp_f, 2);
  s = ( tmp_f[0] != 13.0 || tmp_f[2] != 2.0 || tmp_f[4] != 6.0 || tmp_f[6] != 3.0 );
  gsl_test(s, "gsl_blas_raw_stpmv J");
  status += s;


  gsl_blas_raw_dcopy(4, vector_4_d, 1, tmp_d, 1);
  gsl_blas_raw_dtpmv(CblasUpper, CblasNoTrans, CblasNonUnit, 4, matrix_gen_4_d, tmp_d, 1);
  s = ( tmp_d[0] != 1.0 || tmp_d[1] != 29.0/2.0 || tmp_d[2] != 6.0 || tmp_d[3] != -6.0 );
  gsl_test(s, "gsl_blas_raw_dtpmv A");
  status += s;


  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctpmv(CblasUpper, CblasNoTrans, CblasNonUnit, 4, matrix_gen_4_c, tmp_f, 1);
  s = ( tmp_f[0] != -7.0 || tmp_f[1] != 3.0 || tmp_f[2] != 13.0/2.0 || tmp_f[3] != 75.0/2.0 ||
        tmp_f[4] != -5.0 || tmp_f[5] != 13.0 || tmp_f[6] != -6.0 || tmp_f[7] != -10.0
       );
  gsl_test(s, "gsl_blas_raw_ctpmv A");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctpmv(CblasUpper, CblasNoTrans, CblasUnit, 4, matrix_gen_4_c, tmp_f, 1);
  s = ( tmp_f[0] != -7.0 || tmp_f[1] != 3.0 || tmp_f[2] != 6.0 || tmp_f[3] != 38.0 ||
        tmp_f[4] != 1.0 || tmp_f[5] != 16.0 || tmp_f[6] != 3.0 || tmp_f[7] != 5.0
       );
  gsl_test(s, "gsl_blas_raw_ctpmv B");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctpmv(CblasUpper, CblasTrans, CblasNonUnit, 4, matrix_gen_4_c, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] != 1.0 || tmp_f[2] != -0.5 || tmp_f[3] != 0.5 ||
        tmp_f[4] != -7.0 || tmp_f[5] != -2.0 || tmp_f[6] != -18.0 || tmp_f[7] != -1.0
       );
  gsl_test(s, "gsl_blas_raw_ctpmv C");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctpmv(CblasUpper, CblasTrans, CblasUnit, 4, matrix_gen_4_c, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] != 1.0 || tmp_f[2] != -1.0 || tmp_f[3] != 1.0 ||
        tmp_f[4] != -1.0 || tmp_f[5] != 1.0 || tmp_f[6] != -9.0 || tmp_f[7] != 14.0
       );
  gsl_test(s, "gsl_blas_raw_ctpmv D");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctpmv(CblasLower, CblasNoTrans, CblasUnit, 4, matrix_gen_4_c, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] != 1.0 || tmp_f[2] != -1.0 || tmp_f[3] != 1.0 ||
        tmp_f[4] != -7.0/2.0 || tmp_f[5] != 5.0/2.0 || tmp_f[6] != -13.0 || tmp_f[7] != 12.0
       );
  gsl_test(s, "gsl_blas_raw_ctpmv E");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctpmv(CblasLower, CblasTrans, CblasUnit, 4, matrix_gen_4_c, tmp_f, 1);
  s = ( tmp_f[0] != 5.0 || tmp_f[1] != 32.0 || tmp_f[2] != -11.0 || tmp_f[3] != 17.0/2.0 ||
        tmp_f[4] != 1.0 || tmp_f[5] != 16.0 || tmp_f[6] != 3.0 || tmp_f[7] != 5.0
       );
  gsl_test(s, "gsl_blas_raw_ctpmv F");
  status += s;


  gsl_blas_raw_zcopy(4, vector_4_z, 1, tmp_d, 1);
  gsl_blas_raw_ztpmv(CblasUpper, CblasNoTrans, CblasNonUnit, 4, matrix_gen_4_z, tmp_d, 1);
  s = ( tmp_d[0] != -7.0 || tmp_d[1] != 3.0 || tmp_d[2] != 13.0/2.0 || tmp_d[3] != 75.0/2.0 ||
        tmp_d[4] != -5.0 || tmp_d[5] != 13.0 || tmp_d[6] != -6.0 || tmp_d[7] != -10.0
       );
  gsl_test(s, "gsl_blas_raw_ztpmv A");
  status += s;


  return status;
}


int test_trsv(void)
{
  int status = 0;
  int s;

  float  tmp_f[32];
  double tmp_d[32];

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_strsv(CblasUpper, CblasNoTrans, CblasNonUnit, 4, matrix_t_4_f, 4, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || frac_diff(tmp_f[1], -23.0/6.0) > eps_f ||
        tmp_f[2] !=  1.5 || tmp_f[3] != 3.0
       );
  gsl_test(s, "gsl_blas_raw_strsv A");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 2);
  gsl_blas_raw_strsv(CblasUpper, CblasNoTrans, CblasNonUnit, 4, matrix_t_4_f, 4, tmp_f, 2);
  s = ( tmp_f[0] != -2.0 || frac_diff(tmp_f[2], -23.0/6.0) > eps_f ||
        tmp_f[4] !=  1.5 || tmp_f[6] != 3.0
       );
  gsl_test(s, "gsl_blas_raw_strsv B");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_strsv(CblasUpper, CblasNoTrans, CblasUnit, 4, matrix_t_4_f, 4, tmp_f, 1);
  s = ( tmp_f[0] != -8.0 || tmp_f[1] != 3.5 ||
        tmp_f[2] != -1.5 || tmp_f[3] != 3.0
       );
  gsl_test(s, "gsl_blas_raw_strsv C");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_strsv(CblasLower, CblasNoTrans, CblasNonUnit, 4, matrix_t_4_f, 4, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] != 0.0 ||
        tmp_f[2] != -4.0 || tmp_f[3] != 7.0
       );
  gsl_test(s, "gsl_blas_raw_strsv D");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 2);
  gsl_blas_raw_strsv(CblasLower, CblasNoTrans, CblasNonUnit, 4, matrix_t_4_f, 4, tmp_f, 2);
  s = ( tmp_f[0] != -2.0 || tmp_f[2] != 0.0 ||
        tmp_f[4] != -4.0 || tmp_f[6] != 7.0
       );
  gsl_test(s, "gsl_blas_raw_strsv E");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_strsv(CblasLower, CblasNoTrans, CblasUnit, 4, matrix_t_4_f, 4, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] != 0.0 ||
        tmp_f[2] !=  4.0 || tmp_f[3] != 3.0
       );
  gsl_test(s, "gsl_blas_raw_strsv F");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_strsv(CblasUpper, CblasTrans, CblasNonUnit, 4, matrix_t_4_f, 4, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 ||
        frac_diff(tmp_f[1], -1.0/3.0) > eps_f ||
        frac_diff(tmp_f[2],  7.0/3.0) > eps_f ||
	frac_diff(tmp_f[3], 25.0/6.0) > eps_f
       );
  gsl_test(s, "gsl_blas_raw_strsv G");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_strsv(CblasUpper, CblasTrans, CblasUnit, 4, matrix_t_4_f, 4, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 ||
        tmp_f[1] != -1.0 ||
	tmp_f[2] !=  1.0 ||
	tmp_f[3] !=  5.5
       );
  gsl_test(s, "gsl_blas_raw_strsv H");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 2);
  gsl_blas_raw_strsv(CblasUpper, CblasTrans, CblasUnit, 4, matrix_t_4_f, 4, tmp_f, 2);
  s = ( tmp_f[0] != -2.0 ||
        tmp_f[2] != -1.0 ||
	tmp_f[4] !=  1.0 ||
	tmp_f[6] !=  5.5
       );
  gsl_test(s, "gsl_blas_raw_strsv I");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_strsv(CblasLower, CblasTrans, CblasNonUnit, 4, matrix_t_4_f, 4, tmp_f, 1);
  s = ( frac_diff(tmp_f[0], -53.0/6.0) > eps_f ||
        frac_diff(tmp_f[1],   5.0/3.0) > eps_f ||
        frac_diff(tmp_f[2],    1.5) > eps_f ||
	frac_diff(tmp_f[3],    3.0) > eps_f
       );
  gsl_test(s, "gsl_blas_raw_strsv J");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_strsv(CblasLower, CblasTrans, CblasUnit, 4, matrix_t_4_f, 4, tmp_f, 1);
  s = ( tmp_f[0] != -1.5 ||
        tmp_f[1] != -1.0 ||
	tmp_f[2] != -1.5 ||
	tmp_f[3] !=  3.0
       );
  gsl_test(s, "gsl_blas_raw_strsv K");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 2);
  gsl_blas_raw_strsv(CblasLower, CblasTrans, CblasUnit, 4, matrix_t_4_f, 4, tmp_f, 2);
  s = ( tmp_f[0] != -1.5 ||
        tmp_f[2] != -1.0 ||
	tmp_f[4] != -1.5 ||
	tmp_f[6] !=  3.0
       );
  gsl_test(s, "gsl_blas_raw_strsv L");
  status += s;


  gsl_blas_raw_dcopy(4, vector_4_d, 1, tmp_d, 1);
  gsl_blas_raw_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, 4, matrix_t_4_d, 4, tmp_d, 1);
  s = ( tmp_d[0] != -2.0 || frac_diff(tmp_d[1], -23.0/6.0) > eps ||
        tmp_d[2] !=  1.5 || tmp_d[3] != 3.0
       );
  gsl_test(s, "gsl_blas_raw_dtrsv A");
  status += s;

  gsl_blas_raw_dcopy(4, vector_4_d, 1, tmp_d, 1);
  gsl_blas_raw_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, 4, matrix_t_4_d, 4, tmp_d, 1);
  s = ( tmp_d[0] != -2.0 || tmp_d[1] != 0.0 ||
        tmp_d[2] != -4.0 || tmp_d[3] != 7.0
       );
  gsl_test(s, "gsl_blas_raw_dtrsv B");
  status += s;


  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctrsv(CblasUpper, CblasNoTrans, CblasNonUnit, 4, matrix_gen_4_c, 4, tmp_f, 1);
  s = ( tmp_f[0] !=  47.0/2.0 || tmp_f[1] !=  17.0/2.0 ||
        tmp_f[2] != -47.0/2.0 || tmp_f[3] != -33.0/2.0 ||
	tmp_f[4] !=  23.0/2.0 || tmp_f[5] !=  21.0/2.0 ||
	tmp_f[6] !=  5.0 || tmp_f[7] != -3.0
       );
  gsl_test(s, "gsl_blas_raw_ctrsv A");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 2);
  gsl_blas_raw_ctrsv(CblasUpper, CblasNoTrans, CblasNonUnit, 4, matrix_gen_4_c, 4, tmp_f, 2);
  s = ( tmp_f[0] !=  47.0/2.0 || tmp_f[1] !=  17.0/2.0 ||
        tmp_f[4] != -47.0/2.0 || tmp_f[5] != -33.0/2.0 ||
	tmp_f[8] !=  23.0/2.0 || tmp_f[9] !=  21.0/2.0 ||
	tmp_f[12] !=  5.0 || tmp_f[13] != -3.0
       );
  gsl_test(s, "gsl_blas_raw_ctrsv B");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctrsv(CblasUpper, CblasNoTrans, CblasUnit, 4, matrix_gen_4_c, 4, tmp_f, 1);
  s = ( tmp_f[0] !=  37.0/2.0 || tmp_f[1] != -75.0/2.0 ||
        tmp_f[2] != -67.0     || tmp_f[3] !=  24.0     ||
	tmp_f[4] !=  27.0/2.0 || tmp_f[5] != -17.0/2.0 ||
	tmp_f[6] !=  3.0 || tmp_f[7] != 5.0
       );
  gsl_test(s, "gsl_blas_raw_ctrsv C");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctrsv(CblasLower, CblasNoTrans, CblasNonUnit, 4, matrix_gen_4_c, 4, tmp_f, 1);
  s = ( tmp_f[0] !=  -2.0 || tmp_f[1] != 1.0 ||
        frac_diff(tmp_f[2], 0.05) > eps_f ||
	frac_diff(tmp_f[3], 0.15) > eps_f ||
	frac_diff(tmp_f[4],-51.0/10.0) > eps_f ||
	frac_diff(tmp_f[5],-33.0/10.0) > eps_f ||
	frac_diff(tmp_f[6], 64.0/5.0) > eps_f ||
	frac_diff(tmp_f[7],-28.0/5.0) > eps_f
       );
  gsl_test(s, "gsl_blas_raw_ctrsv D");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 2);
  gsl_blas_raw_ctrsv(CblasLower, CblasNoTrans, CblasNonUnit, 4, matrix_gen_4_c, 4, tmp_f, 2);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] != 1.0 ||
        frac_diff(tmp_f[4], 0.05) > eps_f ||
	frac_diff(tmp_f[5], 0.15) > eps_f ||
	frac_diff(tmp_f[8],-51.0/10.0) > eps_f ||
	frac_diff(tmp_f[9],-33.0/10.0) > eps_f ||
	frac_diff(tmp_f[12], 64.0/5.0) > eps_f ||
	frac_diff(tmp_f[13],-28.0/5.0) > eps_f
       );
  gsl_test(s, "gsl_blas_raw_ctrsv E");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctrsv(CblasLower, CblasNoTrans, CblasUnit, 4, matrix_gen_4_c, 4, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] !=  1.0 ||
        tmp_f[2] !=  0.0 || tmp_f[3] !=  0.5 ||
	tmp_f[4] !=  5.0 || tmp_f[5] !=  4.0 ||
	tmp_f[6] !=  8.5 || tmp_f[7] != -0.5
       );
  gsl_test(s, "gsl_blas_raw_ctrsv F");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctrsv(CblasUpper, CblasTrans, CblasNonUnit, 4, matrix_gen_4_c, 4, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] !=  1.0 ||
        frac_diff(tmp_f[2],  -1.0/5.0) > eps_f  ||
	frac_diff(tmp_f[3],   2.0/5.0) > eps_f  ||
	frac_diff(tmp_f[4],   8.0/5.0) > eps_f  ||
	frac_diff(tmp_f[5], -26.0/5.0) > eps_f  ||
	frac_diff(tmp_f[6],  21.0/5.0) > eps_f  ||
	frac_diff(tmp_f[7],  48.0/5.0) > eps_f
       );
  gsl_test(s, "gsl_blas_raw_ctrsv G");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctrsv(CblasUpper, CblasTrans, CblasUnit, 4, matrix_gen_4_c, 4, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] !=  1.0 ||
        tmp_f[2] != -1.0 || tmp_f[3] !=  1.0 ||
	tmp_f[4] !=  3.0 || tmp_f[5] !=  3.0 ||
	tmp_f[6] != 15.5 || tmp_f[7] != -2.5
       );
  gsl_test(s, "gsl_blas_raw_ctrsv H");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 2);
  gsl_blas_raw_ctrsv(CblasUpper, CblasTrans, CblasUnit, 4, matrix_gen_4_c, 4, tmp_f, 2);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] !=  1.0 ||
        tmp_f[4] != -1.0 || tmp_f[5] !=  1.0 ||
	tmp_f[8] !=  3.0 || tmp_f[9] !=  3.0 ||
	tmp_f[12] != 15.5 || tmp_f[13] != -2.5
       );
  gsl_test(s, "gsl_blas_raw_ctrsv I");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctrsv(CblasLower, CblasTrans, CblasNonUnit, 4, matrix_gen_4_c, 4, tmp_f, 1);
  s = ( frac_diff(tmp_f[0],-213.0/10.0) > eps_f ||
        frac_diff(tmp_f[1], -27.0/5.0) > eps_f  ||
	frac_diff(tmp_f[2],   8.0/5.0) > eps_f  ||
	frac_diff(tmp_f[3], -21.0/5.0) > eps_f  ||
	tmp_f[4] != 5.5 || tmp_f[5] !=  0.5 ||
	tmp_f[6] != 5.0 || tmp_f[7] != -3.0
       );
  gsl_test(s, "gsl_blas_raw_ctrsv J");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctrsv(CblasLower, CblasTrans, CblasUnit, 4, matrix_gen_4_c, 4, tmp_f, 1);
  s = ( tmp_f[0] != -19.0 || tmp_f[1] != -3.0 ||
        tmp_f[2] !=  19.0 || tmp_f[3] != -5.0 ||
	tmp_f[4] !=   3.5 || tmp_f[5] != -2.5 ||
	tmp_f[6] !=   3.0 || tmp_f[7] !=  5.0
       );
  gsl_test(s, "gsl_blas_raw_ctrsv K");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 2);
  gsl_blas_raw_ctrsv(CblasLower, CblasTrans, CblasUnit, 4, matrix_gen_4_c, 4, tmp_f, 2);
  s = ( tmp_f[0] != -19.0 || tmp_f[1] != -3.0 ||
        tmp_f[4] !=  19.0 || tmp_f[5] != -5.0 ||
	tmp_f[8] !=   3.5 || tmp_f[9] != -2.5 ||
	tmp_f[12] !=   3.0 || tmp_f[13] !=  5.0
       );
  gsl_test(s, "gsl_blas_raw_ctrsv L");
  status += s;


  gsl_blas_raw_zcopy(4, vector_4_z, 1, tmp_d, 2);
  gsl_blas_raw_ztrsv(CblasLower, CblasTrans, CblasUnit, 4, matrix_gen_4_z, 4, tmp_d, 2);
  s = ( tmp_d[0] != -19.0 || tmp_d[1] != -3.0 ||
        tmp_d[4] !=  19.0 || tmp_d[5] != -5.0 ||
	tmp_d[8] !=   3.5 || tmp_d[9] != -2.5 ||
	tmp_d[12] !=  3.0 || tmp_d[13] !=  5.0
       );
  gsl_test(s, "gsl_blas_raw_ztrsv A");
  status += s;


  return status;
}


int test_tbsv(void)
{
  int status = 0;
  int s;

  float  tmp_f[32];
  double tmp_d[32];


  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_stbsv(CblasUpper, CblasNoTrans, CblasNonUnit, 4, 3, matrix_t_4_f, 4, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || frac_diff(tmp_f[1], -23.0/6.0) > eps_f ||
        tmp_f[2] !=  1.5 || tmp_f[3] != 3.0
       );
  gsl_test(s, "gsl_blas_raw_stbsv A");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_stbsv(CblasUpper, CblasNoTrans, CblasNonUnit, 4, 2, matrix_t_4_f, 4, tmp_f, 1);
  s = ( tmp_f[0] != 1.0 || frac_diff(tmp_f[1], -23.0/6.0) > eps_f ||
        tmp_f[2] != 1.5 || tmp_f[3] != 3.0
       );
  gsl_test(s, "gsl_blas_raw_stbsv B");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 2);
  gsl_blas_raw_stbsv(CblasUpper, CblasNoTrans, CblasNonUnit, 4, 2, matrix_t_4_f, 4, tmp_f, 2);
  s = ( tmp_f[0] !=  1.0 || frac_diff(tmp_f[2], -23.0/6.0) > eps_f ||
        tmp_f[4] !=  1.5 || tmp_f[6] != 3.0
       );
  gsl_test(s, "gsl_blas_raw_stbsv C");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_stbsv(CblasUpper, CblasTrans, CblasNonUnit, 4, 2, matrix_t_4_f, 4, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 ||
        frac_diff(tmp_f[1], -1.0/3.0) > eps_f ||
	frac_diff(tmp_f[2],  7.0/3.0) > eps_f ||
	frac_diff(tmp_f[3], 13.0/6.0) > eps_f
       );
  gsl_test(s, "gsl_blas_raw_stbsv D");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_stbsv(CblasLower, CblasNoTrans, CblasNonUnit, 4, 2, matrix_t_4_f, 4, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] != 0.0 ||
        tmp_f[2] != -4.0 || tmp_f[3] != 5.0
       );
  gsl_test(s, "gsl_blas_raw_stbsv E");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_stbsv(CblasLower, CblasTrans, CblasNonUnit, 4, 2, matrix_t_4_f, 4, tmp_f, 1);
  s = ( frac_diff(tmp_f[0],-35.0/6.0) > eps_f ||
        frac_diff(tmp_f[1],  5.0/3.0) > eps_f ||
        tmp_f[2] !=  1.5 || tmp_f[3] != 3.0
       );
  gsl_test(s, "gsl_blas_raw_stbsv F");
  status += s;


  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctbsv(CblasUpper, CblasNoTrans, CblasNonUnit, 4, 2, matrix_gen_4_c, 4, tmp_f, 1);
  s = ( tmp_f[0] !=  63.0/2.0 || tmp_f[1] !=  21.0/2.0 ||
        tmp_f[2] != -47.0/2.0 || tmp_f[3] != -33.0/2.0 ||
	tmp_f[4] !=  23.0/2.0 || tmp_f[5] !=  21.0/2.0 ||
	tmp_f[6] !=  5.0 || tmp_f[7] != -3.0
       );
  gsl_test(s, "gsl_blas_raw_ctbsv A");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 2);
  gsl_blas_raw_ctbsv(CblasUpper, CblasNoTrans, CblasNonUnit, 4, 2, matrix_gen_4_c, 4, tmp_f, 2);
  s = ( tmp_f[0] !=  63.0/2.0 || tmp_f[1] !=  21.0/2.0 ||
        tmp_f[4] != -47.0/2.0 || tmp_f[5] != -33.0/2.0 ||
	tmp_f[8] !=  23.0/2.0 || tmp_f[9] !=  21.0/2.0 ||
	tmp_f[12] !=  5.0 || tmp_f[13] != -3.0
       );
  gsl_test(s, "gsl_blas_raw_ctbsv B");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctbsv(CblasUpper, CblasNoTrans, CblasUnit, 4, 2, matrix_gen_4_c, 4, tmp_f, 1);
  s = ( tmp_f[0] !=  33.0/2.0 || tmp_f[1] != -59.0/2.0 ||
        tmp_f[2] != -67.0     || tmp_f[3] !=  24.0     ||
	tmp_f[4] !=  27.0/2.0 || tmp_f[5] != -17.0/2.0 ||
	tmp_f[6] !=  3.0 || tmp_f[7] != 5.0
       );
  gsl_test(s, "gsl_blas_raw_ctbsv C");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctbsv(CblasLower, CblasNoTrans, CblasNonUnit, 4, 2, matrix_gen_4_c, 4, tmp_f, 1);
  s = ( tmp_f[0] !=  -2.0 || tmp_f[1] != 1.0 ||
        frac_diff(tmp_f[2], 0.05) > eps_f ||
	frac_diff(tmp_f[3], 0.15) > eps_f ||
	frac_diff(tmp_f[4],-51.0/10.0) > eps_f ||
	frac_diff(tmp_f[5],-33.0/10.0) > eps_f ||
	frac_diff(tmp_f[6], 59.0/5.0) > eps_f ||
	frac_diff(tmp_f[7],-13.0/5.0) > eps_f
       );
  gsl_test(s, "gsl_blas_raw_ctbsv D");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctbsv(CblasLower, CblasTrans, CblasNonUnit, 4, 2, matrix_gen_4_c, 4, tmp_f, 1);
  s = ( frac_diff(tmp_f[0], -133.0/10.0) > eps_f ||
	frac_diff(tmp_f[1], -17.0/5.0) > eps_f ||
	frac_diff(tmp_f[2],   8.0/5.0) > eps_f ||
	frac_diff(tmp_f[3], -21.0/5.0) > eps_f ||
	tmp_f[4] !=  5.5 || tmp_f[5] !=  0.5 ||
	tmp_f[6] !=  5.0 || tmp_f[7] != -3.0
       );
  gsl_test(s, "gsl_blas_raw_ctbsv E");
  status += s;


  return status;
}


int test_tpsv(void)
{
  int status = 0;
  int s;

  float  tmp_f[32];
  double tmp_d[32];


  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_stpsv(CblasUpper, CblasNoTrans, CblasNonUnit, 4, matrix_t_4_fpu, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || frac_diff(tmp_f[1], -23.0/6.0) > eps_f ||
        tmp_f[2] !=  1.5 || tmp_f[3] != 3.0
       );
  gsl_test(s, "gsl_blas_raw_stpsv A");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 2);
  gsl_blas_raw_stpsv(CblasUpper, CblasNoTrans, CblasNonUnit, 4, matrix_t_4_fpu, tmp_f, 2);
  s = ( tmp_f[0] != -2.0 || frac_diff(tmp_f[2], -23.0/6.0) > eps_f ||
        tmp_f[4] !=  1.5 || tmp_f[6] != 3.0
       );
  gsl_test(s, "gsl_blas_raw_stpsv B");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_stpsv(CblasUpper, CblasNoTrans, CblasUnit, 4, matrix_t_4_fpu, tmp_f, 1);
  s = ( tmp_f[0] != -8.0 || tmp_f[1] != 3.5 ||
        tmp_f[2] != -1.5 || tmp_f[3] != 3.0
       );
  gsl_test(s, "gsl_blas_raw_stpsv C");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_stpsv(CblasLower, CblasNoTrans, CblasNonUnit, 4, matrix_t_4_fpl, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] != 0.0 ||
        tmp_f[2] != -4.0 || tmp_f[3] != 7.0
       );
  gsl_test(s, "gsl_blas_raw_stpsv D");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 2);
  gsl_blas_raw_stpsv(CblasLower, CblasNoTrans, CblasNonUnit, 4, matrix_t_4_fpl, tmp_f, 2);
  s = ( tmp_f[0] != -2.0 || tmp_f[2] != 0.0 ||
        tmp_f[4] != -4.0 || tmp_f[6] != 7.0
       );
  gsl_test(s, "gsl_blas_raw_stpsv E");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_stpsv(CblasLower, CblasNoTrans, CblasUnit, 4, matrix_t_4_fpl, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] != 0.0 ||
        tmp_f[2] !=  4.0 || tmp_f[3] != 3.0
       );
  gsl_test(s, "gsl_blas_raw_stpsv F");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_stpsv(CblasUpper, CblasTrans, CblasNonUnit, 4, matrix_t_4_fpu, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 ||
        frac_diff(tmp_f[1], -1.0/3.0) > eps_f ||
        frac_diff(tmp_f[2],  7.0/3.0) > eps_f ||
	frac_diff(tmp_f[3], 25.0/6.0) > eps_f
       );
  gsl_test(s, "gsl_blas_raw_stpsv G");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_stpsv(CblasUpper, CblasTrans, CblasUnit, 4, matrix_t_4_fpu, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 ||
        tmp_f[1] != -1.0 ||
	tmp_f[2] !=  1.0 ||
	tmp_f[3] !=  5.5
       );
  gsl_test(s, "gsl_blas_raw_stpsv H");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 2);
  gsl_blas_raw_stpsv(CblasUpper, CblasTrans, CblasUnit, 4, matrix_t_4_fpu, tmp_f, 2);
  s = ( tmp_f[0] != -2.0 ||
        tmp_f[2] != -1.0 ||
	tmp_f[4] !=  1.0 ||
	tmp_f[6] !=  5.5
       );
  gsl_test(s, "gsl_blas_raw_stpsv I");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_stpsv(CblasLower, CblasTrans, CblasNonUnit, 4, matrix_t_4_fpl, tmp_f, 1);
  s = ( frac_diff(tmp_f[0], -53.0/6.0) > eps_f ||
        frac_diff(tmp_f[1],   5.0/3.0) > eps_f ||
        frac_diff(tmp_f[2],    1.5) > eps_f ||
	frac_diff(tmp_f[3],    3.0) > eps_f
       );
  gsl_test(s, "gsl_blas_raw_stpsv J");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_stpsv(CblasLower, CblasTrans, CblasUnit, 4, matrix_t_4_fpl, tmp_f, 1);
  s = ( tmp_f[0] != -1.5 ||
        tmp_f[1] != -1.0 ||
	tmp_f[2] != -1.5 ||
	tmp_f[3] !=  3.0
       );
  gsl_test(s, "gsl_blas_raw_stpsv K");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 2);
  gsl_blas_raw_stpsv(CblasLower, CblasTrans, CblasUnit, 4, matrix_t_4_fpl, tmp_f, 2);
  s = ( tmp_f[0] != -1.5 ||
        tmp_f[2] != -1.0 ||
	tmp_f[4] != -1.5 ||
	tmp_f[6] !=  3.0
       );
  gsl_test(s, "gsl_blas_raw_stpsv L");
  status += s;


  gsl_blas_raw_dcopy(4, vector_4_d, 1, tmp_d, 1);
  gsl_blas_raw_dtpsv(CblasUpper, CblasNoTrans, CblasNonUnit, 4, matrix_t_4_dpu, tmp_d, 1);
  s = ( tmp_d[0] != -2.0 || frac_diff(tmp_d[1], -23.0/6.0) > eps ||
        tmp_d[2] !=  1.5 || tmp_d[3] != 3.0
       );
  gsl_test(s, "gsl_blas_raw_dtpsv A");
  status += s;

  gsl_blas_raw_dcopy(4, vector_4_d, 1, tmp_d, 1);
  gsl_blas_raw_dtpsv(CblasLower, CblasNoTrans, CblasNonUnit, 4, matrix_t_4_dpl, tmp_d, 1);
  s = ( tmp_d[0] != -2.0 || tmp_d[1] != 0.0 ||
        tmp_d[2] != -4.0 || tmp_d[3] != 7.0
       );
  gsl_test(s, "gsl_blas_raw_dtpsv B");
  status += s;


  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctpsv(CblasUpper, CblasNoTrans, CblasNonUnit, 4, matrix_gen_4_cpu, tmp_f, 1);
  s = ( tmp_f[0] !=  47.0/2.0 || tmp_f[1] !=  17.0/2.0 ||
        tmp_f[2] != -47.0/2.0 || tmp_f[3] != -33.0/2.0 ||
	tmp_f[4] !=  23.0/2.0 || tmp_f[5] !=  21.0/2.0 ||
	tmp_f[6] !=  5.0 || tmp_f[7] != -3.0
       );
  gsl_test(s, "gsl_blas_raw_ctpsv A");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 2);
  gsl_blas_raw_ctpsv(CblasUpper, CblasNoTrans, CblasNonUnit, 4, matrix_gen_4_cpu, tmp_f, 2);
  s = ( tmp_f[0] !=  47.0/2.0 || tmp_f[1] !=  17.0/2.0 ||
        tmp_f[4] != -47.0/2.0 || tmp_f[5] != -33.0/2.0 ||
	tmp_f[8] !=  23.0/2.0 || tmp_f[9] !=  21.0/2.0 ||
	tmp_f[12] !=  5.0 || tmp_f[13] != -3.0
       );
  gsl_test(s, "gsl_blas_raw_ctpsv B");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctpsv(CblasUpper, CblasNoTrans, CblasUnit, 4, matrix_gen_4_cpu, tmp_f, 1);
  s = ( tmp_f[0] !=  37.0/2.0 || tmp_f[1] != -75.0/2.0 ||
        tmp_f[2] != -67.0     || tmp_f[3] !=  24.0     ||
	tmp_f[4] !=  27.0/2.0 || tmp_f[5] != -17.0/2.0 ||
	tmp_f[6] !=  3.0 || tmp_f[7] != 5.0
       );
  gsl_test(s, "gsl_blas_raw_ctpsv C");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctpsv(CblasLower, CblasNoTrans, CblasNonUnit, 4, matrix_gen_4_cpl, tmp_f, 1);
  s = ( tmp_f[0] !=  -2.0 || tmp_f[1] != 1.0 ||
        frac_diff(tmp_f[2], 0.05) > eps_f ||
	frac_diff(tmp_f[3], 0.15) > eps_f ||
	frac_diff(tmp_f[4],-51.0/10.0) > eps_f ||
	frac_diff(tmp_f[5],-33.0/10.0) > eps_f ||
	frac_diff(tmp_f[6], 64.0/5.0) > eps_f ||
	frac_diff(tmp_f[7],-28.0/5.0) > eps_f
       );
  gsl_test(s, "gsl_blas_raw_ctpsv D");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 2);
  gsl_blas_raw_ctpsv(CblasLower, CblasNoTrans, CblasNonUnit, 4, matrix_gen_4_cpl, tmp_f, 2);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] != 1.0 ||
        frac_diff(tmp_f[4], 0.05) > eps_f ||
	frac_diff(tmp_f[5], 0.15) > eps_f ||
	frac_diff(tmp_f[8],-51.0/10.0) > eps_f ||
	frac_diff(tmp_f[9],-33.0/10.0) > eps_f ||
	frac_diff(tmp_f[12], 64.0/5.0) > eps_f ||
	frac_diff(tmp_f[13],-28.0/5.0) > eps_f
       );
  gsl_test(s, "gsl_blas_raw_ctpsv E");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctpsv(CblasLower, CblasNoTrans, CblasUnit, 4, matrix_gen_4_cpl, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] !=  1.0 ||
        tmp_f[2] !=  0.0 || tmp_f[3] !=  0.5 ||
	tmp_f[4] !=  5.0 || tmp_f[5] !=  4.0 ||
	tmp_f[6] !=  8.5 || tmp_f[7] != -0.5
       );
  gsl_test(s, "gsl_blas_raw_ctpsv F");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctpsv(CblasUpper, CblasTrans, CblasNonUnit, 4, matrix_gen_4_cpu, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] !=  1.0 ||
        frac_diff(tmp_f[2],  -1.0/5.0) > eps_f  ||
	frac_diff(tmp_f[3],   2.0/5.0) > eps_f  ||
	frac_diff(tmp_f[4],   8.0/5.0) > eps_f  ||
	frac_diff(tmp_f[5], -26.0/5.0) > eps_f  ||
	frac_diff(tmp_f[6],  21.0/5.0) > eps_f  ||
	frac_diff(tmp_f[7],  48.0/5.0) > eps_f
       );
  gsl_test(s, "gsl_blas_raw_ctpsv G");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctpsv(CblasUpper, CblasTrans, CblasUnit, 4, matrix_gen_4_cpu, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] !=  1.0 ||
        tmp_f[2] != -1.0 || tmp_f[3] !=  1.0 ||
	tmp_f[4] !=  3.0 || tmp_f[5] !=  3.0 ||
	tmp_f[6] != 15.5 || tmp_f[7] != -2.5
       );
  gsl_test(s, "gsl_blas_raw_ctpsv H");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 2);
  gsl_blas_raw_ctpsv(CblasUpper, CblasTrans, CblasUnit, 4, matrix_gen_4_cpu, tmp_f, 2);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] !=  1.0 ||
        tmp_f[4] != -1.0 || tmp_f[5] !=  1.0 ||
	tmp_f[8] !=  3.0 || tmp_f[9] !=  3.0 ||
	tmp_f[12] != 15.5 || tmp_f[13] != -2.5
       );
  gsl_test(s, "gsl_blas_raw_ctpsv I");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctpsv(CblasLower, CblasTrans, CblasNonUnit, 4, matrix_gen_4_cpl, tmp_f, 1);
  s = ( frac_diff(tmp_f[0],-213.0/10.0) > eps_f ||
        frac_diff(tmp_f[1], -27.0/5.0) > eps_f  ||
	frac_diff(tmp_f[2],   8.0/5.0) > eps_f  ||
	frac_diff(tmp_f[3], -21.0/5.0) > eps_f  ||
	tmp_f[4] != 5.5 || tmp_f[5] !=  0.5 ||
	tmp_f[6] != 5.0 || tmp_f[7] != -3.0
       );
  gsl_test(s, "gsl_blas_raw_ctpsv J");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctpsv(CblasLower, CblasTrans, CblasUnit, 4, matrix_gen_4_cpl, tmp_f, 1);
  s = ( tmp_f[0] != -19.0 || tmp_f[1] != -3.0 ||
        tmp_f[2] !=  19.0 || tmp_f[3] != -5.0 ||
	tmp_f[4] !=   3.5 || tmp_f[5] != -2.5 ||
	tmp_f[6] !=   3.0 || tmp_f[7] !=  5.0
       );
  gsl_test(s, "gsl_blas_raw_ctpsv K");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 2);
  gsl_blas_raw_ctpsv(CblasLower, CblasTrans, CblasUnit, 4, matrix_gen_4_cpl, tmp_f, 2);
  s = ( tmp_f[0] != -19.0 || tmp_f[1] != -3.0 ||
        tmp_f[4] !=  19.0 || tmp_f[5] != -5.0 ||
	tmp_f[8] !=   3.5 || tmp_f[9] != -2.5 ||
	tmp_f[12] !=   3.0 || tmp_f[13] !=  5.0
       );
  gsl_test(s, "gsl_blas_raw_ctpsv L");
  status += s;


  gsl_blas_raw_zcopy(4, vector_4_z, 1, tmp_d, 2);
  gsl_blas_raw_ztpsv(CblasLower, CblasTrans, CblasUnit, 4, matrix_gen_4_zpl, tmp_d, 2);
  s = ( tmp_d[0] != -19.0 || tmp_d[1] != -3.0 ||
        tmp_d[4] !=  19.0 || tmp_d[5] != -5.0 ||
	tmp_d[8] !=   3.5 || tmp_d[9] != -2.5 ||
	tmp_d[12] !=  3.0 || tmp_d[13] !=  5.0
       );
  gsl_test(s, "gsl_blas_raw_ztpsv A");
  status += s;


  return status;
}


int test_symv(void)
{
  int status = 0;
  int s;

  float  tmp_f[32];
  double tmp_d[32];

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_ssymv(CblasUpper, 4, 2.0, matrix_gen_4_f, 4, vector_4_f, 1, 3.0, tmp_f, 1);
  s = ( tmp_f[0] != -4.0 || tmp_f[1] != -3.0 || tmp_f[2] != 1.0 || tmp_f[3] != 3.0 );
  gsl_test(s, "gsl_blas_raw_ssymv A");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_ssymv(CblasLower, 4, 2.0, matrix_gen_4_f, 4, vector_4_f, 1, 3.0, tmp_f, 1);
  s = ( tmp_f[0] != -5.0 || tmp_f[1] != -17.0 || tmp_f[2] != -1.0 || tmp_f[3] != 7.0 );
  gsl_test(s, "gsl_blas_raw_ssymv B");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 2);
  gsl_blas_raw_ssymv(CblasLower, 4, 2.0, matrix_gen_4_f, 4, vector_4_f, 1, 3.0, tmp_f, 2);
  s = ( tmp_f[0] != -5.0 || tmp_f[2] != -17.0 || tmp_f[4] != -1.0 || tmp_f[6] != 7.0 );
  gsl_test(s, "gsl_blas_raw_ssymv C");
  status += s;


  gsl_blas_raw_dcopy(4, vector_4_d, 1, tmp_d, 1);
  gsl_blas_raw_dsymv(CblasUpper, 4, 2.0, matrix_gen_4_d, 4, vector_4_d, 1, 3.0, tmp_d, 1);
  s = ( tmp_d[0] != -4.0 || tmp_d[1] != -3.0 || tmp_d[2] != 1.0 || tmp_d[3] != 3.0 );
  gsl_test(s, "gsl_blas_raw_dsymv A");
  status += s;


  return status;
}

int test_sbmv(void)
{
  int status = 0;
  int s;

  float  tmp_f[32];
  double tmp_d[32];

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_ssbmv(CblasUpper, 4, 2, 2.0, matrix_gen_4_f, 4, vector_4_f, 1, 3.0, tmp_f, 1);
  s = ( tmp_f[0] != -10.0 || tmp_f[1] != -3.0 || tmp_f[2] != 1.0 || tmp_f[3] != 7.0 );
  gsl_test(s, "gsl_blas_raw_ssbmv A");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_ssbmv(CblasUpper, 4, 1, 2.0, matrix_gen_4_f, 4, vector_4_f, 1, 3.0, tmp_f, 1);
  s = ( tmp_f[0] != -10.0 || tmp_f[1] != -9.0 || tmp_f[2] != -7.0 || tmp_f[3] != 9.0 );
  gsl_test(s, "gsl_blas_raw_ssbmv B");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_ssbmv(CblasLower, 4, 2, 2.0, matrix_gen_4_f, 4, vector_4_f, 1, 3.0, tmp_f, 1);
  s = ( tmp_f[0] != -11.0 || tmp_f[1] != -17.0 || tmp_f[2] != -1.0 || tmp_f[3] != 11.0 );
  gsl_test(s, "gsl_blas_raw_ssbmv C");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 2);
  gsl_blas_raw_ssbmv(CblasLower, 4, 2, 2.0, matrix_gen_4_f, 4, vector_4_c, 2, 3.0, tmp_f, 2);
  s = ( tmp_f[0] != -11.0 || tmp_f[2] != -17.0 || tmp_f[4] != -1.0 || tmp_f[6] != 11.0 );
  gsl_test(s, "gsl_blas_raw_ssbmv D");
  status += s;


  gsl_blas_raw_dcopy(4, vector_4_d, 1, tmp_d, 2);
  gsl_blas_raw_dsbmv(CblasUpper, 4, 2, 2.0, matrix_gen_4_d, 4, vector_4_d, 1, 3.0, tmp_d, 2);
  s = ( tmp_d[0] != -10.0 || tmp_d[2] != -3.0 || tmp_d[4] != 1.0 || tmp_d[6] != 7.0 );
  gsl_test(s, "gsl_blas_raw_dsbmv A");
  status += s;


  return status;
}

int test_spmv(void)
{
  int status = 0;
  int s;

  float  tmp_f[32];
  double tmp_d[32];


  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_sspmv(CblasUpper, 4, 2.0, matrix_gen_4_f, vector_4_f, 1, 3.0, tmp_f, 1);
  s = ( tmp_f[0] != -4.0 || tmp_f[1] != 26.0 || tmp_f[2] != 14.0 || tmp_f[3] != -17.0 );
  gsl_test(s, "gsl_blas_raw_sspmv A");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_sspmv(CblasLower, 4, 2.0, matrix_gen_4_f, vector_4_f, 1, 3.0, tmp_f, 1);
  s = ( tmp_f[0] != 20.0 || tmp_f[1] != 7.0 || tmp_f[2] != 7.0 || tmp_f[3] != -25.0 );
  gsl_test(s, "gsl_blas_raw_sspmv B");

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 2);
  gsl_blas_raw_sspmv(CblasLower, 4, 2.0, matrix_gen_4_f, vector_4_f, 1, 3.0, tmp_f, 2);
  s = ( tmp_f[0] != 20.0 || tmp_f[2] != 7.0 || tmp_f[4] != 7.0 || tmp_f[6] != -25.0 );
  gsl_test(s, "gsl_blas_raw_sspmv C");


  gsl_blas_raw_dcopy(4, vector_4_d, 1, tmp_d, 1);
  gsl_blas_raw_dspmv(CblasUpper, 4, 2.0, matrix_gen_4_d, vector_4_d, 1, 3.0, tmp_d, 1);
  s = ( tmp_d[0] != -4.0 || tmp_d[1] != 26.0 || tmp_d[2] != 14.0 || tmp_d[3] != -17.0 );
  gsl_test(s, "gsl_blas_raw_dspmv A");
  status += s;

  return status;
}

int test_ger(void)
{
  int status = 0;
  int s;

  float  tmp_f[32];
  double tmp_d[32];

  gsl_blas_raw_scopy(4*4, matrix_gen_4_f, 1, tmp_f, 1);
  gsl_blas_raw_sger (4, 4, 2.0, vector_4_f, 1, vector_4_c, 1, tmp_f, 4);
  s = ( tmp_f[0] != 9.0 || tmp_f[1] != 4.0 || tmp_f[2] != -2.0 || tmp_f[3] != -11.0 ||
        tmp_f[4] != -7.0/2.0 || tmp_f[5] != 1.0 || tmp_f[6] != 5.0 || tmp_f[7] != 7.0 ||
	tmp_f[8] != 6.0 || tmp_f[9] != 0.0 || tmp_f[10] != -1.0 || tmp_f[11] != -11.0/2.0 ||
        tmp_f[12] != -3.0 || tmp_f[13] != -3.0 || tmp_f[14] != 0.5 || tmp_f[15] != 6.0
	);
  gsl_test(s, "gsl_blas_raw_sger A");
  status += s;

  gsl_blas_raw_scopy(4*4, matrix_gen_4_f, 1, tmp_f, 1);
  gsl_blas_raw_sger (4, 4, 2.0, vector_4_c, 2, vector_4_c, 1, tmp_f, 4);
  s = ( tmp_f[0] != 9.0 || tmp_f[1] != 4.0 || tmp_f[2] != -2.0 || tmp_f[3] != -11.0 ||
        tmp_f[4] != -7.0/2.0 || tmp_f[5] != 1.0 || tmp_f[6] != 5.0 || tmp_f[7] != 7.0 ||
	tmp_f[8] != 6.0 || tmp_f[9] != 0.0 || tmp_f[10] != -1.0 || tmp_f[11] != -11.0/2.0 ||
        tmp_f[12] != -3.0 || tmp_f[13] != -3.0 || tmp_f[14] != 0.5 || tmp_f[15] != 6.0
	);
  gsl_test(s, "gsl_blas_raw_sger B");
  status += s;


  gsl_blas_raw_dcopy(4*4, matrix_gen_4_d, 1, tmp_d, 1);
  gsl_blas_raw_dger(4, 4, 2.0, vector_4_d, 1, vector_4_z, 1, tmp_d, 4);
  s = ( tmp_d[0] != 9.0 || tmp_d[1] != 4.0 || tmp_d[2] != -2.0 || tmp_d[3] != -11.0 ||
        tmp_d[4] != -7.0/2.0 || tmp_d[5] != 1.0 || tmp_d[6] != 5.0 || tmp_d[7] != 7.0 ||
	tmp_d[8] != 6.0 || tmp_d[9] != 0.0 || tmp_d[10] != -1.0 || tmp_d[11] != -11.0/2.0 ||
        tmp_d[12] != -3.0 || tmp_d[13] != -3.0 || tmp_d[14] != 0.5 || tmp_d[15] != 6.0
	);
  gsl_test(s, "gsl_blas_raw_dger A");
  status += s;


  return status;
}

int test_syr(void)
{
  int status = 0;
  int s;

  float  tmp_f[32];
  double tmp_d[32];

  gsl_blas_raw_scopy(4*4, matrix_gen_4_f, 1, tmp_f, 1);
  gsl_blas_raw_ssyr(CblasUpper, 4, 2.0, vector_4_f, 1, tmp_f, 4);
  s = ( tmp_f[0] != 9.0 || tmp_f[1] != 4.0 || tmp_f[2] != -2.0 || tmp_f[3] != -11.0 ||
        tmp_f[5] != 5.0 || tmp_f[6] != 5.0 || tmp_f[7] != -5.0 ||
	tmp_f[10] != -1.0 || tmp_f[11] != 0.5 ||
	tmp_f[15] != 18.0
       );
  gsl_test(s, "gsl_blas_raw_ssyr A");
  status += s;

  gsl_blas_raw_scopy(4*4, matrix_gen_4_f, 1, tmp_f, 1);
  gsl_blas_raw_ssyr(CblasLower, 4, 2.0, vector_4_f, 1, tmp_f, 4);
  s = ( tmp_f[0] != 9.0 ||
        tmp_f[4] != 4.5 || tmp_f[5] !=  5.0 ||
	tmp_f[8] != 2.0 || tmp_f[9] != -2.0 || tmp_f[10] != -1.0 ||
	tmp_f[12] != -11.0 || tmp_f[13] != -7.0 || tmp_f[14] != 0.5 || tmp_f[15] != 18.0
       );
  gsl_test(s, "gsl_blas_raw_ssyr B");
  status += s;

  gsl_blas_raw_scopy(4*4, matrix_gen_4_f, 1, tmp_f, 1);
  gsl_blas_raw_ssyr(CblasLower, 4, 2.0, vector_4_c, 2, tmp_f, 4);
  s = ( tmp_f[0] != 9.0 ||
        tmp_f[4] != 4.5 || tmp_f[5] !=  5.0 ||
	tmp_f[8] != 2.0 || tmp_f[9] != -2.0 || tmp_f[10] != -1.0 ||
	tmp_f[12] != -11.0 || tmp_f[13] != -7.0 || tmp_f[14] != 0.5 || tmp_f[15] != 18.0
       );
  gsl_test(s, "gsl_blas_raw_ssyr B");
  status += s;


  gsl_blas_raw_dcopy(4*4, matrix_gen_4_d, 1, tmp_d, 1);
  gsl_blas_raw_dsyr(CblasUpper, 4, 2.0, vector_4_d, 1, tmp_d, 4);
  s = ( tmp_d[0] != 9.0 || tmp_d[1] != 4.0 || tmp_d[2] != -2.0 || tmp_d[3] != -11.0 ||
        tmp_d[5] != 5.0 || tmp_d[6] != 5.0 || tmp_d[7] != -5.0 ||
	tmp_d[10] != -1.0 || tmp_d[11] != 0.5 ||
	tmp_d[15] != 18.0
       );
  gsl_test(s, "gsl_blas_raw_dsyr A");
  status += s;


  return status;
}

int test_spr(void)
{
  int status = 0;
  int s;

  float  tmp_f[32];
  double tmp_d[32];

  gsl_blas_raw_scopy(4*4, matrix_gen_4_f, 1, tmp_f, 1);
  gsl_blas_raw_sspr(CblasUpper, 4, 2.0, vector_4_f, 1, tmp_f);
  s = ( tmp_f[0] != 9.0 || tmp_f[1] != 4.0 || tmp_f[2] != -2.0 || tmp_f[3] != -11.0 ||
        tmp_f[4] != 5.0/2.0 || tmp_f[5] != 3.0 || tmp_f[6] != -1.0 ||
	tmp_f[7] != 1.0 || tmp_f[8] != 2.0 ||
	tmp_f[9] != 16.0
       );
  gsl_test(s, "gsl_blas_raw_sspr A");
  status += s;

  gsl_blas_raw_scopy(4*4, matrix_gen_4_f, 1, tmp_f, 1);
  gsl_blas_raw_sspr(CblasLower, 4, 2.0, vector_4_f, 1, tmp_f);
  s = ( tmp_f[0] != 9.0 ||
        tmp_f[1] != 4.0 || tmp_f[2] != 0.0 ||
	tmp_f[3] != 1.0 || tmp_f[4] != 0.5 || tmp_f[5] != 3.0 ||
	tmp_f[6] != -7.0 || tmp_f[7] != -5.0 || tmp_f[8] != 2.0 || tmp_f[9] != 16.0
       );
  gsl_test(s, "gsl_blas_raw_sspr B");
  status += s;

  gsl_blas_raw_scopy(4*4, matrix_gen_4_f, 1, tmp_f, 1);
  gsl_blas_raw_sspr(CblasLower, 4, 2.0, vector_4_c, 2, tmp_f);
  s = ( tmp_f[0] != 9.0 ||
        tmp_f[1] != 4.0 || tmp_f[2] != 0.0 ||
	tmp_f[3] != 1.0 || tmp_f[4] != 0.5 || tmp_f[5] != 3.0 ||
	tmp_f[6] != -7.0 || tmp_f[7] != -5.0 || tmp_f[8] != 2.0 || tmp_f[9] != 16.0
       );
  gsl_test(s, "gsl_blas_raw_sspr C");
  status += s;

  gsl_blas_raw_dcopy(4*4, matrix_gen_4_d, 1, tmp_d, 1);
  gsl_blas_raw_dspr(CblasUpper, 4, 2.0, vector_4_d, 1, tmp_d);
  s = ( tmp_d[0] != 9.0 || tmp_d[1] != 4.0 || tmp_d[2] != -2.0 || tmp_d[3] != -11.0 ||
        tmp_d[4] != 5.0/2.0 || tmp_d[5] != 3.0 || tmp_d[6] != -1.0 ||
	tmp_d[7] != 1.0 || tmp_d[8] != 2.0 ||
	tmp_d[9] != 16.0
       );
  gsl_test(s, "gsl_blas_raw_dspr A");
  status += s;


  return status;
}

int test_syr2(void)
{
  int status = 0;
  int s;

  float  tmp_f[32];
  double tmp_d[32];

  gsl_blas_raw_scopy(4*4, matrix_gen_4_f, 1, tmp_f, 1);
  gsl_blas_raw_ssyr2(CblasUpper, 4, 2.0, vector_4_f, 1, vector_4_c, 1, tmp_f, 4);
  s = ( tmp_f[0] != 17.0 || tmp_f[1] != 0.0 || tmp_f[2] != 2.0 || tmp_f[3] != -15.0 ||
        tmp_f[5] != -1.0 || tmp_f[6] != 7.0 || tmp_f[7] != 5.0 ||
	tmp_f[10] != -1.0 || tmp_f[11] != -5.5 ||
	tmp_f[15] != 12.0
       );
  gsl_test(s, "gsl_blas_raw_ssyr2 A");
  status += s;

  gsl_blas_raw_scopy(4*4, matrix_gen_4_f, 1, tmp_f, 1);
  gsl_blas_raw_ssyr2(CblasLower, 4, 2.0, vector_4_f, 1, vector_4_c, 1, tmp_f, 4);
  s = ( tmp_f[0] != 17.0 ||
        tmp_f[4] !=  0.5 || tmp_f[5] != -1.0 ||
	tmp_f[8] !=  6.0 || tmp_f[9] !=  0.0 || tmp_f[10] != -1.0 ||
	tmp_f[12] != -15.0 || tmp_f[13] != 3.0 || tmp_f[14] != -5.5 || tmp_f[15] != 12.0
       );
  gsl_test(s, "gsl_blas_raw_ssyr2 B");
  status += s;

  gsl_blas_raw_scopy(4*4, matrix_gen_4_f, 1, tmp_f, 1);
  gsl_blas_raw_ssyr2(CblasLower, 4, 2.0, vector_4_c, 2, vector_4_c, 1, tmp_f, 4);
  s = ( tmp_f[0] != 17.0 ||
        tmp_f[4] !=  0.5 || tmp_f[5] != -1.0 ||
	tmp_f[8] !=  6.0 || tmp_f[9] !=  0.0 || tmp_f[10] != -1.0 ||
	tmp_f[12] != -15.0 || tmp_f[13] != 3.0 || tmp_f[14] != -5.5 || tmp_f[15] != 12.0
       );
  gsl_test(s, "gsl_blas_raw_ssyr2 C");
  status += s;


  gsl_blas_raw_dcopy(4*4, matrix_gen_4_d, 1, tmp_d, 1);
  gsl_blas_raw_dsyr2(CblasUpper, 4, 2.0, vector_4_d, 1, vector_4_z, 1, tmp_d, 4);
  s = ( tmp_d[0] != 17.0 || tmp_d[1] != 0.0 || tmp_d[2] != 2.0 || tmp_d[3] != -15.0 ||
        tmp_d[5] != -1.0 || tmp_d[6] != 7.0 || tmp_d[7] != 5.0 ||
	tmp_d[10] != -1.0 || tmp_d[11] != -5.5 ||
	tmp_d[15] != 12.0
       );
  gsl_test(s, "gsl_blas_raw_dsyr2 A");
  status += s;


  return status;
}

int test_spr2(void)
{
  int status = 0;
  int s;

  float  tmp_f[32];
  double tmp_d[32];


  gsl_blas_raw_scopy(10, matrix_gen_4_f, 1, tmp_f, 1);
  gsl_blas_raw_sspr2(CblasUpper, 4, 2.0, vector_4_f, 1, vector_4_c, 1, tmp_f);
  s = ( tmp_f[0] != 17.0 || tmp_f[1] !=  0.0 || tmp_f[2] != 2.0 || tmp_f[3] != -15.0 ||
        tmp_f[4] != -3.5 || tmp_f[5] !=  5.0 || tmp_f[6] != 9.0 ||
	tmp_f[7] !=  1.0 || tmp_f[8] != -4.0 ||
	tmp_f[9] != 10.0
       );
  gsl_test(s, "gsl_blas_raw_sspr2 A");
  status += s;

  gsl_blas_raw_scopy(10, matrix_gen_4_f, 1, tmp_f, 1);
  gsl_blas_raw_sspr2(CblasLower, 4, 2.0, vector_4_f, 1, vector_4_c, 1, tmp_f);
  s = ( tmp_f[0] !=  17.0 ||
        tmp_f[1] !=   0.0 || tmp_f[2] != -6.0 ||
	tmp_f[3] !=   5.0 || tmp_f[4] !=  2.5 || tmp_f[5] !=  3.0 ||
	tmp_f[6] != -11.0 || tmp_f[7] !=  5.0 || tmp_f[8] != -4.0 || tmp_f[9] != 10.0
       );
  gsl_test(s, "gsl_blas_raw_sspr2 B");
  status += s;

  gsl_blas_raw_scopy(10, matrix_gen_4_f, 1, tmp_f, 1);
  gsl_blas_raw_sspr2(CblasLower, 4, 2.0, vector_4_f, 1, vector_4_c, 2, tmp_f);
  s = ( tmp_f[0] !=  17.0 ||
        tmp_f[1] !=   8.0 || tmp_f[2] !=   2.0 ||
	tmp_f[3] !=   1.0 || tmp_f[4] !=   0.5 || tmp_f[5] !=  3.0 ||
	tmp_f[6] != -19.0 || tmp_f[7] != -11.0 || tmp_f[8] !=  2.0 || tmp_f[9] != 34.0
       );
  gsl_test(s, "gsl_blas_raw_sspr2 C");
  status += s;


  gsl_blas_raw_dcopy(4*4, matrix_gen_4_d, 1, tmp_d, 1);
  gsl_blas_raw_dspr2(CblasUpper, 4, 2.0, vector_4_d, 1, vector_4_z, 1, tmp_d);
  s = ( tmp_d[0] != 17.0 || tmp_d[1] !=  0.0 || tmp_d[2] != 2.0 || tmp_d[3] != -15.0 ||
        tmp_d[4] != -3.5 || tmp_d[5] !=  5.0 || tmp_d[6] != 9.0 ||
	tmp_d[7] !=  1.0 || tmp_d[8] != -4.0 ||
	tmp_d[9] != 10.0
       );
  gsl_test(s, "gsl_blas_raw_dspr2 A");
  status += s;

  
  return status;
}

int test_hemv(void)
{
  int status = 0;
  int s;

  float  tmp_f[32];
  double tmp_d[32];

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_chemv(CblasUpper, 4, c_2, matrix_her_4_c, 4, vector_4_c, 1, c_3, tmp_f, 1);
  s = ( tmp_f[0] != -20.0 || tmp_f[1] != 9.0 || tmp_f[2] != -35.0 || tmp_f[3] != 51.0 ||
        tmp_f[4] != -25.0 || tmp_f[5] != 38.0 || tmp_f[6] != 35.0 || tmp_f[7] != 38.0
       );
  gsl_test(s, "gsl_blas_raw_chemv A");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_chemv(CblasLower, 4, c_2, matrix_her_4_c, 4, vector_4_c, 1, c_3, tmp_f, 1);
  s = ( tmp_f[0] != -20.0 || tmp_f[1] != 9.0 || tmp_f[2] != -35.0 || tmp_f[3] != 51.0 ||
        tmp_f[4] != -25.0 || tmp_f[5] != 38.0 || tmp_f[6] != 35.0 || tmp_f[7] != 38.0
       );
  gsl_test(s, "gsl_blas_raw_chemv B");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 2);
  gsl_blas_raw_chemv(CblasUpper, 4, c_2, matrix_her_4_c, 4, vector_4_c, 1, c_3, tmp_f, 2);
  s = ( tmp_f[0] != -20.0 || tmp_f[1] != 9.0 || tmp_f[4] != -35.0 || tmp_f[5] != 51.0 ||
        tmp_f[8] != -25.0 || tmp_f[9] != 38.0 || tmp_f[12] != 35.0 || tmp_f[13] != 38.0
       );
  gsl_test(s, "gsl_blas_raw_chemv C");
  status += s;


  gsl_blas_raw_zcopy(4, vector_4_z, 1, tmp_d, 1);
  gsl_blas_raw_zhemv(CblasUpper, 4, z_2, matrix_her_4_z, 4, vector_4_z, 1, z_3, tmp_d, 1);
  s = ( tmp_d[0] != -20.0 || tmp_d[1] != 9.0 || tmp_d[2] != -35.0 || tmp_d[3] != 51.0 ||
        tmp_d[4] != -25.0 || tmp_d[5] != 38.0 || tmp_d[6] != 35.0 || tmp_d[7] != 38.0
       );
  gsl_test(s, "gsl_blas_raw_zhemv A");
  status += s;


  return status;
}

int test_hbmv(void)
{
  int status = 0;
  int s;

  float  tmp_f[32];
  double tmp_d[32];

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_chbmv(CblasUpper, 4, 2, c_2, matrix_her_4_c, 4, vector_4_c, 1, c_3, tmp_f, 1);
  s = ( tmp_f[0] != -16.0 || tmp_f[1] != -7.0 || tmp_f[2] != -35.0 || tmp_f[3] != 51.0 ||
        tmp_f[4] != -25.0 || tmp_f[5] != 38.0 || tmp_f[6] != 37.0 || tmp_f[7] != 32.0
       );
  gsl_test(s, "gsl_blas_raw_chbmv A");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_chbmv(CblasLower, 4, 2, c_2, matrix_her_4_c, 4, vector_4_c, 1, c_3, tmp_f, 1);
  s = ( tmp_f[0] != -16.0 || tmp_f[1] != -7.0 || tmp_f[2] != -35.0 || tmp_f[3] != 51.0 ||
        tmp_f[4] != -25.0 || tmp_f[5] != 38.0 || tmp_f[6] != 37.0 || tmp_f[7] != 32.0
       );
  gsl_test(s, "gsl_blas_raw_chbmv B");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 2);
  gsl_blas_raw_chbmv(CblasUpper, 4, 2, c_2, matrix_her_4_c, 4, vector_4_c, 1, c_3, tmp_f, 2);
  s = ( tmp_f[0] != -16.0 || tmp_f[1] != -7.0 || tmp_f[4] != -35.0 || tmp_f[5] != 51.0 ||
        tmp_f[8] != -25.0 || tmp_f[9] != 38.0 || tmp_f[12] != 37.0 || tmp_f[13] != 32.0
       );
  gsl_test(s, "gsl_blas_raw_chbmv C");
  status += s;


  gsl_blas_raw_zcopy(4, vector_4_z, 1, tmp_d, 1);
  gsl_blas_raw_zhbmv(CblasUpper, 4, 2, z_2, matrix_her_4_z, 4, vector_4_z, 1, z_3, tmp_d, 1);
  s = ( tmp_d[0] != -16.0 || tmp_d[1] != -7.0 || tmp_d[2] != -35.0 || tmp_d[3] != 51.0 ||
        tmp_d[4] != -25.0 || tmp_d[5] != 38.0 || tmp_d[6] != 37.0 || tmp_d[7] != 32.0
       );
  gsl_test(s, "gsl_blas_raw_zhbmv A");
  status += s;


  return status;
}

int test_hpmv(void)
{
  int status = 0;
  int s;

  float  tmp_f[32];
  double tmp_d[32];


  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_chpmv(CblasUpper, 4, c_2, matrix_her_4_cpu, vector_4_c, 1, c_3, tmp_f, 1);
  s = ( tmp_f[0] != -20.0 || tmp_f[1] !=  9.0 || tmp_f[2] != -35.0 || tmp_f[3] != 51.0 ||
        tmp_f[4] != -25.0 || tmp_f[5] != 38.0 || tmp_f[6] !=  35.0 || tmp_f[7] != 38.0
       );
  gsl_test(s, "gsl_blas_raw_chpmv A");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_chpmv(CblasLower, 4, c_2, matrix_her_4_cpl, vector_4_c, 1, c_3, tmp_f, 1);
  s = ( tmp_f[0] != -20.0 || tmp_f[1] !=  9.0 || tmp_f[2] != -35.0 || tmp_f[3] != 51.0 ||
        tmp_f[4] != -25.0 || tmp_f[5] != 38.0 || tmp_f[6] !=  35.0 || tmp_f[7] != 38.0
       );
  gsl_test(s, "gsl_blas_raw_chpmv B");
  status += s;


  gsl_blas_raw_zcopy(4, vector_4_z, 1, tmp_d, 1);
  gsl_blas_raw_zhpmv(CblasUpper, 4, z_2, matrix_her_4_zpu, vector_4_z, 1, z_3, tmp_d, 1);
  s = ( tmp_d[0] != -20.0 || tmp_d[1] !=  9.0 || tmp_d[2] != -35.0 || tmp_d[3] != 51.0 ||
        tmp_d[4] != -25.0 || tmp_d[5] != 38.0 || tmp_d[6] !=  35.0 || tmp_d[7] != 38.0
       );
  gsl_test(s, "gsl_blas_raw_zhpmv A");
  status += s;


  return status;
}

int test_geru(void)
{
  int status = 0;
  int s;

  float  tmp_f[32];
  double tmp_d[32];


  gsl_blas_raw_ccopy(4*4, matrix_gen_4_c, 1, tmp_f, 1);
  gsl_blas_raw_cgeru(4, 4, c_2, vector_4_c, 1, vector_4_c, 1, tmp_f, 4);
  s = ( tmp_f[0] !=  7.0 || tmp_f[1] !=  -8.0 || tmp_f[2] !=   2.0 || tmp_f[3]  !=  -6.0 ||
        tmp_f[4] != -8.0 || tmp_f[5] != -11.0 || tmp_f[6] != -21.0 || tmp_f[7]  != -13.0 ||
	tmp_f[8] !=  2.5 || tmp_f[9] !=  -6.0 || tmp_f[10] !=  3.0 || tmp_f[11] !=  -3.0 ||
        tmp_f[12] != -1.0 || tmp_f[13] != -5.0|| tmp_f[14] != -16.0 || tmp_f[15] != -2.0
	);
  gsl_test(s, "gsl_blas_raw_cgeru A");
  status += s;


  gsl_blas_raw_zcopy(4*4, matrix_gen_4_z, 1, tmp_d, 1);
  gsl_blas_raw_zgeru(4, 4, z_2, vector_4_z, 1, vector_4_z, 1, tmp_d, 4);
  s = ( tmp_d[0] !=  7.0 || tmp_d[1] !=  -8.0 || tmp_d[2] !=   2.0 || tmp_d[3]  !=  -6.0 ||
        tmp_d[4] != -8.0 || tmp_d[5] != -11.0 || tmp_d[6] != -21.0 || tmp_d[7]  != -13.0 ||
	tmp_d[8] !=  2.5 || tmp_d[9] !=  -6.0 || tmp_d[10] !=  3.0 || tmp_d[11] !=  -3.0 ||
        tmp_d[12] != -1.0 || tmp_d[13] != -5.0|| tmp_d[14] != -16.0 || tmp_d[15] != -2.0
	);
  gsl_test(s, "gsl_blas_raw_zgeru A");
  status += s;


  return status;
}

int test_gerc(void)
{
  int status = 0;
  int s;

  float  tmp_f[32];
  double tmp_d[32];

  gsl_blas_raw_ccopy(4*4, matrix_gen_4_c, 1, tmp_f, 1);
  gsl_blas_raw_cgerc(4, 4, c_2, vector_4_c, 1, vector_4_c, 1, tmp_f, 4);
  s = ( tmp_f[0] != 11.0 || tmp_f[1] !=  0.0 || tmp_f[2] !=  6.0 || tmp_f[3] != 2.0 ||
        tmp_f[4] !=  4.0 || tmp_f[5] != 13.0 || tmp_f[6] != -1.0 || tmp_f[7] != 27.0 ||
	tmp_f[8] !=  6.5 || tmp_f[9] != -2.0 || tmp_f[10] != 7.0 || tmp_f[11] !=  1.0 ||
        tmp_f[12] != 11.0 || tmp_f[13] != 7.0|| tmp_f[14] != 4.0 || tmp_f[15] != 18.0
	);
  gsl_test(s, "gsl_blas_raw_cgerc A");
  status += s;


  gsl_blas_raw_zcopy(4*4, matrix_gen_4_z, 1, tmp_d, 1);
  gsl_blas_raw_zgerc(4, 4, z_2, vector_4_z, 1, vector_4_z, 1, tmp_d, 4);
  s = ( tmp_d[0] != 11.0 || tmp_d[1] !=  0.0 || tmp_d[2] !=  6.0 || tmp_d[3] != 2.0 ||
        tmp_d[4] !=  4.0 || tmp_d[5] != 13.0 || tmp_d[6] != -1.0 || tmp_d[7] != 27.0 ||
	tmp_d[8] !=  6.5 || tmp_d[9] != -2.0 || tmp_d[10] != 7.0 || tmp_d[11] !=  1.0 ||
        tmp_d[12] != 11.0 || tmp_d[13] != 7.0|| tmp_d[14] != 4.0 || tmp_d[15] != 18.0
	);
  gsl_test(s, "gsl_blas_raw_zgerc A");
  status += s;


  return status;
}

int test_her(void)
{
  int status = 0;
  int s;

  float  tmp_f[32];
  double tmp_d[32];

  gsl_blas_raw_ccopy(4*4, matrix_her_4_c, 1, tmp_f, 1);
  gsl_blas_raw_cher(CblasUpper, 4, 2.0, vector_4_c, 1, tmp_f, 4);
  s = ( tmp_f[0] != 11.0 || tmp_f[1] !=  0.0 ||
        tmp_f[2] !=  6.0 || tmp_f[3] !=  2.0 ||
        tmp_f[4] !=  4.0 || tmp_f[5] != 13.0 ||
	tmp_f[6] != -1.0 || tmp_f[7] != 27.0 ||
	tmp_f[10] != 7.0 || tmp_f[11] != 0.0
       );
  gsl_test(s, "gsl_blas_raw_cher A");
  status += s;

  gsl_blas_raw_ccopy(4*4, matrix_her_4_c, 1, tmp_f, 1);
  gsl_blas_raw_cher(CblasLower, 4, 2.0, vector_4_c, 1, tmp_f, 4);
  s = ( tmp_f[0] != 11.0 || tmp_f[1] !=  0.0  ||
        tmp_f[8] !=  6.0 || tmp_f[9] != -2.0  ||
        tmp_f[10] !=  7.0 || tmp_f[11] !=  0.0  ||
	tmp_f[16] !=  4.0 || tmp_f[17] != -13.0 ||
	tmp_f[20] != 17.0 || tmp_f[21] !=  0.0
       );
  gsl_test(s, "gsl_blas_raw_cher B");
  status += s;


  gsl_blas_raw_zcopy(4*4, matrix_her_4_z, 1, tmp_d, 1);
  gsl_blas_raw_zher(CblasUpper, 4, 2.0, vector_4_z, 1, tmp_d, 4);
  s = ( tmp_d[0] != 11.0 || tmp_d[1] !=  0.0 ||
        tmp_d[2] !=  6.0 || tmp_d[3] !=  2.0 ||
        tmp_d[4] !=  4.0 || tmp_d[5] != 13.0 ||
	tmp_d[6] != -1.0 || tmp_d[7] != 27.0 ||
	tmp_d[10] != 7.0 || tmp_d[11] != 0.0
       );
  gsl_test(s, "gsl_blas_raw_zher A");
  status += s;


  return status;
}

int test_hpr(void)
{
  int status = 0;
  int s;

  float  tmp_f[32];
  double tmp_d[32];

  gsl_blas_raw_ccopy(10, matrix_her_4_cpu, 1, tmp_f, 1);
  gsl_blas_raw_chpr(CblasUpper, 4, 2.0, vector_4_c, 1, tmp_f);
  s = ( tmp_f[0] != 11.0 || tmp_f[1] !=  0.0 ||
        tmp_f[2] !=  6.0 || tmp_f[3] !=  2.0 ||
        tmp_f[4] !=  4.0 || tmp_f[5] != 13.0 ||
	tmp_f[6] != -1.0 || tmp_f[7] != 27.0 ||
	tmp_f[8] !=  7.0 || tmp_f[9] !=  0.0 ||
	tmp_f[10] != 11.0 || tmp_f[11] != 7.0
       );
  gsl_test(s, "gsl_blas_raw_chpr A");
  status += s;

  gsl_blas_raw_ccopy(10, matrix_her_4_cpl, 1, tmp_f, 1);
  gsl_blas_raw_chpr(CblasLower, 4, 2.0, vector_4_c, 1, tmp_f);
  s = ( tmp_f[0] != 11.0 || tmp_f[1] !=   0.0 ||
        tmp_f[2] !=  6.0 || tmp_f[3] !=  -2.0 ||
        tmp_f[4] !=  7.0 || tmp_f[5] !=   0.0 ||
	tmp_f[6] !=  4.0 || tmp_f[7] != -13.0 ||
	tmp_f[8] != 11.0 || tmp_f[9] !=  -7.0 ||
	tmp_f[10] != 17.0 || tmp_f[11] !=  0.0
	);
  gsl_test(s, "gsl_blas_raw_chpr B");
  status += s;


  gsl_blas_raw_zcopy(10, matrix_her_4_zpu, 1, tmp_d, 1);
  gsl_blas_raw_zhpr(CblasUpper, 4, 2.0, vector_4_z, 1, tmp_d);
  s = ( tmp_d[0] != 11.0 || tmp_d[1] !=  0.0 ||
        tmp_d[2] !=  6.0 || tmp_d[3] !=  2.0 ||
        tmp_d[4] !=  4.0 || tmp_d[5] != 13.0 ||
	tmp_d[6] != -1.0 || tmp_d[7] != 27.0 ||
	tmp_d[8] !=  7.0 || tmp_d[9] !=  0.0 ||
	tmp_d[10] != 11.0 || tmp_d[11] != 7.0
       );
  gsl_test(s, "gsl_blas_raw_zhpr A");
  status += s;


  return status;
}

int test_her2(void)
{
  int status = 0;
  int s;

  float  tmp_f[32];
  double tmp_d[32];


  gsl_blas_raw_ccopy(4*4, matrix_her_4_c, 1, tmp_f, 1);
  gsl_blas_raw_cher2(CblasUpper, 4, c_2, vector_4_c, 1, vector_4_c, 1, tmp_f, 4);
  s = ( tmp_f[0] != 21.0 || tmp_f[1] !=  0.0 || tmp_f[2] != 12.0 || tmp_f[3] !=  4.0 ||
        tmp_f[4] != 10.0 || tmp_f[5] != 25.0 || tmp_f[6] != -3.0 || tmp_f[7] != 53.0 ||
	tmp_f[10] != 11.0 || tmp_f[11] != 0.0 || tmp_f[12] != 17.0 || tmp_f[13] != 13.0 ||
	tmp_f[14] !=  8.0 || tmp_f[15] != 34.0 ||
	tmp_f[20] != 35.0 || tmp_f[21] !=  0.0 || tmp_f[22] != 60.5 || tmp_f[23] != 39.0 ||
	tmp_f[30] != 137.0 || tmp_f[31] !=  0.0
       );
  gsl_test(s, "gsl_blas_raw_cher2 A");
  status += s;

  gsl_blas_raw_ccopy(4*4, matrix_her_4_c, 1, tmp_f, 1);
  gsl_blas_raw_cher2(CblasLower, 4, c_2, vector_4_c, 1, vector_4_c, 1, tmp_f, 4);
  s = ( tmp_f[0] != 21.0 || tmp_f[1] !=   0.0 ||
        tmp_f[8] != 12.0 || tmp_f[9] !=  -4.0 || tmp_f[10] != 11.0 || tmp_f[11] != 0.0 ||
        tmp_f[16] != 10.0 || tmp_f[17] != -25.0 || tmp_f[18] != 17.0 || tmp_f[19] != -13.0 ||
	tmp_f[20] != 35.0 || tmp_f[21] !=   0.0 ||
	tmp_f[24] != -3.0 || tmp_f[25] != -53.0 || tmp_f[26] != 8.0 || tmp_f[27] != -34.0 ||
	tmp_f[28] != 60.5 || tmp_f[29] != -39.0 ||
	tmp_f[30] != 137.0 || tmp_f[31] !=  0.0
       );
  gsl_test(s, "gsl_blas_raw_cher2 B");
  status += s;


  gsl_blas_raw_zcopy(4*4, matrix_her_4_z, 1, tmp_d, 1);
  gsl_blas_raw_zher2(CblasUpper, 4, z_2, vector_4_z, 1, vector_4_z, 1, tmp_d, 4);
  s = ( tmp_d[0] != 21.0 || tmp_d[1] !=  0.0 || tmp_d[2] != 12.0 || tmp_d[3] !=  4.0 ||
        tmp_d[4] != 10.0 || tmp_d[5] != 25.0 || tmp_d[6] != -3.0 || tmp_d[7] != 53.0 ||
	tmp_d[10] != 11.0 || tmp_d[11] != 0.0 || tmp_d[12] != 17.0 || tmp_d[13] != 13.0 ||
	tmp_d[14] !=  8.0 || tmp_d[15] != 34.0 ||
	tmp_d[20] != 35.0 || tmp_d[21] !=  0.0 || tmp_d[22] != 60.5 || tmp_d[23] != 39.0 ||
	tmp_d[30] != 137.0 || tmp_d[31] !=  0.0
       );
  gsl_test(s, "gsl_blas_raw_zher2 A");
  status += s;


  return status;
}

int test_hpr2(void)
{
  int status = 0;
  int s;

  float  tmp_f[32];
  double tmp_d[32];


  gsl_blas_raw_ccopy(10, matrix_her_4_cpu, 1, tmp_f, 1);
  gsl_blas_raw_chpr2(CblasUpper, 4, c_2, vector_4_c, 1, vector_4_c, 1, tmp_f);
  s = ( tmp_f[0] != 21.0 || tmp_f[1] !=  0.0 || tmp_f[2] != 12.0 || tmp_f[3] !=  4.0 ||
        tmp_f[4] != 10.0 || tmp_f[5] != 25.0 || tmp_f[6] != -3.0 || tmp_f[7] != 53.0 ||
        tmp_f[8] != 11.0 || tmp_f[9] !=  0.0 || tmp_f[10] != 17.0 || tmp_f[11] != 13.0 ||
        tmp_f[12] != 8.0 || tmp_f[13] != 34.0 ||
        tmp_f[14] != 35.0 || tmp_f[15] != 0.0
       );
  gsl_test(s, "gsl_blas_raw_chpr2 A");
  status += s;

  gsl_blas_raw_ccopy(10, matrix_her_4_cpl, 1, tmp_f, 1);
  gsl_blas_raw_chpr2(CblasLower, 4, c_2, vector_4_c, 1, vector_4_c, 1, tmp_f);
  s = ( tmp_f[0] !=  21.0 || tmp_f[1] !=   0.0 ||
        tmp_f[2] !=  12.0 || tmp_f[3] !=  -4.0 ||
	tmp_f[4] !=  11.0 || tmp_f[5] !=   0.0 ||
	tmp_f[6] !=  10.0 || tmp_f[7] != -25.0 ||
	tmp_f[8] !=  17.0 || tmp_f[9] != -13.0 ||
	tmp_f[10] != 35.0 || tmp_f[11] !=   0.0 ||
	tmp_f[12] != -3.0 || tmp_f[13] != -53.0 ||
	tmp_f[14] !=  8.0 || tmp_f[15] != -34.0
       );
  gsl_test(s, "gsl_blas_raw_chpr2 B");
  status += s;


  gsl_blas_raw_zcopy(10, matrix_her_4_zpu, 1, tmp_d, 1);
  gsl_blas_raw_zhpr2(CblasUpper, 4, z_2, vector_4_z, 1, vector_4_z, 1, tmp_d);
  s = ( tmp_d[0] != 21.0 || tmp_d[1] !=  0.0 || tmp_d[2] != 12.0 || tmp_d[3] !=  4.0 ||
        tmp_d[4] != 10.0 || tmp_d[5] != 25.0 || tmp_d[6] != -3.0 || tmp_d[7] != 53.0 ||
        tmp_d[8] != 11.0 || tmp_d[9] !=  0.0 || tmp_d[10] != 17.0 || tmp_d[11] != 13.0 ||
        tmp_d[12] != 8.0 || tmp_d[13] != 34.0 ||
        tmp_d[14] != 35.0 || tmp_d[15] != 0.0
       );
  gsl_test(s, "gsl_blas_raw_zhpr2 A");
  status += s;


  return status;
}



int main()
{
  int status = 0;

  status += test_L1();
  status += test_gemv();
  status += test_gbmv();
  status += test_trmv();
  status += test_tbmv();
  status += test_tpmv();

  status += test_trsv();
  status += test_tbsv();
  status += test_tpsv();

  status += test_symv();
  status += test_sbmv();
  status += test_spmv();

  status += test_ger();
  status += test_syr();
  status += test_spr();
  status += test_syr2();
  status += test_spr2();

  status += test_hemv();
  status += test_hbmv();
  status += test_hpmv();

  status += test_geru();
  status += test_gerc();

  status += test_her();
  status += test_hpr();
  status += test_her2();
  status += test_hpr2();

  return status;
}
