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
  gsl_blas_raw_sswap(4, tmp_a_d, 1, tmp_b_d, 1);
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
  gsl_test(s, "gsl_blas_raw_scopy");
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
  gsl_test(s, "gsl_blas_raw_saxpy");
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


  /* FIXME: tests for rot stuff */


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


  gsl_blas_raw_dcopy(4, vector_4_d, 1, tmp_d, 1);
  gsl_blas_raw_dgemv(CblasNoTrans, 4, 4, 2.0, matrix_gen_4_d, 4, vector_4_d, 1, 3.0, tmp_d, 1);
  s = ( tmp_d[0] != -4.0 || tmp_d[1] != -5.0 || tmp_d[2] != -1.0 || tmp_d[3] != 7.0 );
  gsl_test(s, "gsl_blas_raw_dgemv A");
  status += s;

  return status;
}

int test_gbmv(void)
{
  int status = 0;
  int s;

/* FIXME */

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


  gsl_blas_raw_dcopy(4, vector_4_d, 1, tmp_d, 1);
  gsl_blas_raw_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit, 4, matrix_gen_4_d, 4, tmp_d, 1);
  s = ( tmp_d[0] != 1.0 || tmp_d[1] != 0.0 || tmp_d[2] != 1.5 || tmp_d[3] != 0.0 );
  gsl_test(s, "gsl_blas_raw_dtrmv A");
  status += s;


  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctrmv(CblasUpper, CblasNoTrans, CblasNonUnit, 4, matrix_gen_4_c, 4, tmp_f, 1);
  s = ( tmp_f[0] != -7.0 || tmp_f[1] != 2.0 || tmp_f[2] != -17.0 || tmp_f[3] != 23.0 );
  gsl_test(s, "gsl_blas_raw_ctrmv A");
  status += s;
/*
  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctrmv(CblasUpper, CblasNoTrans, CblasUnit, 4, matrix_gen_4_c, 4, tmp_f, 1);
  s = ( tmp_f[0] != 1.0 || tmp_f[1] != 2.0 || tmp_f[2] != 1.5 || tmp_f[3] != 3.0 );
  gsl_test(s, "gsl_blas_raw_ctrmv B");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctrmv(CblasUpper, CblasTrans, CblasNonUnit, 4, matrix_gen_4_c, 4, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] != -3.0 || tmp_f[2] != -1.0 || tmp_f[3] != -3.0 );
  gsl_test(s, "gsl_blas_raw_ctrmv C");
  status += s;

  gsl_blas_raw_ccopy(4, vector_4_c, 1, tmp_f, 1);
  gsl_blas_raw_ctrmv(CblasUpper, CblasTrans, CblasUnit, 4, matrix_gen_4_c, 4, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] != -1.0 || tmp_f[2] != -1.0 || tmp_f[3] != 0.0 );
  gsl_test(s, "gsl_blas_raw_ctrmv D");
  status += s;
*/

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
  s = ( tmp_f[0] != -2.0 || tmp_f[1] != 0.0 || tmp_f[2] != -1.0 || tmp_f[3] != -1.0 );
  gsl_test(s, "gsl_blas_raw_stbmv C");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_stbmv(CblasUpper, CblasTrans, CblasUnit, 4, 2, matrix_gen_4_f, 4, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] != 1.0 || tmp_f[2] != -1.0 || tmp_f[3] != 2.0 );
  gsl_test(s, "gsl_blas_raw_stbmv D");
  status += s;

/* FIXME */

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

/*
  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_stpmv(CblasUpper, CblasTrans, CblasNonUnit, 4, matrix_gen_4_f, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] != -0.5 || tmp_f[2] != 1.0 || tmp_f[3] != -13.0 );
  gsl_test(s, "gsl_blas_raw_stpmv C");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_stpmv(CblasUpper, CblasTrans, CblasUnit, 4, matrix_gen_4_f, tmp_f, 1);
  s = ( tmp_f[0] != -2.0 || tmp_f[1] != -1.0 || tmp_f[2] != 1.0 || tmp_f[3] != -4.0 );
  gsl_test(s, "gsl_blas_raw_stpmv D");
  status += s;
*/

/* FIXME */

  return status;
}


/* FIXME: tests for trsv, tbsv, tpsv */


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

  return status;
}

int test_spmv(void)
{
  int status = 0;
  int s;

  float  tmp_f[32];
  double tmp_d[32];

/*
  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_sspmv(CblasUpper, 4, 2.0, matrix_gen_4_f, vector_4_f, 1, 3.0, tmp_f, 1);
  s = ( tmp_f[0] != -4.0 || tmp_f[1] != 26.0 || tmp_f[2] != 14.0 || tmp_f[3] != -17.0 );
  gsl_test(s, "gsl_blas_raw_sspmv A");
  status += s;

  gsl_blas_raw_scopy(4, vector_4_f, 1, tmp_f, 1);
  gsl_blas_raw_sspmv(CblasLower, 4, 2.0, matrix_gen_4_f, vector_4_f, 1, 3.0, tmp_f, 1);
  s = ( tmp_f[0] != 6.0 || tmp_f[1] != -11.0 || tmp_f[2] != -20.0 || tmp_f[3] != 11.0 );
  gsl_test(s, "gsl_blas_raw_sspmv B");
*/
  return status;
}

int test_ger(void)
{
  int status = 0;
  int s;

/* FIXME */

  return status;
}

int test_syr(void)
{
  int status = 0;
  int s;

/* FIXME */

  return status;
}

int test_spr(void)
{
  int status = 0;
  int s;

/* FIXME */

  return status;
}

int test_syr2(void)
{
  int status = 0;
  int s;

/* FIXME */

  return status;
}

int test_spr2(void)
{
  int status = 0;
  int s;

/* FIXME */

  return status;
}

int test_hemv(void)
{
  int status = 0;
  int s;

/* FIXME */

  return status;
}

int test_hbmv(void)
{
  int status = 0;
  int s;

/* FIXME */

  return status;
}

int test_hpmv(void)
{
  int status = 0;
  int s;

/* FIXME */

  return status;
}

int test_geru(void)
{
  int status = 0;
  int s;

/* FIXME */

  return status;
}

int test_gerc(void)
{
  int status = 0;
  int s;

/* FIXME */

  return status;
}

int test_her(void)
{
  int status = 0;
  int s;

/* FIXME */

  return status;
}

int test_hpr(void)
{
  int status = 0;
  int s;

/* FIXME */

  return status;
}

int test_her2(void)
{
  int status = 0;
  int s;

/* FIXME */

  return status;
}

int test_hpr2(void)
{
  int status = 0;
  int s;

/* FIXME */

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
