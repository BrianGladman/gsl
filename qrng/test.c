/* Author: G. Jungman
 * RCS:    $Id$
 */
#include "gsl_qrng.h"


int test_sobol(void)
{
  int status = 0;
  double v[3];
  /* int i; */

  /* test in dimension 2 */
  gsl_qrng * g = gsl_qrng_alloc(gsl_qrng_sobol, 2);
  gsl_qrng_get(g, v);
  gsl_qrng_get(g, v);
  gsl_qrng_get(g, v);
  status += ( v[0] != 0.25 || v[1] != 0.75 );
  gsl_qrng_get(g, v);
  status += ( v[0] != 0.375 || v[1] != 0.375 );
  gsl_qrng_free(g);

  /* test in dimension 3 */
  g = gsl_qrng_alloc(gsl_qrng_sobol, 3);
  gsl_qrng_get(g, v);
  gsl_qrng_get(g, v);
  gsl_qrng_get(g, v);
  status += ( v[0] != 0.25 || v[1] != 0.75 || v[2] != 0.25 );
  gsl_qrng_get(g, v);
  status += ( v[0] != 0.375 || v[1] != 0.375 || v[2] != 0.625 );
  gsl_qrng_free(g);

  return status;
}


int test_nied2(void)
{
  double v[3];
  int i;

  gsl_qrng * g = gsl_qrng_alloc(gsl_qrng_niederreiter_2, 3);

  for(i=0; i<16*1024; i++) {
    gsl_qrng_get(g, v);
    printf("%g %g %g\n", v[0], v[1], v[2]);
  }

  gsl_qrng_free(g);

  return 0;
}


int main()
{
  int status = 0;

  status += test_sobol(); 
  /* status += test_nied2(); */

  return status;
}
