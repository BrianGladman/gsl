/* Author: G. Jungman
 * RCS:    $Id$
 */
#include "gsl_qrng.h"


int test_sobol(void)
{
  double v[3];
  int i;

  gsl_qrng * g = gsl_qrng_alloc(gsl_qrng_sobol, 3);

  for(i=0; i<16*1024; i++) {
    gsl_qrng_get(g, v);
    /* printf("%g %g %g\n", v[0], v[1], v[2]); */
  }

  gsl_qrng_free(g);

  return 0;
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
