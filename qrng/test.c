/* Author: G. Jungman
 * RCS:    $Id$
 */
#include "gsl_qrng.h"


int test_sobol(void)
{
  double v[3];
  int i;

  gsl_qrng * g = gsl_qrng_alloc(gsl_qrng_sobol, 2);

  for(i=0; i<1024; i++) {
    gsl_qrng_get(g, v);
    /* printf("%g %g\n", v[0], v[1]); */
  }

  gsl_qrng_free(g);

  return 0;
}


int main()
{
  int status = 0;

  status += test_sobol(); 

  return status;
}
