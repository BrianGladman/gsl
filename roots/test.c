#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_roots.h>

#include <gsl/gsl_ieee_utils.h>

#include "test.h"

int
main (void)
{
  gsl_ieee_env_setup();

  test_macros();
  test_roots();

  /* now summarize the results */

  return gsl_test_summary ();
}

