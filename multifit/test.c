/* These tests are based on the NIST Statistical Reference Datasets
   See http://www.nist.gov/itl/div898/strd/index.html for more
   information. */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_multifit.h>

#include <gsl/gsl_ieee_utils.h>

#include "test_longley.c"
#include "test_filip.c"
#include "test_pontius.c"
#include "test_fn.c"

int
main (void)
{
  gsl_ieee_env_setup();

  test_longley();
  test_filip();
  test_pontius();

  test_lmder();

  /* now summarize the results */

  return gsl_test_summary ();
}
