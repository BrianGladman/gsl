#include <gsl_math.h>
#include <gsl_test.h>
#include <gsl_roots.h>

#include "test.h"

int
main (void)
{
  test_macros();
  test_roots();

  /* now summarize the results */

  return gsl_test_summary ();
}

