#include <gsl_math.h>
#include <gsl_test.h>

#include "test.h"

int
main (void)
{
  test_poly();

  /* now summarize the results */

  return gsl_test_summary ();
}


