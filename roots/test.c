/* test -- program to test root finding functions */


/* config headers */
#include <config.h>

/* standard headers */
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/* gsl headers */
#include <gsl_errno.h>

/* gsl/roots headers */
#include <gsl_roots.h>
#include "test.h"


/* administer all the tests */
int
main (void)
{
  int num_tests = 0, num_passed = 0, num_failed = 0;

  printf( "Testing root finding functions (GSL version %s)...\n", VERSION );

  /* turn off error handling -- we will take care of everything */
  gsl_set_error_handler (gsl_no_error_handler);

  /* run tests */
  test_bisection (&num_tests, &num_passed, &num_failed); 
  /* etc. */

  /* now summarize the results */
  printf( "Testing complete.\n" );
  if (num_tests == num_passed)
    printf( "All tests passed.\n" );
  else {
    printf( "%d of %d tests passed, %d failed.\n", num_passed, num_tests,
            num_failed );
    printf( "The integrity of this package may be suspect!\n" );
  }
}

void
test_bisection (int * num_tests, int * num_passed, int * num_failed)
{
    printf( "Testing bisection: " );
}
