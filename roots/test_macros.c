#include <gsl_test.h>
#include <gsl_math.h>
#include <gsl_roots.h>

#include "test.h"

/* Test certain macros. */
void
test_macros (void)
{
  int result;
  double inf, nan ;

  /* 1.0 is real */
  result = GSL_IS_REAL (1.0);
  gsl_test (result != 1, "GSL_IS_REAL(1.0) is 1");

  inf = 1.0 / (sqrt(1.0) - 1) ;

  /* 1.0/0.0 == Inf is not real */
  result = GSL_IS_REAL (inf);
  gsl_test (result != 0, "GSL_IS_REAL(Inf) is 0");

  nan = inf - inf ;

  /* 0.0/0.0 == NaN is not real */
  result = GSL_IS_REAL (nan);
  gsl_test (result != 0, "GSL_IS_REAL(NaN) is 0");
}
