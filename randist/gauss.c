

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gsl_ran.h"
#include "gsl_randist.h"

double 
gsl_ran_gaussian (void)
{
  double x, y, rr;

  do
    {
      /* choose x,y in uniform square */
      x = -1.0 + 2.0 * gsl_ran_uniform ();
      y = -1.0 + 2.0 * gsl_ran_uniform ();
      rr = x * x + y * y;
      /* see if it is in unit circle */
    }
  while (rr > 1.0 && !(x == 0 && y == 0));

  rr = sqrt (-2.0 * log (rr) / rr);	/* Box-Muller transform */
  /* Ignore one of the random deviates: x*rr;

   * This is a factor of two inefficiency, and if anybody wants to
   * do anything about it, they are welcome to.  In fact, I
   * originally _did_ do something about it, but the code was _very_
   * ugly!
   */
  return y * rr;		/* return the other */
}
