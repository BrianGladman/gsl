#include <math.h>

#include <roots.h>


/* the epsilon used is in this static variable; there are get/set
 * routines to modify it which are available to the outside.
 */
static double newton_epsilon = 0.0001;

void gsl_set_newton_epsilon(double new_val)
{
  if (new_val > 0) {
    newton_epsilon = new_val;
  } else {
    /* FIXME: should give some nasty error message!! */
  }
}

double gsl_get_newton_epsilon()
{
  return newton_epsilon;
}

/* finds the root, given fn, dfn and guess.  uses the current epsilon */
double gsl_newton1D(double (*fn)(double x), double (*dfn)(double x),
		    double guess)
{
  double x, old_x;

  old_x = x = guess;

  while (fabs(fn(x)) > newton_epsilon) {
    x = old_x - fn(old_x)/dfn(old_x);
    old_x = x;
  }

  return x;
}
