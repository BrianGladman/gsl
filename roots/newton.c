#include <math.h>

#include <roots.h>

#define NEWTON1D_EPSILON 0.0001

double gsl_newton1D(double (*fn)(double x), double (*dfn)(double x),
		    double guess)
{
  double x, old_x;

  old_x = x = guess;

  while (fabs(fn(x)) > NEWTON1D_EPSILON) {
    x = old_x - fn(old_x)/dfn(old_x);
    old_x = x;
  }

  return x;
}
