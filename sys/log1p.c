#include <config.h>
#include <math.h>

double gsl_log1p (const double x);

double gsl_log1p (const double x)
{
  volatile double y;
  y = 1 + x;
  return log(y) - ((y-1)-x)/y ;  /* cancels errors with IEEE arithmetic */
}
