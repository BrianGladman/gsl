#include <math.h>
#include "gsl_specfunc.h"


double pow_int(double x, int n)
{
  double value = 1;

  if(abs(n) > 50) return pow(x, n);  /* Defer for large powers. */
  if(x == 0 || n == 0) return 1;     /* Trap. */

  if(n < 0) {
    x = 1./x;
    n = -n;
  }

  for(; n > 0; n--){
    value *= x;
  }
  return value;
}
