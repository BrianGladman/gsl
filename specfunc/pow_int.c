#include <math.h>
#include <gsl_math.h>
#include "gsl_sf_pow_int.h"


double gsl_sf_pow_int(double x, int n)
{
  double value = 1;

  if(abs(n) > 50) return pow(x, n);  /* Defer for large powers. */
  
  /* Trap. */
  if(x == 0) {
    if(n >= 0) {
      return 0.;
    }
    else {
      return GSL_INF;
    }
  }
  else {
    if(n == 0) return 1.;
  } 

  if(n < 0) {
    x = 1./x;
    n = -n;
  }

  for(; n > 0; n--){
    value *= x;
  }
  return value;
}
