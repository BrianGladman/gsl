/* Created: [GJ] Fri Jan 17 19:02:12 EST 1997
 */
#include <math.h>
#include "poly_solve.h"
#include "constants.h"
#include "sorting.h"


int gsl_root_solve_quadratic(double a, double b, double c, double x[])
{
  double disc = b*b - 4.*a*c;
  if(disc >= 0.) {
    if(b == 0.) {
      x[0] = -fabs(0.5 * sqrt(disc) / a);
      x[1] = -x[0];
      return 2;
    }
    else {
      double sgnb  = (b > 0. ? 1.0 : -1.0);
      double temp = -0.5 * (b + sgnb * sqrt(disc));
      x[0] = temp/a;
      x[1] = c/temp;
      if(x[1] < x[0]) {
	temp = x[0];
	x[0] = x[1];
	x[1] = temp;
      }
      return 2;
    }
  }
  else {
    return 0;
  }
}
