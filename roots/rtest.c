/* tests the root finding library in GSL */

#include <math.h>

#include <roots.h>

int main()
{
  gsl_newton1D(sin, cos, 0.2);

  return 0;
}
