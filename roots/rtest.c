/* tests the root finding library in GSL */

#include <stdio.h>
#include <math.h>

#include <roots.h>

int main()
{
  double epsilon;

  printf("testing newton's methtod\n");
  for (epsilon = 0.1; epsilon > 1.0e-8; epsilon /= 10) {
    gsl_set_newton_epsilon(epsilon);
    printf("    fn=sin, dfn=cos, guess=2.3, epsilon=%g, root is %.20g\n",
	   gsl_get_newton_epsilon(), gsl_newton1D(sin, cos, 2.3));
  }
  return 0;
}
