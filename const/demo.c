#include <stdio.h>
#include <gsl/gsl_const_mksa.h>

int
main (void)
{
  double c  = GSL_CONST_MKS_SPEED_OF_LIGHT;
  double au = GSL_CONST_MKS_ASTRONOMICAL_UNIT;
  double min = GSL_CONST_MKS_MINUTE;

  double r_earth = 1.00 * au;
  double r_mars  = 1.52 * au;

  double tmin, tmax;

  tmin = (r_mars - r_earth)/c;
  tmax = (r_mars + r_earth)/c;

  printf("light travel time from Earth to Mars:\n");
  printf("minimum = %.1f minutes\n", tmin/min);
  printf("maximum = %.1f minutes\n", tmax/min);
}

