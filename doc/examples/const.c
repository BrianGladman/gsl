#include <stdio.h>
#include <gsl/gsl_const_mks.h>

int
main (void)
{
  double c  = GSL_CONST_MKS_SPEED_OF_LIGHT;
  double au = GSL_CONST_MKS_ASTRONOMICAL_UNIT;
  double minutes = GSL_CONST_MKS_MINUTE;

  /* distance stored in meters */
  double r_earth = 1.00 * au;  
  double r_mars  = 1.52 * au;

  double t_min, t_max;

  t_min = (r_mars - r_earth) / c;
  t_max = (r_mars + r_earth) / c;

  printf ("light travel time from Earth to Mars:\n");
  printf ("minimum = %.1f minutes\n", t_min / minutes);
  printf ("maximum = %.1f minutes\n", t_max / minutes);

  return 0;
}
