#include <stdio.h>
#include "gsl_sf.h"


int main(int argc, char * argv[])
{
  double x;
  double xmin = 0.;
  double xmax = 10.;
  double dx   = 0.1;

  for(x=xmin; x<xmax; x += dx){

    double y = gsl_sf_bessel_J0(x);

    printf("%19.15g  %19.15g\n", x, y);
  }

  exit(0);
}
