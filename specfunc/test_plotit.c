#include <stdio.h>
#include "gsl_sf.h"

extern double gsl_sf_lngamma_test(double);


int main(int argc, char * argv[])
{
  double x;
  double xmin = 0.5;
  double xmax = 50.;
  double dx   = 0.5;

  double nu = 100.;
  int n = 2;

  for(x=xmin; x<xmax; x += dx){

    double xi = 3.;
    double yr, yi;
    printf("%22.17g    %22.17g\n",
	    x,
	    gsl_sf_lngamma(x)
	    );
  }

  exit(0);
}
