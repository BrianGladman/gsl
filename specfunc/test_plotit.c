#include <stdio.h>
#include "gsl_sf.h"


int main(int argc, char * argv[])
{
  double x;
  double xmin = 0.;
  double xmax = 10.;
  double dx   = 0.1;

  double nu = 100.;
  int n = 2;

  for(x=xmin; x<xmax; x += dx){

    double y;
    int status = gsl_sf_bessel_Jnu_taylor_e(nu, x, 10, &y);

    printf("%20.15g  %20.15g", x, y);
    if(status) printf("     *");
    printf("\n");
  }

  exit(0);
}
