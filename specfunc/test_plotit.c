#include <stdio.h>
#include "gsl_sf.h"

extern double gsl_sf_lngamma_test(double);
extern int gsl_sf_conical_sph_reg_0_impl(double, double, double, double*);
extern int gsl_sf_conical_sph_reg_array_impl(int, double, double, double *, double *);

extern int gsl_sf_conical_sph_reg_array_impl(int, double, double, double *, double *);


int main(int argc, char * argv[])
{
  double x;
  double xmin =  -10.0;
  double xmax =   10.0;
  double dx   = 0.1;
  double lambda = 0.5;
  double y;
  double y_array[1024];
  int i;

  double nu = 100.;
  int n = 2;

  for(x=xmin; x<xmax; x += dx) {
    y = gsl_sf_Shi(x);
    printf("%22.17g    %22.17g\n",
	    x,
	    y
	    );
  }

/*
  for(x=xmin; x<xmax; x += dx){
    gsl_sf_conical_sph_reg_0_impl(lambda, 1-x, 1+x, &y);
    printf("%22.17g    %22.17g\n",
	    x,
	    y
	    );
  }
*/

/*
  gsl_sf_conical_sph_reg_array_impl(100, lambda, 0.5, &y, y_array);
  
  for(i=0; i<=100; i++) {
  printf("%3d   %22.17g    %22.17g\n",
	    i,
	    x,
	    y_array[i]
	    );
  }
*/

  
  exit(0);
}
