#include <stdio.h>
#include "gsl_sf.h"

extern double gsl_sf_lngamma_test(double);
extern int gsl_sf_legendre_sph_con_reg_0_impl(double, double, double, double*);
extern int gsl_sf_legendre_sph_con_reg_array_impl(int, double, double, double *, double *);

extern int gsl_sf_hyper_array_impl(int, double, double, double *, double *);


int main(int argc, char * argv[])
{
  double x;
  double xmin = -2.0;
  double xmax =  2.0;
  double dx   = 0.1;
  double lambda = 0.5;
  double y;
  double y_array[1024];
  int i;

  double nu = 100.;
  int n = 2;

  
  /*
  for(x=xmin; x<xmax; x += dx){
    gsl_sf_legendre_sph_con_reg_0_impl(lambda, x, &y);
    printf("%22.17g    %22.17g\n",
	    x,
	    y
	    );
  }
  */


  gsl_sf_hyper_array_impl(1000, lambda, 2., &y, y_array);
  
  for(i=0; i<=1000; i++) {
  printf("%3d   %22.17g    %22.17g\n",
	    i,
	    x,
	    y_array[i]
	    );
  }

  exit(0);
}
