#include <stdio.h>
#include <gsl_math.h>
#include "gsl_sf.h"

#include "bessel.h"
#include "bessel_Jn_impl.h"

int main(int argc, char * argv[])
{
  double x;
  double xmin =  -10.0;
  double xmax =   10.0;
  double dx   = 0.1;
  double lambda = 0.5;
  double y;
  double y_array[4096];
  int i;

  double nu = 100.;
  int n = 2;

  for(n=1; n<=1000; n++) {
    x = 50.0 * n;
    gsl_sf_bessel_Inu_scaled_asymp_unif_impl((double)n, x, &y);
    printf("%3d    %22.17g   %22.17g    %22.17g\n",
	    n,
	    x,
	    y,
	    (n - pow(n, 0.34))/n
	    );
  }


/*
  for(x=xmin; x<xmax; x += dx) {
    y = gsl_sf_Shi(x);
    printf("%22.17g    %22.17g\n",
	    x,
	    y
	    );
  }
*/

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
  gsl_sf_conical_sph_reg_array_impl(60, lambda, 0.51324, &y, y_array);
  
  for(i=0; i<=60; i++) {
  printf("%3d   %22.17g    %22.17g\n",
	    i,
	    x,
	    y_array[i]
	    );
  }
*/

/*
  gsl_sf_conical_sph_reg_array_e(100, 1., 50., y_array);
  
  for(i=0; i<=100; i++) {
  printf("%3d   %22.17g    %22.17g\n",
	    i,
	    x,
	    y_array[i]
	    );
  }
*/

/*
  gsl_sf_bessel_In_scaled_array_e(100, 1., y_array);
  
  for(i=0; i<=100; i++) {
  printf("%3d   %22.17g    %22.17g\n",
	    i,
	    x,
	    y_array[i]
	    );
  }
*/

/*
  bessel_In_scaled(10, 3., &y, y_array);
  printf("%d  %22.17g     %22.17g   %22.17g\n", 10, 3., y, y_array[0]);
*/
  exit(0);
}
