#include <stdio.h>
#include <gsl_math.h>
#include "gsl_sf.h"

#include "airy_impl.h"
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

plotto();
exit(0);

/*
  for(x=-5; x < 10.; x += .05) {
    printf("%20.16g    %20.16g %20.16g     %20.16g %20.16g\n",
           x,
	   gsl_sf_airy_Ai(x),
	   gsl_sf_airy_Bi(x),
	   gsl_sf_airy_Ai_scaled(x),
	   gsl_sf_airy_Bi_scaled(x)
           );
  }
  exit(0);
*/

/*
  struct gsl_sf_ChebSeries * cs = gsl_sf_cheb_new(sin, -M_PI, M_PI, 15);
  
  for(x=-M_PI; x<M_PI; x += 0.1) {
    printf("%20.16g    %20.16g %20.16g     %20.16g %20.16g    %20.16g %20.16g\n",
           x,
	   gsl_sf_cheb_eval(x, cs),
	   sin(x),
	   gsl_sf_cheb_eval_deriv(x, cs),
	   cos(x),
	   gsl_sf_cheb_eval_integ(x, cs),
	   -cos(x)
	   );
  }
  exit(0);
*/
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
