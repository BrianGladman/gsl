/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl_test.h>
#include <gsl_errno.h>
#include <gsl_interp.h>


double test_func(double x)
{
  return sin(x);
}


int main(int argc, char ** argv)
{
  int status;
  int N = 20000;
  int i;

  double x;
  double xmin = 0.0;
  double xmax = 2.0*M_PI;
  double dx = 5.e-6 * xmax;

  double * xa = (double *) malloc(N * sizeof(double));
  double * ya = (double *) malloc(N * sizeof(double));
  
  gsl_interp_factory f = gsl_interp_factory_cspline_periodic;
  gsl_interp_accel * a = gsl_interp_accel_new();
  gsl_interp_obj * interp;

  for(i=0; i<N; i++) {
    xa[i] = xmin + (xmax - xmin) * (double)i/(double)(N-1);
    ya[i] = test_func(xa[i]);
  }

  interp = f.create(xa, ya, N);

  for(x=xmin; x<xmax; x += dx) {
    double y_func = test_func(x);
    double y;
    int eval_stat = gsl_interp_eval_impl(interp, xa, ya, x, a, &y);
    
    printf("%g   %g    %g    %g    %s\n",
                x,
		y_func,
		y,
		(y - y_func)/(y + y_func),
		gsl_strerror(eval_stat)
		);
    
  }

  return gsl_test_summary();
}
