#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#define N 8

int
main (void)
{
  int i;
  double xi, yi;
  double x[N] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 };
  double y[N] = { 0.0, 0.0, 2.0, 0.0,-1.0,-2.0,-1.0, 0.0 };

  printf ("#m=0,S=2\n");

  for (i = 0; i < N; i++)
    {
      printf ("%g %g\n", x[i], y[i]);
    }

  printf ("#m=1,S=0\n");

  {
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline = gsl_spline_alloc (gsl_interp_akima, N);
    gsl_spline_init (spline, x, y, N);

    for (xi = x[0]; xi < x[N-1]; xi += 0.01)
      {
        double yi = gsl_spline_eval (spline, xi, acc);
        printf ("%g %g\n", xi, yi);
      }
    gsl_spline_free (spline);
    gsl_interp_accel_free(acc);
  }
}
