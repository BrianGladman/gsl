#include <stdio.h>
#include <gsl_errno.h>
#include <gsl_math.h>
#include <gsl_roots.h>

/* Compile with:

   gcc -I. -I.. -I../err demo.c libgslroots.a ../err/libgslerr.a  -lm 

 */


struct quadratic_params
  {
    double a, b, c;
  };

double quadratic (double x, void *params);
double quadratic_deriv (double x, void *params);
void quadratic_fdf (double x, void *params, double *y, double *dy);

int
main ()
{
  int status;
  int iterations = 0, max_iterations = 100;
  gsl_root_fsolver *s;
  gsl_interval x = {0.0, 5.0};
  gsl_function F;
  struct quadratic_params params = {1.0, 0.0, -5.0};

  F.function = &quadratic;
  F.params = &params;

  s = gsl_root_fsolver_alloc (gsl_root_fsolver_brent, &F, x);

  do
    {
      iterations++;
      gsl_root_fsolver_iterate (s);
      status = gsl_root_test_interval (s->interval, 0.001, 0.001);
      printf ("%5d  %.7f [%.7f,%.7f]\n",
	      iterations, s->root, s->interval.lower, s->interval.upper);
    }
  while (status == GSL_CONTINUE && iterations < max_iterations);

  printf ("best estimate of root = %.7f\n", s->root);
  printf ("          actual root = %.7f\n", sqrt (5.0));
  printf ("interval bounds = [%.7f,%.7f]\n", 
	  s->interval.lower, s->interval.upper);

}


double
quadratic (double x, void *params)
{
  struct quadratic_params *p = (struct quadratic_params *) params;

  double a = p->a;
  double b = p->b;
  double c = p->c;

  return (a * x + b) * x + c;
}

double
quadratic_deriv (double x, void *params)
{
  struct quadratic_params *p = (struct quadratic_params *) params;

  double a = p->a;
  double b = p->b;
  double c = p->c;

  return 2.0 * a * x + b;
}

void
quadratic_fdf (double x, void *params, double *y, double *dy)
{
  struct quadratic_params *p = (struct quadratic_params *) params;

  double a = p->a;
  double b = p->b;
  double c = p->c;

  *y = (a * x + b) * x + c;
  *dy = 2.0 * a * x + b;
}
