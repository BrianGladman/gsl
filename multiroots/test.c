#include <stdio.h>
#include <gsl_vector.h>
#include <gsl_multiroots.h>

#include "test_funcs.h"

int 
main (void)
{
  int status;

  gsl_vector *x = gsl_vector_alloc (2);

  gsl_multiroot_fdfsolver *s;

  gsl_multiroot_function_fdf function;

  function.f = rosenbrock_f;
  function.df = rosenbrock_df;
  function.fdf = rosenbrock_fdf;
  function.n = 2;
  function.params = 0;

  rosenbrock_initpt (x);

  s = gsl_multiroot_fdfsolver_alloc (gsl_multiroot_fdfsolver_newton,
				     &function, x);

  gsl_vector_fprintf (stdout, x, "%g");
  gsl_vector_fprintf (stdout, s->x, "%g");
  gsl_vector_fprintf (stdout, s->f, "%g");

  do
    {
      status = gsl_multiroot_fdfsolver_iterate (s);
      printf("x "); gsl_vector_fprintf (stdout, s->x, "%g");
      printf("f "); gsl_vector_fprintf (stdout, s->f, "%g");
    }
  while (status == GSL_CONTINUE);

  gsl_multiroot_fdfsolver_free (s);

  return 0;
}
