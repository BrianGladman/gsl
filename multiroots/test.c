#include <stdio.h>
#include <gsl_vector.h>
#include <gsl_multiroots.h>

#include "test_funcs.h"

int 
main (void)
{
  int status;

  gsl_vector *x = gsl_vector_alloc (4);

  gsl_multiroot_fdfsolver *s;

  gsl_multiroot_function_fdf function;

#ifdef JUNK
  function.f = rosenbrock_f;
  function.df = rosenbrock_df;
  function.fdf = rosenbrock_fdf;
  function.n = 2;
  function.params = 0;

  rosenbrock_initpt (x);
#endif

  function.f = powellsingular_f;
  function.df = powellsingular_df;
  function.fdf = powellsingular_fdf;
  function.n = 4;
  function.params = 0;

  powellsingular_initpt (x);


  s = gsl_multiroot_fdfsolver_alloc (gsl_multiroot_fdfsolver_newton,
				     &function, x);

  gsl_vector_fprintf (stdout, x, "%g");
  printf("x "); gsl_vector_fprintf (stdout, s->x, "%g"); printf("\n");
  printf("f "); gsl_vector_fprintf (stdout, s->f, "%g"); printf("\n");
  printf("J "); gsl_matrix_fprintf (stdout, s->J, "%g"); printf("\n");

  do
    {
      status = gsl_multiroot_fdfsolver_iterate (s);
      printf("x "); gsl_vector_fprintf (stdout, s->x, "%g"); printf("\n");
      printf("f "); gsl_vector_fprintf (stdout, s->f, "%g"); printf("\n");
      printf("J "); gsl_matrix_fprintf (stdout, s->J, "%g"); printf("\n");

      status = gsl_multiroot_test_residual (s->f, 0.0001);
    }
  while (status == GSL_CONTINUE);

  gsl_multiroot_fdfsolver_free (s);
  gsl_vector_free(x);

  return 0;
}
