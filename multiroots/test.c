#include <stdio.h>
#include <gsl_vector.h>
#include <gsl_test.h>
#include <gsl_multiroots.h>

#include "test_funcs.h"

int test_fdf (const char * desc, gsl_multiroot_function_fdf * function, initpt_function initpt, const gsl_multiroot_fdfsolver_type * T);


int 
main (void)
{
  const gsl_multiroot_fdfsolver_type * fdfsolvers[2] ;

  const gsl_multiroot_fdfsolver_type ** T ;

  fdfsolvers[0] = gsl_multiroot_fdfsolver_newton;
  fdfsolvers[1] = gsl_multiroot_fdfsolver_mtrnewton;
  fdfsolvers[2] = 0 ;

  T = fdfsolvers ;

  while (*T != 0) 
    {
      test_fdf ("Rosenbrock", &rosenbrock, rosenbrock_initpt, *T);
      test_fdf ("Powell singular", &powellsing, powellsing_initpt, *T);
      test_fdf ("Powell badly scaled", &powellscal, powellscal_initpt, *T);
      test_fdf ("Wood", &wood, wood_initpt, *T);
      T++;
    }

  return 0;
}

int
test_fdf (const char * desc, gsl_multiroot_function_fdf * function, 
          initpt_function initpt,
          const gsl_multiroot_fdfsolver_type * T)
{
  int status;

  gsl_vector *x = gsl_vector_alloc (function->n);

  gsl_multiroot_fdfsolver *s;

  (*initpt) (x);

  s = gsl_multiroot_fdfsolver_alloc (T, function, x);

  printf("x "); gsl_vector_fprintf (stdout, s->x, "%g"); printf("\n");
  /* printf("f "); gsl_vector_fprintf (stdout, s->f, "%g"); printf("\n");
  printf("J "); gsl_matrix_fprintf (stdout, s->J, "%g"); printf("\n"); */

  do
    {
      status = gsl_multiroot_fdfsolver_iterate (s);
      /* printf("x "); gsl_vector_fprintf (stdout, s->x, "%g"); printf("\n");
      printf("f "); gsl_vector_fprintf (stdout, s->f, "%g"); printf("\n");
      printf("J "); gsl_matrix_fprintf (stdout, s->J, "%g"); printf("\n"); */

      status = gsl_multiroot_test_residual (s->f, 0.0000001);
    }
  while (status == GSL_CONTINUE);

  printf("x "); gsl_vector_fprintf (stdout, s->x, "%g"); printf("\n");
  printf("f "); gsl_vector_fprintf (stdout, s->f, "%g"); printf("\n");

  gsl_multiroot_fdfsolver_free (s);
  gsl_vector_free(x);

  gsl_test(status, "%s with %s", desc, T->name);

  return status;
}
