#include <stdio.h>
#include <gsl_vector.h>
#include <gsl_test.h>
#include <gsl_multiroots.h>

#include "test_funcs.h"
int test_fdf (const char * desc, gsl_multiroot_function_fdf * function, initpt_function initpt, const gsl_multiroot_fdfsolver_type * T);
int test_f (const char * desc, gsl_multiroot_function_fdf * fdf, initpt_function initpt, const gsl_multiroot_fsolver_type * T);


int 
main (void)
{
  const gsl_multiroot_fsolver_type * fsolvers[3] ;
  const gsl_multiroot_fsolver_type ** T1 ;

  const gsl_multiroot_fdfsolver_type * fdfsolvers[3] ;
  const gsl_multiroot_fdfsolver_type ** T2 ;

  fsolvers[0] = gsl_multiroot_fsolver_dnewton;
  fsolvers[1] = gsl_multiroot_fsolver_broyden;
  fsolvers[2] = 0;

  fdfsolvers[0] = gsl_multiroot_fdfsolver_newton;
  fdfsolvers[1] = gsl_multiroot_fdfsolver_gnewton;
  fdfsolvers[2] = 0 ;

  T1 = fsolvers ;
  T2 = fdfsolvers ;

  while (*T1 != 0) 
    {
      test_f ("Rosenbrock", &rosenbrock, rosenbrock_initpt, *T1); 
      test_f ("Powell singular", &powellsing, powellsing_initpt, *T1);  
      test_f ("Powell badly scaled", &powellscal, powellscal_initpt, *T1); 
      test_f ("Wood", &wood, wood_initpt, *T1); 
      test_f ("Helical", &helical, helical_initpt, *T1); 
      T1++;
    }

  while (*T2 != 0) 
    {
      test_fdf ("Rosenbrock", &rosenbrock, rosenbrock_initpt, *T2);
      test_fdf ("Powell singular", &powellsing, powellsing_initpt, *T2);
      test_fdf ("Powell badly scaled", &powellscal, powellscal_initpt, *T2);
      test_fdf ("Wood", &wood, wood_initpt, *T2);
      test_fdf ("Helical", &helical, helical_initpt, *T2);
      T2++;
    }

  return 0;
}

int
test_fdf (const char * desc, gsl_multiroot_function_fdf * function, 
          initpt_function initpt,
          const gsl_multiroot_fdfsolver_type * T)
{
  int status;
  double residual = 0;
  size_t i, n = function->n, iter = 0;
  
  gsl_vector *x = gsl_vector_alloc (n);

  gsl_multiroot_fdfsolver *s;

  (*initpt) (x);

  s = gsl_multiroot_fdfsolver_alloc (T, function, x);

  do
    {
      iter++;
      status = gsl_multiroot_fdfsolver_iterate (s);
      status = gsl_multiroot_test_residual (s->f, 0.0000001);
    }
  while (status == GSL_CONTINUE);

#ifdef DEBUG
  printf("x "); gsl_vector_fprintf (stdout, s->x, "%g"); printf("\n");
  printf("f "); gsl_vector_fprintf (stdout, s->f, "%g"); printf("\n");
#endif

  for (i = 0; i < n ; i++)
    {
      residual += fabs(gsl_vector_get(s->f, i));
    }

  gsl_multiroot_fdfsolver_free (s);
  gsl_vector_free(x);

  gsl_test(status, "%s with %s, %u iterations, residual = %g", desc, T->name, iter, residual);

  return status;
}


int
test_f (const char * desc, gsl_multiroot_function_fdf * fdf, 
        initpt_function initpt,
        const gsl_multiroot_fsolver_type * T)
{
  int status;
  size_t i, n = fdf->n, iter = 0;
  double residual = 0;

  gsl_vector *x;

  gsl_multiroot_fsolver *s;
  gsl_multiroot_function function;

  function.f = fdf->f;
  function.params = fdf->params;
  function.n = n ;

  x = gsl_vector_alloc (n);

  (*initpt) (x);

  s = gsl_multiroot_fsolver_alloc (T, &function, x);

/*   printf("x "); gsl_vector_fprintf (stdout, s->x, "%g"); printf("\n"); */
/*   printf("f "); gsl_vector_fprintf (stdout, s->f, "%g"); printf("\n"); */

  do
    {
      iter++;
      status = gsl_multiroot_fsolver_iterate (s);
      status = gsl_multiroot_test_residual (s->f, 0.0000001);

/*     printf("x "); gsl_vector_fprintf (stdout, s->x, "%g"); printf("\n");  */
/*     printf("f "); gsl_vector_fprintf (stdout, s->f, "%g"); printf("\n");   */

    }
  while (status == GSL_CONTINUE);

#ifdef DEBUG
  printf("x "); gsl_vector_fprintf (stdout, s->x, "%g"); printf("\n");
  printf("f "); gsl_vector_fprintf (stdout, s->f, "%g"); printf("\n");
#endif

  for (i = 0; i < n ; i++)
    {
      residual += fabs(gsl_vector_get(s->f, i));
    }

  gsl_multiroot_fsolver_free (s);
  gsl_vector_free(x);

  gsl_test(status, "%s with %s, %u iterations, residual = %g", desc, T->name, iter, residual);

  return status;
}
