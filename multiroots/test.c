#include <gsl_multiroots.h>

#include "test_funcs.h"

int main (void)
{
  gsl_vector * x = gsl_vector_alloc(2);

  gsl_multiroot_fdfsolver * s ;

  gsl_multiroot_function_fdf function ;

  function.f = rosenbrock_f;
  function.df = rosenbrock_df;
  function.fdf = rosenbrock_fdf;
  function.n = 2;
  function.params = 0 ;

  s = gsl_multiroot_fdfsolver_alloc(gsl_multiroot_fdfsolver_newton, 
                                    &function, x);

  gsl_multiroot_fdfsolver_iterate(s);

  gsl_multiroot_fdfsolver_free(s);

  return 0;
}

