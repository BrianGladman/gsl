#include <gsl/gsl_test.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_min.h>


#include "test_funcs.h"

/* stopping parameters for line search */

const double EPSABS_LINE = 1e-2 ;
const double EPSREL_LINE = 1e-2 ;
/*const double EPSABS_LINE = GSL_SQRT_DBL_EPSILON ;
  const double EPSREL_LINE = GSL_SQRT_DBL_EPSILON ;*/
const double EPSABS = GSL_DBL_EPSILON ;

const unsigned int MAX_ITERATIONS_LINE = 100;
const unsigned int MAX_ITERATIONS = 100000;


int
test_fdf(const char * desc, gsl_multimin_function_fdf *f,
	 initpt_function initpt,
	 const gsl_multimin_fdf_minimizer_type *T,
	 int restarting_period)
{
  size_t iterations = 0;
  size_t iterations_line = 0;
  size_t total_i_line = 0;
  int status;
  int just_started = 1;
  double minimum,a,b;
  gsl_interval bracket;
  
  gsl_vector *x = gsl_vector_alloc (f->n);

  gsl_multimin_fdf_minimizer *s;

  (*initpt) (x);

  s = gsl_multimin_fdf_minimizer_alloc(T,f,x,
				       gsl_min_find_bracket,
				       gsl_min_fminimizer_brent);
#ifdef DEBUG
  printf("x "); gsl_vector_fprintf (stdout, s->history->x, "%g"); 
  printf("g "); gsl_vector_fprintf (stdout, s->history->g, "%g"); 
#endif
  do 
    {
      iterations++;
      status = gsl_multimin_fdf_minimizer_next_direction(s);
      status = gsl_multimin_fdf_minimizer_bracket(s,10.0,50);
      if (status == GSL_FAILURE) 
	{
	  if (just_started)
	    {
	      GSL_ERROR ("Can't find bracketing interval after restarting", GSL_FAILURE);
	    }
	  else
	    {
#ifdef DEBUG
	      printf("%i: automatic restart\n",iterations);
#endif
	      gsl_multimin_fdf_minimizer_restart(s);
	      just_started = 1;
	      status = GSL_CONTINUE;
	      continue;
	    }
	}
      else 
	{
	  iterations_line = 0;
	  do 
	    {
	      iterations_line++;
	      status = gsl_multimin_fdf_minimizer_iterate(s);
	      
	      minimum = gsl_min_fminimizer_minimum(s->line_search);
	      bracket = gsl_min_fminimizer_interval(s->line_search);
	      
	      a = bracket.lower;
	      b = bracket.upper;
	      
#ifdef DEBUG
	      printf("%.12f %.18f %.12f %.18f %.12f %.18f\n", 
		     a, s->line_search->f_lower, minimum,s->line_search->f_minimum, b,s->line_search->f_upper);
#endif
	      status = gsl_min_test_interval (bracket, EPSABS_LINE, EPSREL_LINE);
	    }
	  while (status == GSL_CONTINUE && iterations_line < MAX_ITERATIONS_LINE);
	  total_i_line += iterations_line;
	  gsl_multimin_fdf_minimizer_best_step(s);
	}

#ifdef DEBUG
      printf("%i: \n",iterations);
      printf("x "); gsl_vector_fprintf (stdout, s->history->x, "%g"); 
      printf("old x "); gsl_vector_fprintf (stdout, s->history->x1, "%g"); 
      printf("g "); gsl_vector_fprintf (stdout, s->history->g, "%g"); 
      printf("old g "); gsl_vector_fprintf (stdout, s->history->g1, "%g"); 
      printf("f(x) %g\n",s->line_search->f_minimum);
      printf("\n");
#endif
      /* This is not mandatory */
      if (iterations%restarting_period == 0)
	{
	  gsl_multimin_fdf_minimizer_restart(s);
	  just_started = 1;
	}
      else
	{
	  just_started = 0;
	}
      status = gsl_multimin_test_gradient_sqr_norm(s->history,EPSABS);
    }
  while (iterations <= MAX_ITERATIONS && status == GSL_CONTINUE);
  gsl_test(status,
	   "%s (restarting every %d), on %s: %i iterations (%d), f(x)=%g",
	   T->name,restarting_period,desc,iterations,total_i_line,
	   s->line_search->f_minimum);
  gsl_multimin_fdf_minimizer_free(s);
  gsl_vector_free(x);
}

int
main (void)
{
  const gsl_multimin_fdf_minimizer_type *fdfminimizers[5];
  const gsl_multimin_fdf_minimizer_type ** T;

  fdfminimizers[0] = gsl_multimin_fdf_minimizer_steepest_descent;
  fdfminimizers[1] = gsl_multimin_fdf_minimizer_conjugate_pr;
  fdfminimizers[2] = gsl_multimin_fdf_minimizer_conjugate_fr;
  fdfminimizers[3] = gsl_multimin_fdf_minimizer_vector_bfgs;
  fdfminimizers[4] = 0;

  T = fdfminimizers;

  while (*T != 0) 
    {
      test_fdf("Roth", &roth, roth_initpt,*T,5*roth.n);
      test_fdf("Wood", &wood, wood_initpt,*T,5*wood.n);
      test_fdf("Rosenbrock", &rosenbrock, rosenbrock_initpt,*T,5*rosenbrock.n);
      T++;
    }
  
}
