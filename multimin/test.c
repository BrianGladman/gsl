#include <gsl_multimin.h>
#include <gsl_min.h>

/* stopping parameters for line search */

const double EPSABS_LINE = 0.00001 ;
const double EPSREL_LINE = 0.00001 ;
const double EPSABS = GSL_DBL_EPSILON ;

const unsigned int MAX_ITERATIONS_LINE = 100;

f_rosenbrock(const gsl_vector * x, void * params)
{
  double u,v;

  u=gsl_vector_get(x,0);
  v=gsl_vector_get(x,1);
  return (u-1)*(u-1)+10*(u*u-v)*(u*u-v);
}

void
f_rosenbrock_df(const gsl_vector * x, void * params,gsl_vector * df)
{
  double u,v;

  u=gsl_vector_get(x,0);
  v=gsl_vector_get(x,1);
  gsl_vector_set(df,0,2*(u-1)+40*u*(u*u-v));
  gsl_vector_set(df,1,-2*(u*u-v));
}

void
f_rosenbrock_fdf(const gsl_vector * x, void * params,double *f,gsl_vector * df)
{
  f_rosenbrock_df(x,params,df);
  *f=f_rosenbrock(x,params);
}

double
f_wood(const gsl_vector * x, void * params)
{
  double u1 = gsl_vector_get(x,0);
  double u2 = gsl_vector_get(x,1);
  double u3 = gsl_vector_get(x,2);
  double u4 = gsl_vector_get(x,3);

  return 100 * (u1 * u1 - u2) * (u1 * u1 - u2) + (1 - u1) * (1 - u1)
    + 90 * (u3 * u3 - u4) * (u3 * u3 - u4) + (1 - u3) * (1 - u3)
    + 10.1 * ( (1 - u2) * (1 - u2) + (1 - u4) * (1 - u4) )
    + 19.8 * (1 - u2) * (1 - u4);

}

void
f_wood_df(const gsl_vector * x, void * params,gsl_vector * df)
{
  double u1 = gsl_vector_get(x,0);
  double u2 = gsl_vector_get(x,1);
  double u3 = gsl_vector_get(x,2);
  double u4 = gsl_vector_get(x,3);

  gsl_vector_set(df,0, 400 * u1 * (u1 * u1 - u2) - 2 * (1 - u1) );
  gsl_vector_set(df,1, -200 * (u1 * u1 - u2) - 20.2 * (1 - u2) - 19.8 * (1 - u4) );
  gsl_vector_set(df,2, 360 * u3 * (u3 * u3 - u4) - 2 * (1 - u3) );
  gsl_vector_set(df,3, -180 * (u3 * u3 - u4) - 20.2 * (1 - u4) - 19.8 * (1 - u2) );
}

void
f_wood_fdf(const gsl_vector * x, void * params,double *f,gsl_vector * df)
{
  f_wood_df(x,params,df);
  *f=f_wood(x,params);
}

int
main (void)
{
  gsl_multimin_function_fdf f;
  gsl_multimin_fdf_minimizer *s;
  gsl_vector * x;
  size_t iterations = 0;
  size_t iterations_line = 0;
  int status;
  int just_started = 1;
  gsl_min_fminimizer * min_one_dim=NULL;

  double minimum,f_minimum,f_lower,f_upper,a,b;
  gsl_interval bracket;

  f.f=f_wood;
  f.df=f_wood_df;
  f.fdf=f_wood_fdf;
  f.n=4;

  x=gsl_vector_calloc(f.n);
  gsl_vector_set(x,0,-3.0);
  gsl_vector_set(x,1,-1.0);
  gsl_vector_set(x,2,-3.0);
  gsl_vector_set(x,3,-1.0);

  s=gsl_multimin_fdf_minimizer_alloc(gsl_multimin_fdf_minimizer_conjugate_fr,
				     &f,
				     x);
/*  s=gsl_multimin_fdf_minimizer_alloc(gsl_multimin_fdf_minimizer_steepest_descent,
				     &f,
				     x);*/
  do 
    {
      iterations++;
      status = gsl_multimin_fdf_minimizer_next_direction(s);
      bracket.lower = 0;
      bracket.upper = 2;
      f_lower = s->history->f;
      f_upper = GSL_FN_EVAL(s->f_directional,bracket.upper);
      status = gsl_min_find_bracket(s->f_directional,&minimum,&f_minimum,
				    &bracket,&f_lower,&f_upper,10);
      if (status == GSL_FAILURE) 
	{
	  if (f_minimum < f_lower)
	    {
	      gsl_multimin_fdf_minimizer_step_with_value(s,minimum,f_minimum);
	    }
	  else
	    {
	      if (just_started)
		{
		  GSL_ERROR ("Can't find bracketing interval after restarting", GSL_FAILURE);
		}
	      else
		{
		  printf("%i: automatic restart\n",iterations);
		  gsl_multimin_fdf_minimizer_restart(s);
		  just_started = 1;
		  continue;
		}
	    }
	}
      else 
	{
	  if (min_one_dim == NULL) 
	    {
	      min_one_dim 
		= gsl_min_fminimizer_alloc_with_values (gsl_min_fminimizer_brent,
							s->f_directional,
							minimum,f_minimum,
							bracket,f_lower,f_upper);
	    } 
	  else 
	    {
	      gsl_min_fminimizer_set_with_values(min_one_dim,
						 s->f_directional,
						 minimum,f_minimum,
						 bracket,f_lower,f_upper);
	    }
	  iterations_line = 0;
	  do 
	    {
	      iterations_line++;
	      status = gsl_min_fminimizer_iterate (min_one_dim);
	      
	      minimum = gsl_min_fminimizer_minimum(min_one_dim);
	      bracket = gsl_min_fminimizer_interval(min_one_dim);
	      
	      a = bracket.lower;
	      b = bracket.upper;
	      
	      /*	  printf("%.12f %.18f %.12f %.18f %.12f %.18f\n", 
			  a, min_one_dim->f_lower, minimum, min_one_dim->f_minimum, b,min_one_dim->f_upper);*/
	      status = gsl_min_test_interval (bracket, EPSABS_LINE, EPSREL_LINE);
	    }
	  while (status == GSL_CONTINUE && iterations_line < MAX_ITERATIONS_LINE);
	  gsl_multimin_fdf_minimizer_step_with_value(s,
						     min_one_dim->minimum,
						     min_one_dim->f_minimum);
	}
      printf("%i: \n",iterations);
      printf("x "); gsl_vector_fprintf (stdout, s->history->x, "%g"); 
      printf("f(x) %g\n",min_one_dim->f_minimum);
      printf("\n");
      if (iterations%(2*f.n) == 0)
	{
	  gsl_multimin_fdf_minimizer_restart(s);
	  just_started = 1;
	}
      else
	{
	  just_started = 0;
	}
    }
  while (iterations <= 100 
	 && gsl_multimin_test_gradient_sqr_norm(s->history,EPSABS) == GSL_CONTINUE);
}
