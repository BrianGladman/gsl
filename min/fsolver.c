#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>

#include "min.h"

static int 
compute_f_values (gsl_function * f, double minimum, double * f_minimum,
                  gsl_interval x, double * f_lower, double * f_upper);


static int 
compute_f_values (gsl_function * f, double minimum, double * f_minimum,
                  gsl_interval x, double * f_lower, double * f_upper)
{
  SAFE_FUNC_CALL(f, x.lower, f_lower);
  SAFE_FUNC_CALL(f, x.upper, f_upper);
  SAFE_FUNC_CALL(f, minimum, f_minimum);

  return GSL_SUCCESS;
}

gsl_min_fminimizer *
gsl_min_fminimizer_alloc (const gsl_min_fminimizer_type * T, 
			 gsl_function * f, double minimum, gsl_interval x)
{
  int status ;

  gsl_min_fminimizer * s;

  double f_minimum, f_lower, f_upper;

  status = compute_f_values (f, minimum, &f_minimum, x, &f_lower, &f_upper);

  if (status != GSL_SUCCESS)
    {
      GSL_ERROR_RETURN ("bad function value", GSL_EBADFUNC, 0);
    }
  
  s = gsl_min_fminimizer_alloc_with_values (T, f, minimum, f_minimum, 
                                            x, f_lower, f_upper);

  return s;
}

int
gsl_min_fminimizer_set (gsl_min_fminimizer * s, 
                        gsl_function * f, double minimum, gsl_interval x)
{
  int status ;

  double f_minimum, f_lower, f_upper;

  status = compute_f_values (f, minimum, &f_minimum, x, &f_lower, &f_upper);

  if (status != GSL_SUCCESS)
    {
      return status ;
    }
  
  status = gsl_min_fminimizer_set_with_values (s, f, minimum, f_minimum, 
                                               x, f_lower, f_upper);
  return status;
}


gsl_min_fminimizer *
gsl_min_fminimizer_alloc_with_values (const gsl_min_fminimizer_type * T, 
                                      gsl_function * f, 
                                      double minimum, double f_minimum,
                                      gsl_interval x, 
                                      double f_lower, double f_upper)
{
  int status;

  gsl_min_fminimizer * s = 
    (gsl_min_fminimizer *) malloc (sizeof (gsl_min_fminimizer));

  if (s == 0)
    {
      GSL_ERROR_RETURN ("failed to allocate space for minimizer struct",
			GSL_ENOMEM, 0);
    };

  s->state = malloc (T->size);

  if (s->state == 0)
    {
      free (s);		/* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for minimizer state",
			GSL_ENOMEM, 0);
    };

  s->type = T ;

  status = gsl_min_fminimizer_set_with_values (s, f, minimum, f_minimum,
                                               x, f_lower, f_upper); 

  if (status != GSL_SUCCESS)
    {
      free (s->state);
      free (s);		/* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to set minimizer", status, 0);
    };

  return s;
}

int
gsl_min_fminimizer_set_with_values (gsl_min_fminimizer * s, gsl_function * f, 
                                    double minimum, double f_minimum, 
                                    gsl_interval x, 
                                    double f_lower, double f_upper)
{
  s->function = f;
  s->minimum = minimum;
  s->interval = x;

  if (x.lower > x.upper)
    {
      GSL_ERROR ("invalid interval (lower > upper)", GSL_EINVAL);
    }

  if (minimum >= x.upper || minimum <= x.lower) 
    {
      GSL_ERROR ("minimum must lie inside interval (lower < x < upper)",
                 GSL_EINVAL);
    }

  s->f_lower = f_lower;
  s->f_upper = f_upper;
  s->f_minimum = f_minimum;

  if (f_minimum >= f_lower || f_minimum >= f_upper)
    {
      GSL_ERROR ("endpoints do not enclose a minimum", GSL_EINVAL);
    }

  return (s->type->set) (s->state, s->function, 
                         minimum, f_minimum, 
                         x, f_lower, f_upper);
}


int
gsl_min_fminimizer_iterate (gsl_min_fminimizer * s)
{
  return (s->type->iterate) (s->state, s->function, 
                             &(s->minimum), &(s->f_minimum),
                             &(s->interval), &(s->f_lower), &(s->f_upper));
}

void
gsl_min_fminimizer_free (gsl_min_fminimizer * s)
{
  free (s->state);
  free (s);
}

const char *
gsl_min_fminimizer_name (const gsl_min_fminimizer * s)
{
  return s->type->name;
}

double
gsl_min_fminimizer_minimum (const gsl_min_fminimizer * s)
{
  return s->minimum;
}

gsl_interval
gsl_min_fminimizer_interval (const gsl_min_fminimizer * s)
{
  return s->interval;
}

