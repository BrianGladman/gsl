#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>

int
gsl_min_test_interval (const gsl_interval x, double epsabs, double epsrel)
{
  const double lower = x.lower;
  const double upper = x.upper;

  const double abs_lower = fabs(lower) ;
  const double abs_upper = fabs(upper) ;

  double min_abs, tolerance;

  if (epsrel < 0.0)
    GSL_ERROR ("relative tolerance is negative", GSL_EBADTOL);
  
  if (epsabs < 0.0)
    GSL_ERROR ("absolute tolerance is negative", GSL_EBADTOL);

  if (lower > upper)
    GSL_ERROR ("lower bound larger than upper_bound", GSL_EINVAL);

  if ((lower > 0 && upper > 0) || (lower < 0 && upper < 0)) 
    {
      min_abs = GSL_MIN_DBL(abs_lower, abs_upper) ;
    }
  else
    {
      min_abs = 0;
    }

  tolerance = epsabs + epsrel * min_abs  ;
  
  if (fabs(upper - lower) < tolerance)
    return GSL_SUCCESS;
  
  return GSL_CONTINUE ;
}

