#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include <gsl_roots.h>

int
gsl_root_test_interval (const gsl_interval x, 
			double rel_epsilon, double abs_epsilon)
{
  const double lower = x.lower;
  const double upper = x.upper;

  const double abs_lower = fabs(lower) ;
  const double abs_upper = fabs(upper) ;

  const double max_abs = GSL_MAX_DBL(abs_lower, abs_upper) ;
  const double min_abs = GSL_MIN_DBL(abs_lower, abs_upper) ;

  const double tolerance = abs_epsilon + rel_epsilon * min_abs  ;

  if (rel_epsilon < 0.0 || abs_epsilon < 0.0)
    GSL_ERROR ("relative or absolute tolerance negative", GSL_EBADTOL);
  
  if (rel_epsilon < GSL_DBL_EPSILON * GSL_ROOT_EPSILON_BUFFER)
    GSL_ERROR ("relative tolerance too small", GSL_EBADTOL);

  if (lower > upper)
    GSL_ERROR ("lower bound larger than upper_bound", GSL_EINVAL);
  
  if (tolerance < max_abs * GSL_DBL_EPSILON * GSL_ROOT_EPSILON_BUFFER)
    GSL_ERROR("tolerances too small for this context", GSL_EBADTOL); 
  
  if (fabs(upper - lower) < tolerance)
    return GSL_SUCCESS;
  
  /* GSL_ERROR ("exceeded maximum number of iterations", GSL_EMAXITER); */

  return GSL_CONTINUE ;
}

int
gsl_root_test_delta (double x1, double x0, 
			double rel_epsilon, double abs_epsilon)
{
  const double abs_x1 = fabs(x1) ;
  const double abs_x0 = fabs(x0) ;

  const double max_abs = GSL_MAX_DBL(abs_x0, abs_x1) ;
  const double min_abs = GSL_MIN_DBL(abs_x0, abs_x1) ;

  const double tolerance = abs_epsilon + rel_epsilon * min_abs  ;

  if (rel_epsilon < 0.0 || abs_epsilon < 0.0)
    GSL_ERROR ("relative or absolute tolerance negative", GSL_EBADTOL);
  
  if (rel_epsilon < GSL_DBL_EPSILON * GSL_ROOT_EPSILON_BUFFER)
    GSL_ERROR ("relative tolerance too small", GSL_EBADTOL);

  if (tolerance < max_abs * GSL_DBL_EPSILON * GSL_ROOT_EPSILON_BUFFER)
    GSL_ERROR("tolerances too small for this context", GSL_EBADTOL); 
  
  if (fabs(x1 - x0) < tolerance)
    return GSL_SUCCESS;
  
  /* GSL_ERROR ("exceeded maximum number of iterations", GSL_EMAXITER); */

  return GSL_CONTINUE ;
}

int
gsl_root_test_residual (double f, double abs_epsilon)
{
  if (abs_epsilon < 0.0)
    GSL_ERROR ("absolute tolerance is negative", GSL_EBADTOL);
 
  if (fabs(f) < abs_epsilon)
    return GSL_SUCCESS;
  
  /* GSL_ERROR ("exceeded maximum number of iterations", GSL_EMAXITER); */

  return GSL_CONTINUE ;
}

