#include <gsl_math.h>
#include <gsl_min.h>

#include "test.h"

gsl_function create_function (simple_function * f) 
{
  gsl_function F ;
  F.function = eval_function ;
  F.params = (void *)f ;
  return F ;
}

double eval_function (double x, void * params) 
{
  simple_function *f = (simple_function *)params ;
  return (*f)(x) ;
}


/* f(x) = x^4 - 1 */
/* minimum at x = 0 */

double
func1 (double x)
{
  return pow (x, 4.0) - 1;
}

/* f(x) = sqrt(|x|) */
/* minimum at x = 0 */

double
func2 (double x)
{
  return sqrt(fabs(x));
}


