#include <math.h>
#include <stdlib.h>
#include <gsl_math.h>

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

gsl_fdf create_fdf (simple_function * f, simple_function * df) 
{
  gsl_fdf FDF ;
  function_pair * fp = malloc(sizeof(function_pair));
  fp->f = f ; 
  fp->df = df ;
  FDF.f = eval_fdf_f ;
  FDF.df = eval_fdf_df ;
  FDF.fdf = eval_fdf ;
  FDF.params = fp ;
  return FDF ;
}

double eval_fdf_f (double x, void * params) 
{
  function_pair *fp = (function_pair *)params ;
  return (*(fp->f))(x) ;
}

double eval_fdf_df (double x, void * params) 
{
  function_pair *fp = (function_pair *)params ;
  return (*(fp->df))(x) ;
}

void eval_fdf (double x, void * params, double * y1, double * y2) 
{
  function_pair *fp = (function_pair *)params ;
  *y1 = (*(fp->f))(x) ;
  *y2 = (*(fp->df))(x) ;
}

/* f(x) = x^{20} - 1 */
/* f'(x) = 20x^{19} */
/* zero at x = 1 or -1 */

double
func1 (double x)
{
  return pow (x, 20.0) - 1;
}

double
func1_df (double x)
{
  return 20.0 * pow (x, 19.0);
}

void
func1_fdf (double x, double *y, double *yprime)
{
  *y = func1 (x);
  *yprime = 20.0 * pow (x, 19.0);
}

/* f(x) = sqrt(abs(x))*sgn(x) */
/* f'(x) = 1 / sqrt(abs(x) */
/* zero at x = 0 */
double
func2 (double x)
{
  double delta;

  if (x > 0)
    delta = 1.0;
  else if (x < 0)
    delta = -1.0;
  else
    delta = 0.0;

  return sqrt (fabs (x)) * delta;
}

double
func2_df (double x)
{
  return 1 / sqrt (fabs (x));
}

void
func2_fdf (double x, double *y, double *yprime)
{
  *y = func2 (x);
  *yprime = 1 / sqrt (fabs (x));
}


/* f(x) = x^2 - 1e-8 */
/* f'(x) = 2x */
/* zero at x = sqrt(1e-8) or -sqrt(1e-8) */
double
func3 (double x)
{
  return pow (x, 2.0) - 1e-8;
}

double
func3_df (double x)
{
  return 2 * x;
}

void
func3_fdf (double x, double *y, double *yprime)
{
  *y = func3 (x);
  *yprime = 2 * x;
}

/* f(x) = x exp(-x) */
/* f'(x) = exp(-x) - x exp(-x) */
/* zero at x = 0 */
double
func4 (double x)
{
  return x * exp (-x);
}

double
func4_df (double x)
{
  return exp (-x) - x * exp (-x);
}

void
func4_fdf (double x, double *y, double *yprime)
{
  *y = func4 (x);
  *yprime = exp (-x) - x * exp (-x);
}

/* f(x) = 1/(1+exp(x)) */
/* f'(x) = -exp(x) / (1 + exp(x))^2 */
/* no roots! */
double
func5 (double x)
{
  return 1 / (1 + exp (x));
}

double
func5_df (double x)
{
  return -exp (x) / pow (1 + exp (x), 2.0);
}

void
func5_fdf (double x, double *y, double *yprime)
{
  *y = func5 (x);
  *yprime = -exp (x) / pow (1 + exp (x), 2.0);
}

/* f(x) = (x - 1)^7 */
/* f'(x) = 7 * (x - 1)^6 */
/* zero at x = 1 */
double
func6 (double x)
{
  return pow (x - 1, 7.0);
}

double
func6_df (double x)
{
  return 7.0 * pow (x - 1, 6.0);
}

void
func6_fdf (double x, double *y, double *yprime)
{
  *y = func6 (x);
  *yprime = 7.0 * pow (x - 1, 6.0);
}

/* sin(x) packaged up nicely. */
double
sin_f (double x)
{
  return sin (x);
}

double
sin_df (double x)
{
  return cos (x);
}

void
sin_fdf (double x, double *y, double *yprime)
{
  *y = sin (x);
  *yprime = cos (x);
}

/* cos(x) packaged up nicely. */
double
cos_f (double x)
{
  return cos (x);
}

double
cos_df (double x)
{
  return -sin (x);
}

void
cos_fdf (double x, double *y, double *yprime)
{
  *y = cos (x);
  *yprime = -sin (x);
}
