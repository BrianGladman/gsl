/* test -- program to test root finding functions
   Interfaces correctly with `make check'. */

/* config headers */
#include <config.h>

/* standard headers */
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>

/* gsl headers */
#include <gsl_errno.h>
#include <gsl_math.h>
#include <gsl_roots.h>

/* roots headers */
#include "roots.h"

/* test headers */
#include "test.h"
#include <gsl_test.h>

/* stopping parameters */
const double REL_EPSILON = (10 * DBL_EPSILON * GSL_ROOT_EPSILON_BUFFER);
const double ABS_EPSILON = (10 * DBL_EPSILON * GSL_ROOT_EPSILON_BUFFER);
const unsigned int MAX_ITERATIONS = 100;

void my_error_handler (const char *reason, const char *file,
		       int line, int err);

int
main (void)
{

  gsl_set_error_handler (&my_error_handler);


  test_macros ();

  test_bisection ("gsl_root_bisection, sin(x) [3, 4]",
		  sin, 3.0, 4.0, M_PI);
  test_bisection ("gsl_root_bisection, sin(x) [-4, -3]",
		  sin, -4.0, -3.0, -M_PI);
  test_bisection ("gsl_root_bisection, sin(x) [-1/3, 1]",
		  sin, -1.0 / 3.0, 1.0, 0.0);
  test_bisection ("gsl_root_bisection, cos(x) [0, 3]",
		  cos, 0.0, 3.0, M_PI / 2.0);
  test_bisection ("gsl_root_bisection, cos(x) [-3, 0]",
		  cos, -3.0, 0.0, -M_PI / 2.0);
  test_bisection ("gsl_root_bisection, x^20 - 1 [0.1, 2]",
		  func1, 0.1, 2.0, 1.0);
  test_bisection ("gsl_root_bisection, sqrt(|x|)*sgn(x)",
		  func2, -1.0 / 3.0, 1.0, 0.0);
  test_bisection ("gsl_root_bisection, x^2 - 1e-8 [0, 1]",
		  func3, 0.0, 1.0, sqrt (1e-8));
  test_bisection ("gsl_root_bisection, x exp(-x) [-1/3, 2]",
		  func4, -1.0 / 3.0, 2.0, 0.0);
  test_bisection ("gsl_root_bisection, (x - 1)^7 [0.1, 2]",
		  func6, 0.1, 2.0, 1.0);

  test_bisection_failure ("gsl_root_bisection, invalid range check [4, 0]",
			  sin, 4.0, 0.0, M_PI);
  test_bisection_failure ("gsl_root_bisection, invalid range check [1, 1]",
			  sin, 1.0, 1.0, M_PI);

  /* Test false position. */

  test_falsepos ("gsl_root_falsepos, sin(x) [3, 3.2]",
		 sin, 3.0, 3.2, M_PI);
  test_falsepos ("gsl_root_falsepos, sin(x) [-4, -3]",
		 sin, -4.0, -3.0, -M_PI);
  test_falsepos ("gsl_root_falsepos, sin(x) [-1/3, 1]",
		 sin, -1.0 / 3.0, 1.0, 0.0);
  test_falsepos ("gsl_root_falsepos, cos(x) [0, 3]",
		 cos, 0.0, 3.0, M_PI / 2.0);
  test_falsepos ("gsl_root_falsepos, cos(x) [-3, 0]",
		 cos, -3.0, 0.0, -M_PI / 2.0);
  test_falsepos ("gsl_root_falsepos, x^{20} - 1 [0.99, 1.01]",
		 func1, 0.99, 1.01, 1.0);
  test_falsepos ("gsl_root_falsepos, sqrt(|x|)*sgn(x) [-1/3, 1]",
		 func2, -1.0 / 3.0, 1.0, 0.0);
  test_falsepos ("gsl_root_falsepos, x^2 - 1e-8 [5e-5, 3e-4]",
		 func3, 5e-5, 3.0e-4, sqrt (1e-8));
  test_falsepos ("gsl_root_falsepos, x exp(-x) [-1/3, 2]",
		 func4, -1.0 / 3.0, 2.0, 0.0);
  test_falsepos ("gsl_root_falsepos, (x - 1)^7 [0.8, 1.2]",
		 func6, 0.8, 1.2, 1.0);


  /* Test secant method. */

  test_secant ("gsl_root_secant, sin(x) {3.3, 3.4}",
	       sin, 3.3, 3.4, M_PI);
  test_secant ("gsl_root_secant, sin(x) {-3.3, -3.4}",
	       sin, -3.3, -3.4, -M_PI);
  test_secant ("gsl_root_secant, sin(x) {0.4, 0.5}",
	       sin, 0.4, 0.5, 0.0);
  test_secant ("gsl_root_secant, cos(x) {0.5, 0.6}",
	       cos, 0.5, 0.6, M_PI / 2.0);
  test_secant ("gsl_root_secant, cos(x) {-2.5, -3.0}",
	       cos, -2.5, -3.0, -M_PI / 2.0);
  test_secant ("gsl_root_secant, x^20 - 1 {0.9, 0.91}",
	       func1, 0.9, 0.91, 1.0);
  test_secant ("gsl_root_secant, x^20 - 1 {1.1, 1.11}",
	       func1, 1.1, 1.11, 1.0);
  test_secant ("gsl_root_secant, sqrt(abs(x)) * sgn(x) {1, 1.01}",
	       func2, 1.0, 1.01, 0.0);
  test_secant ("gsl_root_secant, x^2 - 1e-8 {1, 1.01}",
	       func3, 1.0, 1.01, sqrt (1e-8));

  test_secant_failure ("gsl_root_secant, max iterations x -> +Inf, x exp(-x) {2, 2.01}",
		       func4, 2.0, 2.01, 0.0);
  test_secant_failure ("gsl_root_secant, max iterations x -> -Inf, 1/(1 + exp(-x)) {0, 0.01}",
		       func5, 0.0, 0.01, 0.0);
  test_secant_failure ("gsl_root_secant, max iterations towards root, (x - 1)^7 {0, 0.01}",
		       func6, 0.0, 0.01, 1.0);

  /* Test Newton's Method. */

  test_newton ("gsl_root_newton, sin(x) {3.4}",
	       sin_f, sin_df, sin_fdf, 3.4, M_PI);
  test_newton ("gsl_root_newton, sin(x) {-3.3}",
	       sin_f, sin_df, sin_fdf, -3.3, -M_PI);
  test_newton ("gsl_root_newton, sin(x) {0.5}",
	       sin_f, sin_df, sin_fdf, 0.5, 0.0);
  test_newton ("gsl_root_newton, cos(x) {0.6}",
	       cos_f, cos_df, cos_fdf, 0.6, M_PI / 2.0);
  test_newton ("gsl_root_newton, cos(x) {-2.5}",
	       cos_f, cos_df, cos_fdf, -2.5, -M_PI / 2.0);
  test_newton ("gsl_root_newton, x^{20} - 1 {0.9}",
	       func1, func1_df, func1_fdf, 0.9, 1.0);
  test_newton ("gsl_root_newton, x^{20} - 1 {1.1}",
	       func1, func1_df, func1_fdf, 1.1, 1.0);
  test_newton ("gsl_root_newton, sqrt(|x|)*sgn(x) {1.001}",
	       func2, func2_df, func2_fdf, 0.001, 0.0);
  test_newton ("gsl_root_newton, x^2 - 1e-8 {1}",
	       func3, func3_df, func3_fdf, 1.0, sqrt (1e-8));
  test_newton ("gsl_root_newton, x exp(-x) {-2}",
	       func4, func4_df, func4_fdf, -2.0, 0.0);

  test_newton_failure ("gsl_root_newton, max iterations x -> +Inf, x exp(-x) {2}",
		       func4, func4_df, func4_fdf, 2.0, 0.0);
  test_newton_failure ("gsl_root_newton, max iterations x -> -Inf, 1/(1 + exp(-x)) {0}",
		       func5, func5_df, func5_fdf, 0.0, 0.0);
  
  test_newton_failure ("gsl_root_newton, max iterations towards root, (x - 1)^7 {0}",
		       func6, func6_df, func6_fdf, 0.0, 1.0);


  /* now summarize the results */

  return gsl_test_summary ();
}


/* Test certain macros. */
void
test_macros (void)
{
  int result;
  double inf, nan ;

  /* 1.0 is real */
  result = GSL_ISREAL (1.0);
  gsl_test (result != 1, "GSL_ISREAL(1.0) is 1");

  inf = 1.0 / (sqrt(1.0) - 1) ;

  /* 1.0/0.0 == Inf is not real */
  result = GSL_ISREAL (inf);
  gsl_test (result != 0, "GSL_ISREAL(Inf) is 0");

  nan = inf - inf ;

  /* 0.0/0.0 == NaN is not real */
  result = GSL_ISREAL (nan);
  gsl_test (result != 0, "GSL_ISREAL(NaN) is 0");
}

/* Using gsl_root_bisection, find the root of the function pointed to by f,
   using the interval [lower_bound, upper_bound]. Check if f succeeded and
   that it was accurate enough. */

void
test_bisection (const char *description,
		double (*f) (double),
		double lower_bound, double upper_bound,
		double correct_root)
{
  int status;
  double root;

  status = gsl_root_bisection (&root, f, &lower_bound, &upper_bound,
			       REL_EPSILON, ABS_EPSILON,
			       MAX_ITERATIONS);

  gsl_test (status, description, root - correct_root);

  /* check the validity of the returned result */

  if (!WITHIN_TOL (root, correct_root, REL_EPSILON, ABS_EPSILON))
    {
      status = 1; /* failed */ ;
      gsl_test (status, "precision incorrectly reported");
    }
}

void
test_bisection_failure (const char *description,
			double (*f) (double),
			double lower_bound, double upper_bound,
			double correct_root)
{
  int status;
  double root;

  status = gsl_root_bisection (&root, f, &lower_bound, &upper_bound,
			       REL_EPSILON, ABS_EPSILON,
			       MAX_ITERATIONS);

  gsl_test (!status, description, root - correct_root);
}


/* Using gsl_root_falsepos, find the root of the function pointed to by f,
   using the interval [lower_bound, upper_bound]. Check if f succeeded and
   that it was accurate enough. */

void
test_falsepos (const char *description,
	       double (*f) (double),
	       double lower_bound, double upper_bound,
	       const double correct_root)
{
  int status;
  double root;

  status = gsl_root_falsepos (&root, f, &lower_bound, &upper_bound,
			      REL_EPSILON, ABS_EPSILON,
			      MAX_ITERATIONS);

  gsl_test (status, description);

  /* check the validity of the returned result */

  if (!WITHIN_TOL (root, correct_root, REL_EPSILON, ABS_EPSILON))
    {
      status = 1; /* failed */ ;
      gsl_test (status, "precision incorrectly reported");
    }
}

void
test_falsepos_failure (const char *description,
		       double (*f) (double),
		       double lower_bound, double upper_bound,
		       double correct_root)
{
  int status;
  double root;

  status = gsl_root_falsepos (&root, f, &lower_bound, &upper_bound,
			      REL_EPSILON, ABS_EPSILON,
			      MAX_ITERATIONS);

  gsl_test (!status, description, root - correct_root);
}


/* Using gsl_root_secant, find the root of the function pointed to by f, with
   guesses guess1 and guess2. Check if f succeeded and that it was accurate
   enough. */

void
test_secant (const char *description,
	     double (*f) (double),
	     double lower_bound, double upper_bound,
	     double correct_root)
{
  int status;
  double root;

  status = gsl_root_secant (&root, f, &lower_bound, &upper_bound,
			    REL_EPSILON, ABS_EPSILON,
			    MAX_ITERATIONS);

  /* check the validity of the returned result */

  gsl_test (status, description, root - correct_root);

  if (!WITHIN_TOL(root, correct_root, REL_EPSILON, ABS_EPSILON))
    {
      status = 1; /* failed */ ;
      gsl_test (status, "precision incorrectly reported");
    }
}

void
test_secant_failure (const char *description,
		     double (*f) (double),
		     double lower_bound, double upper_bound,
		     double correct_root)
{
  int status;
  double root;

  status = gsl_root_secant (&root, f, &lower_bound, &upper_bound,
			    REL_EPSILON, ABS_EPSILON,
			    MAX_ITERATIONS);

  gsl_test (!status, description, root - correct_root);
}


/* Using gsl_root_newton, find the root of the function pointed to by fdf,
   with guess guess. Check if f succeeded and that it was accurate enough. */
void
  test_newton (const char *description,
	       double (*f) (double),
	       double (*df) (double),
	       void (*fdf) (double, double *, double *),
	       double guess, double correct_root)
{
  int status;
  double root;

  status = gsl_root_newton (&root, f, df, fdf, &guess, 
			    REL_EPSILON, ABS_EPSILON,
			    MAX_ITERATIONS);

  gsl_test (status, description, root - correct_root);

  /* check the validity of the returned result */

  if (!WITHIN_TOL (root, correct_root, REL_EPSILON, ABS_EPSILON))
    {
      status = 1; /* failed */ ;
      gsl_test (status, "precision incorrectly reported (%g obs vs %g expected)", root, correct_root);
    }

}

void
  test_newton_failure (const char *description,
		       double (*f) (double),
		       double (*df) (double),
		       void (*fdf) (double, double *, double *),
		       double guess, double correct_root)
{
  int status;
  double root;

  status = gsl_root_newton (&root, f, df, fdf, &guess, 
			    REL_EPSILON, ABS_EPSILON,
			    MAX_ITERATIONS);

  gsl_test (!status, description, root - correct_root);
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
  return - sin (x);
}

void
cos_fdf (double x, double *y, double *yprime)
{
  *y = cos (x);
  *yprime = -sin (x);
}


void
my_error_handler (const char *reason, const char *file, int line, int err)
{
  if (0)
    printf ("(caught [%s:%d: %s (%d)])\n", file, line, reason, err);
}


