#include <gsl_math.h>
#include <gsl_roots.h>
#include <gsl_errno.h>
#include <gsl_test.h>

#include "roots.h"
#include "test.h"

/* stopping parameters */
const double REL_EPSILON = (10 * GSL_DBL_EPSILON * GSL_ROOT_EPSILON_BUFFER);
const double ABS_EPSILON = (10 * GSL_DBL_EPSILON * GSL_ROOT_EPSILON_BUFFER);
const unsigned int MAX_ITERATIONS = 100;


void my_error_handler (const char *reason, const char *file,
		       int line, int err);

void
test_roots (void)
{
  gsl_function F_sin, F_cos, F_func1, F_func2, F_func3, F_func4,
    F_func5, F_func6;
  
  gsl_fdf FDF_sin, FDF_cos, FDF_func1, FDF_func2, FDF_func3, FDF_func4,
    FDF_func5, FDF_func6;

  F_sin = create_function (sin) ;
  F_cos = create_function (cos) ; 
  F_func1 = create_function (func1) ;
  F_func2 = create_function (func2) ;
  F_func3 = create_function (func3) ;
  F_func4 = create_function (func4) ;
  F_func5 = create_function (func5) ;
  F_func6 = create_function (func6) ;

  FDF_sin = create_fdf (sin_f, sin_df) ;
  FDF_cos = create_fdf (cos_f, cos_df) ;
  FDF_func1 = create_fdf (func1, func1_df) ;
  FDF_func2 = create_fdf (func2, func2_df) ;
  FDF_func3 = create_fdf (func3, func3_df) ;
  FDF_func4 = create_fdf (func4, func4_df) ;
  FDF_func5 = create_fdf (func5, func5_df) ;
  FDF_func6 = create_fdf (func6, func6_df) ;

  gsl_set_error_handler (&my_error_handler);

  test_macros ();

  test_bisection ("gsl_root_bisection, sin(x) [3, 4]",
		  &F_sin, 3.0, 4.0, M_PI);
  test_bisection ("gsl_root_bisection, sin(x) [-4, -3]",
		  &F_sin, -4.0, -3.0, -M_PI);
  test_bisection ("gsl_root_bisection, sin(x) [-1/3, 1]",
		  &F_sin, -1.0 / 3.0, 1.0, 0.0);
  test_bisection ("gsl_root_bisection, cos(x) [0, 3]",
		  &F_cos, 0.0, 3.0, M_PI / 2.0);
  test_bisection ("gsl_root_bisection, cos(x) [-3, 0]",
		  &F_cos, -3.0, 0.0, -M_PI / 2.0);
  test_bisection ("gsl_root_bisection, x^20 - 1 [0.1, 2]",
		  &F_func1, 0.1, 2.0, 1.0);
  test_bisection ("gsl_root_bisection, sqrt(|x|)*sgn(x)",
		  &F_func2, -1.0 / 3.0, 1.0, 0.0);
  test_bisection ("gsl_root_bisection, x^2 - 1e-8 [0, 1]",
		  &F_func3, 0.0, 1.0, sqrt (1e-8));
  test_bisection ("gsl_root_bisection, x exp(-x) [-1/3, 2]",
		  &F_func4, -1.0 / 3.0, 2.0, 0.0);
  test_bisection ("gsl_root_bisection, (x - 1)^7 [0.1, 2]",
		  &F_func6, 0.1, 2.0, 1.0);

  test_bisection_failure ("gsl_root_bisection, invalid range check [4, 0]",
			  &F_sin, 4.0, 0.0, M_PI);
  test_bisection_failure ("gsl_root_bisection, invalid range check [1, 1]",
			  &F_sin, 1.0, 1.0, M_PI);
  test_bisection_failure ("gsl_root_bisection, invalid range check [0.1, 0.2]",
			  &F_sin, 0.1, 0.2, M_PI);

  test_brent ("gsl_root_brent, sin(x) [3, 4]",
		  &F_sin, 3.0, 4.0, M_PI);
  test_brent ("gsl_root_brent, sin(x) [-4, -3]",
		  &F_sin, -4.0, -3.0, -M_PI);
  test_brent ("gsl_root_brent, sin(x) [-1/3, 1]",
		  &F_sin, -1.0 / 3.0, 1.0, 0.0);
  test_brent ("gsl_root_brent, cos(x) [0, 3]",
		  &F_cos, 0.0, 3.0, M_PI / 2.0);
  test_brent ("gsl_root_brent, cos(x) [-3, 0]",
		  &F_cos, -3.0, 0.0, -M_PI / 2.0);
  test_brent ("gsl_root_brent, x^20 - 1 [0.1, 2]",
		  &F_func1, 0.1, 2.0, 1.0);
  test_brent ("gsl_root_brent, sqrt(|x|)*sgn(x)",
		  &F_func2, -1.0 / 3.0, 1.0, 0.0);
  test_brent ("gsl_root_brent, x^2 - 1e-8 [0, 1]",
		  &F_func3, 0.0, 1.0, sqrt (1e-8));
  test_brent ("gsl_root_brent, x exp(-x) [-1/3, 2]",
		  &F_func4, -1.0 / 3.0, 2.0, 0.0);
  test_brent ("gsl_root_brent, (x - 1)^7 [0.1, 2]",
		  &F_func6, 0.1, 2.0, 1.0);

  test_brent_failure ("gsl_root_brent, invalid range check [4, 0]",
			  &F_sin, 4.0, 0.0, M_PI);
  test_brent_failure ("gsl_root_brent, invalid range check [1, 1]",
			  &F_sin, 1.0, 1.0, M_PI);
  test_brent_failure ("gsl_root_brent, invalid range check [0.1, 0.2]",
			  &F_sin, 0.1, 0.2, M_PI);

  /* Test false position. */

  test_falsepos ("gsl_root_falsepos, sin(x) [3, 3.2]",
		 &F_sin, 3.0, 3.2, M_PI);
  test_falsepos ("gsl_root_falsepos, sin(x) [-4, -3]",
		 &F_sin, -4.0, -3.0, -M_PI);
  test_falsepos ("gsl_root_falsepos, sin(x) [-1/3, 1]",
		 &F_sin, -1.0 / 3.0, 1.0, 0.0);
  test_falsepos ("gsl_root_falsepos, cos(x) [0, 3]",
		 &F_cos, 0.0, 3.0, M_PI / 2.0);
  test_falsepos ("gsl_root_falsepos, cos(x) [-3, 0]",
		 &F_cos, -3.0, 0.0, -M_PI / 2.0);
  test_falsepos ("gsl_root_falsepos, x^{20} - 1 [0.99, 1.01]",
		 &F_func1, 0.99, 1.01, 1.0);
  test_falsepos ("gsl_root_falsepos, sqrt(|x|)*sgn(x) [-1/3, 1]",
		 &F_func2, -1.0 / 3.0, 1.0, 0.0);
  test_falsepos ("gsl_root_falsepos, x^2 - 1e-8 [5e-5, 3e-4]",
		 &F_func3, 5e-5, 3.0e-4, sqrt (1e-8));
  test_falsepos ("gsl_root_falsepos, x exp(-x) [-1/3, 2]",
		 &F_func4, -1.0 / 3.0, 2.0, 0.0);
  test_falsepos ("gsl_root_falsepos, (x - 1)^7 [0.8, 1.2]",
		 &F_func6, 0.8, 1.2, 1.0);


  /* Test secant method. */

  test_secant ("gsl_root_secant, sin(x) {3.3, 3.4}",
	       &F_sin, 3.3, 3.4, M_PI);
  test_secant ("gsl_root_secant, sin(x) {-3.3, -3.4}",
	       &F_sin, -3.3, -3.4, -M_PI);
  test_secant ("gsl_root_secant, sin(x) {0.4, 0.5}",
	       &F_sin, 0.4, 0.5, 0.0);
  test_secant ("gsl_root_secant, cos(x) {0.5, 0.6}",
	       &F_cos, 0.5, 0.6, M_PI / 2.0);
  test_secant ("gsl_root_secant, cos(x) {-2.5, -3.0}",
	       &F_cos, -2.5, -3.0, -M_PI / 2.0);
  test_secant ("gsl_root_secant, x^20 - 1 {0.9, 0.91}",
	       &F_func1, 0.9, 0.91, 1.0);
  test_secant ("gsl_root_secant, x^20 - 1 {1.1, 1.11}",
	       &F_func1, 1.1, 1.11, 1.0);
  test_secant ("gsl_root_secant, sqrt(abs(x)) * sgn(x) {1, 1.01}",
	       &F_func2, 1.0, 1.01, 0.0);
  test_secant ("gsl_root_secant, x^2 - 1e-8 {1, 1.01}",
	       &F_func3, 1.0, 1.01, sqrt (1e-8));

  test_secant_failure ("gsl_root_secant, max iterations x -> +Inf, x exp(-x) {2, 2.01}",
		       &F_func4, 2.0, 2.01, 0.0);
  test_secant_failure ("gsl_root_secant, max iterations x -> -Inf, 1/(1 + exp(-x)) {0, 0.01}",
		       &F_func5, 0.0, 0.01, 0.0);
  test_secant_failure ("gsl_root_secant, max iterations towards root, (x - 1)^7 {0, 0.01}",
		       &F_func6, 0.0, 0.01, 1.0);

  /* Test Newton's Method. */

  test_newton ("gsl_root_newton, sin(x) {3.4}",
	       &FDF_sin, 3.4, M_PI);
  test_newton ("gsl_root_newton, sin(x) {-3.3}",
	       &FDF_sin, -3.3, -M_PI);
  test_newton ("gsl_root_newton, sin(x) {0.5}",
	       &FDF_sin, 0.5, 0.0);
  test_newton ("gsl_root_newton, cos(x) {0.6}",
	       &FDF_cos, 0.6, M_PI / 2.0);
  test_newton ("gsl_root_newton, cos(x) {-2.5}",
	       &FDF_cos, -2.5, -M_PI / 2.0);
  test_newton ("gsl_root_newton, x^{20} - 1 {0.9}",
	       &FDF_func1, 0.9, 1.0);
  test_newton ("gsl_root_newton, x^{20} - 1 {1.1}",
	       &FDF_func1, 1.1, 1.0);
  test_newton ("gsl_root_newton, sqrt(|x|)*sgn(x) {1.001}",
	       &FDF_func2, 0.001, 0.0);
  test_newton ("gsl_root_newton, x^2 - 1e-8 {1}",
	       &FDF_func3, 1.0, sqrt (1e-8));
  test_newton ("gsl_root_newton, x exp(-x) {-2}",
	       &FDF_func4, -2.0, 0.0);

  test_newton_failure ("gsl_root_newton, max iterations x -> +Inf, x exp(-x) {2}",
		       &FDF_func4, 2.0, 0.0);
  test_newton_failure ("gsl_root_newton, max iterations x -> -Inf, 1/(1 + exp(-x)) {0}",
		       &FDF_func5, 0.0, 0.0);
  
  test_newton_failure ("gsl_root_newton, max iterations towards root, (x - 1)^7 {0}",
		       &FDF_func6, 0.0, 1.0);
}


/* Using gsl_root_bisection, find the root of the function pointed to by f,
   using the interval [lower_bound, upper_bound]. Check if f succeeded and
   that it was accurate enough. */

void
test_bisection (const char *description,
		const gsl_function *f,
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
			const gsl_function *f,
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


void
test_brent (const char *description,
		const gsl_function *f,
		double lower_bound, double upper_bound,
		double correct_root)
{
  int status;
  double root;

  /* The Brent algorithm can be up to 3 times slower than bisection on
     pathological cases, so we'll allow more iterations here. The
     factor of 3 is a guaranteed bound. */

  status = gsl_root_brent (&root, f, &lower_bound, &upper_bound,
			   REL_EPSILON, ABS_EPSILON,
			   3 * MAX_ITERATIONS);

  gsl_test (status, description, root - correct_root);

  /* check the validity of the returned result */

  if (!WITHIN_TOL (root, correct_root, REL_EPSILON, ABS_EPSILON))
    {
      status = 1; /* failed */ ;
      gsl_test (status, "precision incorrectly reported (%.18g vs %.18g)", 
		root, correct_root);
    }
}

void
test_brent_failure (const char *description,
		    const gsl_function *f,
		    double lower_bound, double upper_bound,
		    double correct_root)
{
  int status;
  double root;

  status = gsl_root_brent (&root, f, &lower_bound, &upper_bound,
			   REL_EPSILON, ABS_EPSILON,
			   MAX_ITERATIONS);
  
  gsl_test (!status, description, root - correct_root);
}


/* Using gsl_root_falsepos, find the root of the function pointed to by f,
   using the interval [lower_bound, upper_bound]. Check if f succeeded and
   that it was accurate enough. */

void
test_falsepos (const char *description,
	       const gsl_function *f,
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
		       const gsl_function *f,
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
	     const gsl_function *f,
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
		     const gsl_function *f,
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
	       const gsl_fdf *fdf,
	       double guess, double correct_root)
{
  int status;
  double root;

  status = gsl_root_newton (&root, fdf, &guess, 
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
		     const gsl_fdf *fdf,
		     double guess, double correct_root)
{
  int status;
  double root;

  status = gsl_root_newton (&root, fdf, &guess, 
			    REL_EPSILON, ABS_EPSILON,
			    MAX_ITERATIONS);

  gsl_test (!status, description, root - correct_root);
}

void
my_error_handler (const char *reason, const char *file, int line, int err)
{
  if (0)
    printf ("(caught [%s:%d: %s (%d)])\n", file, line, reason, err);
}


