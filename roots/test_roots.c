#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_test.h>

#include "roots.h"
#include "test.h"

/* stopping parameters */
const double EPSREL = (10 * GSL_DBL_EPSILON);
const double EPSABS = (10 * GSL_DBL_EPSILON);
const unsigned int MAX_ITERATIONS = 150;

void my_error_handler (const char *reason, const char *file,
		       int line, int err);

#define WITHIN_TOL(a, b, epsrel, epsabs) \
 ((fabs((a) - (b)) < (epsrel) * GSL_MIN(fabs(a), fabs(b)) + (epsabs)))

void
test_roots (void)
{
  gsl_function F_sin, F_cos, F_func1, F_func2, F_func3, F_func4,
    F_func5, F_func6;
  
  gsl_function_fdf FDF_sin, FDF_cos, FDF_func1, FDF_func2, FDF_func3, FDF_func4,
    FDF_func5, FDF_func6;

  const gsl_root_fsolver_type * fsolver[4] ;
  const gsl_root_fdfsolver_type * fdfsolver[4] ;

  const gsl_root_fsolver_type ** T;
  const gsl_root_fdfsolver_type ** S;


  fsolver[0] = gsl_root_fsolver_bisection;
  fsolver[1] = gsl_root_fsolver_brent;
  fsolver[2] = gsl_root_fsolver_falsepos;
  fsolver[3] = 0;

  fdfsolver[0] = gsl_root_fdfsolver_newton;
  fdfsolver[1] = gsl_root_fdfsolver_secant;
  fdfsolver[2] = gsl_root_fdfsolver_steffenson;
  fdfsolver[3] = 0;

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

  for (T = fsolver ; *T != 0 ; T++)
    {
      test_f (*T, "sin(x) [3, 4]", &F_sin, 3.0, 4.0, M_PI);
      test_f (*T, "sin(x) [-4, -3]", &F_sin, -4.0, -3.0, -M_PI);
      test_f (*T, "sin(x) [-1/3, 1]", &F_sin, -1.0 / 3.0, 1.0, 0.0);
      test_f (*T, "cos(x) [0, 3]", &F_cos, 0.0, 3.0, M_PI / 2.0);
      test_f (*T, "cos(x) [-3, 0]", &F_cos, -3.0, 0.0, -M_PI / 2.0);
      test_f (*T, "x^20 - 1 [0.1, 2]", &F_func1, 0.1, 2.0, 1.0);
      test_f (*T, "sqrt(|x|)*sgn(x)", &F_func2, -1.0 / 3.0, 1.0, 0.0);
      test_f (*T, "x^2 - 1e-8 [0, 1]", &F_func3, 0.0, 1.0, sqrt (1e-8));
      test_f (*T, "x exp(-x) [-1/3, 2]", &F_func4, -1.0 / 3.0, 2.0, 0.0);
      test_f (*T, "(x - 1)^7 [0.9995, 1.0002]", &F_func6, 0.9995, 1.0002, 1.0);
      
      test_f_e (*T, "invalid range check [4, 0]", &F_sin, 4.0, 0.0, M_PI);
      test_f_e (*T, "invalid range check [1, 1]", &F_sin, 1.0, 1.0, M_PI);
      test_f_e (*T, "invalid range check [0.1, 0.2]", &F_sin, 0.1, 0.2, M_PI);
    }

  for (S = fdfsolver ; *S != 0 ; S++)
    {
      test_fdf (*S,"sin(x) {3.4}", &FDF_sin, 3.4, M_PI);
      test_fdf (*S,"sin(x) {-3.3}", &FDF_sin, -3.3, -M_PI);
      test_fdf (*S,"sin(x) {0.5}", &FDF_sin, 0.5, 0.0);
      test_fdf (*S,"cos(x) {0.6}", &FDF_cos, 0.6, M_PI / 2.0);
      test_fdf (*S,"cos(x) {-2.5}", &FDF_cos, -2.5, -M_PI / 2.0);
      test_fdf (*S,"x^{20} - 1 {0.9}", &FDF_func1, 0.9, 1.0);
      test_fdf (*S,"x^{20} - 1 {1.1}", &FDF_func1, 1.1, 1.0);
      test_fdf (*S,"sqrt(|x|)*sgn(x) {1.001}", &FDF_func2, 0.001, 0.0);
      test_fdf (*S,"x^2 - 1e-8 {1}", &FDF_func3, 1.0, sqrt (1e-8));
      test_fdf (*S,"x exp(-x) {-2}", &FDF_func4, -2.0, 0.0);
      test_fdf_e (*S,"max iterations x -> +Inf, x exp(-x) {2}", &FDF_func4, 2.0, 0.0);
      test_fdf_e (*S,"max iterations x -> -Inf, 1/(1 + exp(-x)) {0}", &FDF_func5, 0.0, 0.0);
    }

  test_fdf (gsl_root_fdfsolver_steffenson,
	    "(x - 1)^7 {0.9}", &FDF_func6, 0.9, 1.0);    

  


}

/* Using gsl_root_bisection, find the root of the function pointed to by f,
   using the interval [lower_bound, upper_bound]. Check if f succeeded and
   that it was accurate enough. */

void
test_f (const gsl_root_fsolver_type * T, const char * description, gsl_function *f,
	double lower_bound, double upper_bound, double correct_root)
{
  int status;
  size_t iterations = 0;
  double r, a, b;
  gsl_interval x;
  gsl_root_fsolver * s;

  x.lower = lower_bound;
  x.upper = upper_bound;

  s = gsl_root_fsolver_alloc(T, f, x) ;
  
  do 
    {
      iterations++ ;

      gsl_root_fsolver_iterate (s);

      r = gsl_root_fsolver_root(s);
      x = gsl_root_fsolver_interval(s);
      
      a = x.lower;
      b = x.upper;

      if (a > b)
	gsl_test (GSL_FAILURE, "interval is invalid (%g,%g)", a, b);

      if (r < a || r > b)
	gsl_test (GSL_FAILURE, "r lies outside interval %g (%g,%g)", r, a, b);

      status = gsl_root_test_interval (x, EPSABS, EPSREL);
    }
  while (status == GSL_CONTINUE && iterations < MAX_ITERATIONS);


  gsl_test (status, "%s, %s (%g obs vs %g expected) ", 
	    gsl_root_fsolver_name(s), description, 
	    gsl_root_fsolver_root(s), correct_root);

  if (iterations == MAX_ITERATIONS)
    {
      gsl_test (GSL_FAILURE, "exceeded maximum number of iterations");
    }

  /* check the validity of the returned result */

  if (!WITHIN_TOL (r, correct_root, EPSREL, EPSABS))
    {
      gsl_test (GSL_FAILURE, "incorrect precision (%g obs vs %g expected)", 
		r, correct_root);

    }
}

void
test_f_e (const gsl_root_fsolver_type * T, 
	  const char * description, gsl_function *f,
	  double lower_bound, double upper_bound, double correct_root)
{
  int status;
  size_t iterations = 0;
  gsl_interval x;
  gsl_root_fsolver * s;

  x.lower = lower_bound;
  x.upper = upper_bound;

  s = gsl_root_fsolver_alloc(T, f, x) ;

  if (s == 0) 
    {
      gsl_test (s != 0, "%s, %s", T->name, description);
      return ;
    }

  do 
    {
      iterations++ ;
      gsl_root_fsolver_iterate (s);
      status = gsl_root_test_interval (gsl_root_fsolver_interval(s), 
				      EPSABS, EPSREL);
    }
  while (status == GSL_CONTINUE && iterations < MAX_ITERATIONS);

  gsl_test (!status, "%s, %s", gsl_root_fsolver_name(s), description, 
	    gsl_root_fsolver_root(s) - correct_root);

}

void
test_fdf (const gsl_root_fdfsolver_type * T, const char * description, 
	gsl_function_fdf *fdf, double root, double correct_root)
{
  int status;
  size_t iterations = 0;
  double prev = 0 ;

  gsl_root_fdfsolver * s = gsl_root_fdfsolver_alloc(T, fdf, root) ;

  do 
    {
      iterations++ ;
      prev = gsl_root_fdfsolver_root(s);
      gsl_root_fdfsolver_iterate (s);
      status = gsl_root_test_delta(gsl_root_fdfsolver_root(s), prev, 
				   EPSABS, EPSREL);
    }
  while (status == GSL_CONTINUE && iterations < MAX_ITERATIONS);

  gsl_test (status, "%s, %s (%g obs vs %g expected) ", 
	    gsl_root_fdfsolver_name(s), description, 
	    gsl_root_fdfsolver_root(s), correct_root);

  if (iterations == MAX_ITERATIONS)
    {
      gsl_test (GSL_FAILURE, "exceeded maximum number of iterations");
    }

  /* check the validity of the returned result */

  if (!WITHIN_TOL (gsl_root_fdfsolver_root(s), correct_root, 
		   EPSREL, EPSABS))
    {
      gsl_test (GSL_FAILURE, "incorrect precision (%g obs vs %g expected)", 
		gsl_root_fdfsolver_root(s), correct_root);

    }
}

void
test_fdf_e (const gsl_root_fdfsolver_type * T, 
	    const char * description, gsl_function_fdf *fdf,
	    double root, double correct_root)
{
  int status;
  size_t iterations = 0;
  double prev = 0 ;

  gsl_root_fdfsolver * s = gsl_root_fdfsolver_alloc(T, fdf, root) ;

  if (s == 0) 
    {
      gsl_test (s != 0, "%s, %s", T->name, description);
      return ;
    }

  do 
    {
      iterations++ ;
      prev = gsl_root_fdfsolver_root(s);
      gsl_root_fdfsolver_iterate (s);
      status = gsl_root_test_delta(gsl_root_fdfsolver_root(s), prev, 
				   EPSABS, EPSREL);
    }
  while (status == GSL_CONTINUE && iterations < MAX_ITERATIONS);

  gsl_test (!status, "%s, %s", gsl_root_fdfsolver_name(s), 
	    description, gsl_root_fdfsolver_root(s) - correct_root);
}

void
my_error_handler (const char *reason, const char *file, int line, int err)
{
  if (0)
    printf ("(caught [%s:%d: %s (%d)])\n", file, line, reason, err);
}


