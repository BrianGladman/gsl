#include <gsl_math.h>
#include <gsl_min.h>
#include <gsl_errno.h>
#include <gsl_test.h>

#include "test.h"

/* stopping parameters */

const double EPSABS = 0.0001 ;
const double EPSREL = 0.0001 ;

const unsigned int MAX_ITERATIONS = 100;

void my_error_handler (const char *reason, const char *file,
		       int line, int err);

#define WITHIN_TOL(a, b, epsrel, epsabs) \
 (fabs((a) - (b)) < (epsrel) * GSL_MIN(fabs(a), fabs(b)) + (epsabs))

int
main (void)
{
  gsl_function F_cos;
  
  const gsl_min_fsolver_type * fsolver[4] ;
  const gsl_min_fsolver_type ** T;

  fsolver[0] = gsl_min_fsolver_bisection;
  fsolver[1] = 0;

  F_cos = create_function (cos) ;

  gsl_set_error_handler (&my_error_handler);

  for (T = fsolver ; *T != 0 ; T++)
    {
      test_f (*T, "cos(x) [0 (3) 6]", &F_cos, 0.0, 3.0, 6.0, M_PI);

#ifdef JUNK
      test_f_e (*T, "invalid range check [4, 0]", &F_sin, 4.0, 0.0, M_PI);
      test_f_e (*T, "invalid range check [1, 1]", &F_sin, 1.0, 1.0, M_PI);
      test_f_e (*T, "invalid range check [0.1, 0.2]", &F_sin, 0.1, 0.2, M_PI);
#endif
    }

  return gsl_test_summary ();
}

void
test_f (const gsl_min_fsolver_type * T, 
        const char * description, gsl_function *f,
	double lower_bound, double middle, double upper_bound, 
        double correct_minimum)
{
  int status;
  size_t iterations = 0;
  double m, a, b;
  gsl_interval x;
  gsl_min_fsolver * s;

  x.lower = lower_bound;
  x.upper = upper_bound;

  s = gsl_min_fsolver_alloc(T, f, middle, x) ;
  
  do 
    {
      iterations++ ;

      gsl_min_fsolver_iterate (s);

      m = gsl_min_fsolver_minimum(s);
      x = gsl_min_fsolver_interval(s);
      
      a = x.lower;
      b = x.upper;

      if (a > b)
	gsl_test (GSL_FAILURE, "interval is invalid (%g,%g)", a, b);

      if (m < a || m > b)
	gsl_test (GSL_FAILURE, "r lies outside interval %g (%g,%g)", m, a, b);

      status = gsl_min_test_interval (x, EPSABS, EPSREL);
    }
  while (status == GSL_CONTINUE && iterations < MAX_ITERATIONS);

  gsl_test (status, "%s, %s (%g obs vs %g expected) ", 
	    gsl_min_fsolver_name(s), description, 
	    gsl_min_fsolver_minimum(s), correct_minimum);

  /* check the validity of the returned result */

  if (!WITHIN_TOL (m, correct_minimum, EPSREL, EPSABS))
    {
      gsl_test (GSL_FAILURE, "incorrect precision (%g obs vs %g expected)", 
		m, correct_minimum);
    }
}

#ifdef JUNK
void
test_f_e (const gsl_min_fsolver_type * T, 
	  const char * description, gsl_function *f,
	  double lower_bound, double upper_bound, double correct_root)
{
  int status;
  size_t iterations = 0;
  gsl_interval x;
  gsl_min_fsolver * s;

  x.lower = lower_bound;
  x.upper = upper_bound;

  s = gsl_min_fsolver_alloc(T, f, x) ;

  if (s == 0) 
    {
      gsl_test (s != 0, "%s, %s", T->name, description);
      return ;
    }

  do 
    {
      iterations++ ;
      gsl_min_fsolver_iterate (s);
      status = gsl_min_test_interval (gsl_min_fsolver_interval(s), 
				      EPSABS, EPSREL);
    }
  while (status == GSL_CONTINUE && iterations < MAX_ITERATIONS);

  gsl_test (!status, "%s, %s", gsl_min_fsolver_name(s), description, 
	    gsl_min_fsolver_root(s) - correct_root);

}
#endif

void
my_error_handler (const char *reason, const char *file, int line, int err)
{
  if (0)
    printf ("(caught [%s:%d: %s (%d)])\n", file, line, reason, err);
}
