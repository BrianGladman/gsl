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
#include <gsl_roots.h>

/* roots headers */
#include "roots.h"

/* test headers */
#include "test.h"
#include <gsl_test.h>

/* stopping parameters */
#define REL_EPSILON    (10 * DBL_EPSILON * GSL_ROOT_EPSILON_BUFFER)
#define ABS_EPSILON    (10 * DBL_EPSILON * GSL_ROOT_EPSILON_BUFFER)
#define MAX_ITERATIONS 100
#define MAX_DELTAY     2000000.0
#define MAX_STEP_SIZE  100.0

void my_error_handler (const char *reason, const char *file,
		       int line, int err);

/* Print usage instructions. */
void
usage(void)
{
  printf("Usage:  test <srcdir> <tests>\n"
"where <srcdir> is the path to the directory with the source code\n"
"where <tests> is a string indicating which tests to run. It can contain the\n"
"following characters:\n"
"\n"
"  m -- test macros\n"
"  b -- test gsl_root_bisection\n"
"  f -- test gsl_root_falsepos\n"
"  s -- test gsl_root_secant\n"
"  n -- test gsl_root_newton\n"
"\n"
"Example:\n"
"\n"
"  test mb\n"
"\n"
"tests macros and gsl_root_bisection.\n");
} ;

int
main(int argc, char ** argv)
{
  int status ;
  int num_passed = 0, num_failed = 0;
  char *srcdir = NULL;

  /* Do some basic checking of the arguments. */
  {
    char c;
   
    /* The user must specifiy one argument. */
    if (argc != 3) {
      usage();
      return 1;
    }

    srcdir = (char *) malloc(strlen(argv[1])*sizeof(char) + 1);
    strcpy(srcdir, argv[1]);
    printf("running tests in the source directory %s\n", srcdir);

    /* argv[2] can contain only the characters "mbfsn". */
    for (c = ' '; c <= '~'; c++)
      if (strchr(argv[2], c) && !strchr("mbfsn", c)) {
        usage();
        return 1;
      }
  }

  gsl_set_error_handler(&my_error_handler);

  

  /* Test macros if so instructed. */
  if (strchr(argv[2], 'm')) {
    test_macros() ;
  }

  /* Test bisection if so instructed. */

  if (strchr(argv[2], 'b')) 
    {
      test_bisection("gsl_root_bisection, sin(x) [3, 4]", 
		     sin, 3.0, 4.0, M_PI); 
      test_bisection("gsl_root_bisection, sin(x) [-4, -3]", 
		     sin, -4.0, -3.0, -M_PI); 
      test_bisection("gsl_root_bisection, sin(x) [-1/3, 1]",
		     sin, -1.0/3.0, 1.0, 0.0);
      test_bisection("gsl_root_bisection, cos(x) [0, 3]",
		     cos, 0.0, 3.0, M_PI/2.0);
      test_bisection("gsl_root_bisection, cos(x) [-3, 0]",
		     cos, -3.0, 0.0, -M_PI/2.0);
      test_bisection("gsl_root_bisection, x^20 - 1 [0.1, 2]",
		     test_hairy_1, 0.1, 2.0, 1.0); 
      test_bisection("gsl_root_bisection, sqrt(|x|)*sgn(x)",
		     test_hairy_2, -1.0/3.0, 1.0, 0.0);
      test_bisection("gsl_root_bisection, x^2 - 1e-8 [0, 1]",
		     test_hairy_3, 0.0, 1.0, sqrt(1e-8));
      test_bisection("gsl_root_bisection, x exp(-x) [-1/3, 2]",
		     test_hairy_4, -1.0/3.0, 2.0, 0.0);
      test_bisection("gsl_root_bisection, (x - 1)^7 [0.1, 2]",
		     test_hairy_6, 0.1, 2.0, 1.0);
      
      test_bisection_failure("gsl_root_bisection, invalid range check [4, 0]",
			     sin, 4.0, 0.0, M_PI);
      test_bisection_failure("gsl_root_bisection, invalid range check [1, 1]",
		     sin, 1.0, 1.0, M_PI);
    }
  
  /* Test false position if so instructed. */
  if (strchr(argv[2], 'f')) 
    {
      test_falsepos("gsl_root_falsepos, sin(x) [3, 3.2]",
		    sin, 3.0, 3.2, M_PI); 
      test_falsepos("gsl_root_falsepos, sin(x) [-4, -3]", 
		    sin, -4.0, -3.0, -M_PI); 
      test_falsepos("gsl_root_falsepos, sin(x) [-1/3, 1]",
		    sin, -1.0/3.0, 1.0, 0.0); 
      test_falsepos("gsl_root_falsepos, cos(x) [0, 3]",
		    cos, 0.0, 3.0, M_PI/2.0); 
      test_falsepos("gsl_root_falsepos, cos(x) [-3, 0]", 
		    cos, -3.0, 0.0, -M_PI/2.0); 
      test_falsepos("gsl_root_falsepos, x^{20} - 1 [0.99, 1.01]", 
		    test_hairy_1, 0.99, 1.01, 1.0); 
      test_falsepos("gsl_root_falsepos, sqrt(|x|)*sgn(x) [-1/3, 1]", 
		    test_hairy_2, -1.0/3.0, 1.0, 0.0);
      test_falsepos("gsl_root_falsepos, x^2 - 1e-8 [0, 1]", 
		    test_hairy_3, 0.0, 1.0, sqrt(1e-8));
      test_falsepos("gsl_root_falsepos, x exp(-x) [-1/3, 2]", 
		    test_hairy_4, -1.0/3.0, 2.0, 0.0) ;
      test_falsepos("gsl_root_falsepos, (x - 1)^7 [pi/10, 2]", 
		    test_hairy_6, M_PI/10.0, 2.0, 1.0);

  }

  /* Test secant method if so instructed. */
  if (strchr(argv[2], 's')) 
    {
      test_secant("gsl_root_secant, sin(x) {3.3, 3.4}", 
		  sin, 3.3, 3.4, M_PI); 
      test_secant("gsl_root_secant, sin(x) {-3.3, -3.4}",
		  sin, -3.3, -3.4, -M_PI); 
      test_secant("gsl_root_secant, sin(x) {0.4, 0.5}",
		  sin, 0.4, 0.5, 0.0); 
      test_secant("gsl_root_secant, cos(x) {0.5, 0.6}",
		  cos, 0.5, 0.6, M_PI/2.0); 
      test_secant("gsl_root_secant, cos(x) {-2.5, -3.0}",
		  cos, -2.5, -3.0, -M_PI/2.0);
      test_secant("gsl_root_secant, x^20 - 1 {0.9, 0.91}",
		  test_hairy_1, 0.9, 0.91, 1.0); 
      test_secant("gsl_root_secant, x^20 - 1 {1.1, 1.11}",
		  test_hairy_1, 1.1, 1.11, 1.0); 
      test_secant("gsl_root_secant, sqrt(abs(x)) * sgn(x) {1, 1.01}",
		  test_hairy_2, 1.0, 1.01, 0.0);
      test_secant("gsl_root_secant, x^2 - 1e-8 {1, 1.01}",
		  test_hairy_3, 1.0, 1.01, sqrt(1e-8)) ;

      test_secant_failure("gsl_root_secant, max iterations x -> +Inf, x exp(-x) {2, 2.01}",
		  test_hairy_4, 2.0, 2.01, 0.0);
      test_secant_failure("gsl_root_secant, max iterations x -> -Inf, 1/(1 + exp(-x)) {0, 0.01}",
			  test_hairy_5, 0.0, 0.01, 0.0);
      test_secant_failure("gsl_root_secant, max iterations towards root, (x - 1)^7 {0, 0.01}",
		  test_hairy_6, 0.0, 0.01, 1.0);

  }

  /* Test Newton's Method if so instructed. */
  if (strchr(argv[2], 'n')) {
    test_newton("gsl_root_newton, sin(x) {3.4}",
		sin_fdf, 3.4, M_PI); 
    test_newton("gsl_root_newton, sin(x) {-3.3}",
		sin_fdf, -3.3, -M_PI);
    test_newton("gsl_root_newton, sin(x) {0.5}",
		sin_fdf, 0.5, 0.0); 
    test_newton("gsl_root_newton, cos(x) {0.6}",
		cos_fdf, 0.6, M_PI/2.0);
    test_newton("gsl_root_newton, cos(x) {-2.5}",
		cos_fdf, -2.5, -M_PI/2.0);
    test_newton("gsl_root_newton, x^{20} - 1 {0.9}",
		test_hairy_1_fdf, 0.9, 1.0); 
    test_newton("gsl_root_newton, x^{20} - 1 {1.1}",
		test_hairy_1_fdf, 1.1, 1.0); 
    test_newton("gsl_root_newton, sqrt(|x|)*sgn(x) {1.001}",
		test_hairy_2_fdf, 0.001, 0.0);
    test_newton("gsl_root_newton, x^2 - 1e-8 {1}",
		test_hairy_3_fdf, 1.0, sqrt(1e-8));
    test_newton("gsl_root_newton, x exp(-x) {-2}",
		test_hairy_4_fdf, -2.0, 0.0);

    test_newton_failure("gsl_root_newton, max iterations x -> +Inf, x exp(-x) {2}",
			test_hairy_4_fdf, 2.0, 0.0);
    test_newton_failure("gsl_root_newton, max iterations x -> -Inf, 1/(1 + exp(-x)) {0}",
		test_hairy_5_fdf, 0.0, 0.0);

    test_newton_failure("gsl_root_newton, max iterations towards root, (x - 1)^7 {0}",
		test_hairy_6_fdf, 0.0, 1.0);

  }

  /* now summarize the results */

  return gsl_test_summary() ;
}


/* Test certain macros. */
void
test_macros () 
{
  int result ;

  /* 1.0 is real */
  result = GSL_ISREAL(1.0) ;
  gsl_test(result != 1,"GSL_ISREAL(1.0) is 1");

  /* 1.0/0.0 == Inf is not real */
  result = GSL_ISREAL(1.0/0.0) ;
  gsl_test(result != 0, "GSL_ISREAL(Inf) is 0");

  /* 0.0/0.0 == NaN is not real */
  result = GSL_ISREAL(0.0/0.0) ;
  gsl_test(result != 0,"GSL_ISREAL(NaN) is 0") ;
}

/* Using gsl_root_bisection, find the root of the function pointed to by f,
   using the interval [lower_bound, upper_bound]. Check if f succeeded and
   that it was accurate enough. */

int
test_bisection (const char * description, 
		double (* f)(double), 
		double lower_bound, double upper_bound, 
		double correct_root)
{
  int status ;
  double root;

  status = gsl_root_bisection(&root, f, &lower_bound, &upper_bound, 
			      REL_EPSILON, ABS_EPSILON, 
			      MAX_ITERATIONS, MAX_DELTAY);

  gsl_test (status, description, root - correct_root);

  /* check the validity of the returned result */
  
  if (!_WITHIN_TOL(root, correct_root, REL_EPSILON, ABS_EPSILON)) 
    {
      status = 1 ; /* failed */ ;
      gsl_test (status, "precision incorrectly reported");
    }
  
  return status ;

} 

int
test_bisection_failure (const char * description, 
			double (* f)(double), 
			double lower_bound, double upper_bound, 
			double correct_root)
{
  int status ;
  double root;
  
  status = gsl_root_bisection(&root, f, &lower_bound, &upper_bound, 
			      REL_EPSILON, ABS_EPSILON, 
			      MAX_ITERATIONS, MAX_DELTAY);
  
  gsl_test (!status, description, root - correct_root);
  
}


/* Using gsl_root_falsepos, find the root of the function pointed to by f,
   using the interval [lower_bound, upper_bound]. Check if f succeeded and
   that it was accurate enough. */

int
test_falsepos (const char * description, 
		double (* f)(double), 
		double lower_bound, double upper_bound, 
		double correct_root)
{
  int status ;
  double root;

  status = gsl_root_falsepos (&root, f, &lower_bound, &upper_bound, 
			      REL_EPSILON, ABS_EPSILON, 
			      MAX_ITERATIONS, MAX_DELTAY);

  gsl_test (status, description, root - correct_root);
  
  /* check the validity of the returned result */
  
  if (!_WITHIN_TOL(root, correct_root, REL_EPSILON, ABS_EPSILON)) 
    {
      status = 1 ; /* failed */ ;
      gsl_test (status, "precision incorrectly reported");
    }

  return status ;

} 

int
test_falsepos_failure (const char * description, 
		       double (* f)(double), 
		       double lower_bound, double upper_bound, 
		       double correct_root)
{
  int status ;
  double root;

  status = gsl_root_falsepos (&root, f, &lower_bound, &upper_bound, 
			      REL_EPSILON, ABS_EPSILON, 
			      MAX_ITERATIONS, MAX_DELTAY);
  
  gsl_test (!status, description, root - correct_root);
}


/* Using gsl_root_secant, find the root of the function pointed to by f, with
   guesses guess1 and guess2. Check if f succeeded and that it was accurate
   enough. */

int
test_secant (const char * description, 
	     double (* f)(double), 
	     double lower_bound, double upper_bound, 
	     double correct_root)
{
  int status ;
  double root;

  status = gsl_root_secant(&root, f, &lower_bound, &upper_bound, 
			      REL_EPSILON, ABS_EPSILON, 
			      MAX_ITERATIONS, MAX_STEP_SIZE);
  
  /* check the validity of the returned result */
  
  gsl_test (status, description, root - correct_root);

  if (!_WITHIN_TOL(root, correct_root, REL_EPSILON, ABS_EPSILON)) 
    {
      status = 1 ; /* failed */ ;
      gsl_test (status, "precision incorrectly reported");
    }
  
  return status ;

} 

int
test_secant_failure (const char * description, 
		double (* f)(double), 
		double lower_bound, double upper_bound, 
		double correct_root)
{
  int status ;
  double root;
  
  status = gsl_root_secant(&root, f, &lower_bound, &upper_bound, 
			      REL_EPSILON, ABS_EPSILON, 
			      MAX_ITERATIONS, MAX_STEP_SIZE);
  
  gsl_test (!status, description, root - correct_root);
}


/* Using gsl_root_newton, find the root of the function pointed to by fdf,
   with guess guess. Check if f succeeded and that it was accurate enough. */
int
test_newton (const char * description,
	     void (* fdf)(double *, double *, double, int, int),
	     double guess, double correct_root)
{
  int status;
  double root;
  
  status = gsl_root_newton(&root, fdf, &guess, REL_EPSILON, ABS_EPSILON,
			   MAX_ITERATIONS, MAX_STEP_SIZE);

  gsl_test (status, description, root - correct_root);

  /* check the validity of the returned result */
  
  if (!_WITHIN_TOL(root, correct_root, REL_EPSILON, ABS_EPSILON)) 
    {
      status = 1 ; /* failed */ ;
      gsl_test (status, "precision incorrectly reported (%g obs vs %g expected)",root, correct_root);
    }

  return status ;
} 

int
test_newton_failure (const char * description, 
		     void (* fdf)(double *, double *, double, int, int),
		     double guess, double correct_root)
{
  int status ;
  double root;

  status = gsl_root_newton(&root, fdf, &guess, REL_EPSILON, ABS_EPSILON,
			   MAX_ITERATIONS, MAX_STEP_SIZE);
  
  gsl_test (!status, description, root - correct_root);
}



/* f(x) = x^{20} - 1 */
/* f'(x) = 20x^{19} */
/* zero at x = 1 or -1 */
double
test_hairy_1(double x)
{
  return pow(x, 20.0) - 1;
}

void
test_hairy_1_fdf(double * y, double * yprime, double x, int y_wanted,
                          int yprime_wanted)
{
  *y = test_hairy_1(x);
  *yprime = 20.0 * pow(x, 19.0);
}

/* f(x) = sqrt(abs(x))*sgn(x) */
/* f'(x) = 1 / sqrt(abs(x) */
/* zero at x = 0 */
double
test_hairy_2(double x)
{
  double delta;

  if (x > 0)
    delta = 1.0;
  else if (x < 0)
    delta = -1.0;
  else
    delta = 0.0;

  return sqrt(fabs(x))*delta;
}

void
test_hairy_2_fdf(double * y, double * yprime, double x, int y_wanted,
                          int yprime_wanted)
{
  *y = test_hairy_2(x);
  *yprime = 1 / sqrt(fabs(x));
}


/* f(x) = x^2 - 1e-8 */
/* f'(x) = 2x */
/* zero at x = sqrt(1e-8) or -sqrt(1e-8) */
double
test_hairy_3(double x)
{
  return pow(x, 2.0) - 1e-8;
}

void
test_hairy_3_fdf(double * y, double * yprime, double x, int y_wanted,
                          int yprime_wanted)
{
  *y = test_hairy_3(x);
  *yprime = 2 * x;
}

/* f(x) = x exp(-x) */
/* f'(x) = exp(-x) - x exp(-x) */
/* zero at x = 0 */
double
test_hairy_4(double x)
{
  return x * exp(-x);
}

void
test_hairy_4_fdf(double * y, double * yprime, double x, int y_wanted,
                          int yprime_wanted)
{
  *y = test_hairy_4(x);
  *yprime = exp(-x) - x * exp(-x);
}

/* f(x) = 1/(1+exp(x)) */
/* f'(x) = -exp(x) / (1 + exp(x))^2 */
/* no roots! */
double
test_hairy_5(double x)
{
  return 1 / (1 + exp(x));
}

void
test_hairy_5_fdf(double * y, double * yprime, double x, int y_wanted,
                          int yprime_wanted)
{
  *y = test_hairy_5(x);
  *yprime = -exp(x) / pow(1 + exp(x), 2.0);
}

/* f(x) = (x - 1)^7 */
/* f'(x) = 7 * (x - 1)^6 */
/* zero at x = 1 */
double
test_hairy_6(double x)
{
  return pow(x - 1, 7.0);
}

void
test_hairy_6_fdf(double * y, double * yprime, double x, int y_wanted,
                          int yprime_wanted)
{
  *y = test_hairy_6(x);
  *yprime = 7.0 * pow(x - 1, 6.0);
}

/* sin(x) packaged up nicely. */
void
sin_fdf(double * y, double * yprime, double x, int y_wanted,
        int yprime_wanted)
{
  *y = sin(x);
  *yprime = cos(x);
}

/* cos(x) packaged up nicely. */
void
cos_fdf(double * y, double * yprime, double x, int y_wanted,
        int yprime_wanted)
{
  *y = cos(x);
  *yprime = -sin(x);
}


void
my_error_handler (const char *reason, const char *file, int line, int err)
{
  if(0) printf ("(caught [%s:%d: %s (%d)])\n", file, line, reason, err) ;
}
