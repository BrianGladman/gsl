/* test -- program to test root finding functions
   Interfaces correctly with `make check'. */
/* $Id$ */


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


/* stopping parameters */
#define TEST_REL_EPSILON    10 * DBL_EPSILON * GSL_ROOT_EPSILON_BUFFER
/* #define TEST_REL_EPSILON    0.0000001 */
#define TEST_ABS_EPSILON    10 * DBL_EPSILON * GSL_ROOT_EPSILON_BUFFER
/* #define TEST_ABS_EPSILON    0.0000001 */
#define TEST_MAX_ITERATIONS 100
#define TEST_MAX_DELTAY     2000000.0
#define TEST_MAX_STEP_SIZE  100.0

#define FAIL 1
#define SUCC 0

/* Some globals to store data about errors that occur. Any call to
   error_handler will overwrite these! We set the strings to NULL so we can
   pass them to free() safely. */
char * g_reason = NULL; /* the reason for an error */
char * g_file = NULL;   /* the file where an error occured */
int g_line;             /* the line where an error occured */


/* administer all the tests */
int
main(int argc, char ** argv)
{
  int num_passed = 0, num_failed = 0;

  /* Do some basic checking of the arguments. */
  {
    char c;
   
    /* The user must specifiy one argument. */
    if (argc != 2) {
      print_usage();
      return 1;
    }

    /* It can contain only the characters "mbfsn". */
    for (c = ' '; c <= '~'; c++)
      if (strchr(argv[1], c) && !strchr("mbfsn", c)) {
        print_usage();
        return 1;
      }
  }

  gsl_set_error_handler(error_handler);


  /* Test macros if so instructed. */
  if (strchr(argv[1], 'm')) {
    printf("Testing root finding macros:\n");
    gsl_test_root_macros(&num_passed, &num_failed);
  }

  /* Test bisection if so instructed. */
  if (strchr(argv[1], 'b')) {
    printf("Testing `gsl_root_bisection': \n");
    printf("  sin(x) [3, 4] ... ");
    gsl_test_root_bisection(sin, 3.0, 4.0, M_PI, &num_passed, &num_failed); 
    printf("  sin(x) [-4, -3] ... ");
    gsl_test_root_bisection(sin, -4.0, -3.0, -M_PI, &num_passed, &num_failed); 
    printf("  sin(x) [-1/3, 1] ... ");
    gsl_test_root_bisection(sin, -1.0/3.0, 1.0, 0.0, &num_passed, &num_failed);
    printf("  cos(x) [0, 3] ... ");
    gsl_test_root_bisection(cos, 0.0, 3.0, M_PI/2.0, &num_passed, &num_failed);
    printf("  cos(x) [-3, 0] ... ");
    gsl_test_root_bisection(cos, -3.0, 0.0, -M_PI/2.0, &num_passed,
                            &num_failed); 
    printf("  x^20 - 1 [0.1, 2] ... ");
    gsl_test_root_bisection(gsl_test_root_hairy_1, 0.1, 2.0, 1.0, &num_passed,
                            &num_failed); 
    printf("  sqrt(abs(x)) * sgn(x) [-1/3, 1] ... ");
    gsl_test_root_bisection(gsl_test_root_hairy_2, -1.0/3.0, 1.0, 0.0,
                            &num_passed, &num_failed);
    printf("  x^2 - 1e-8 [0, 1] ... ");
    gsl_test_root_bisection(gsl_test_root_hairy_3, 0.0, 1.0, sqrt(1e-8),
                            &num_passed, &num_failed);
    printf("  x exp(-x) [-1/3, 2] ... ");
    gsl_test_root_bisection(gsl_test_root_hairy_4, -1.0/3.0, 2.0, 0.0,
                            &num_passed, &num_failed);
    printf("  (x - 1)^7 [pi/10, 2] ... ");
    gsl_test_root_bisection(gsl_test_root_hairy_6, M_PI/10.0, 2.0, 1.0,
                            &num_passed, &num_failed);
  }

  /* Test false position if so instructed. */
  if (strchr(argv[1], 'f')) {
    printf("Testing `gsl_root_falsepos': \n");
    printf("  sin(x) [3, 4] ... ");
    gsl_test_root_falsepos(sin, 3.0, 4.0, M_PI, &num_passed, &num_failed); 
    printf("  sin(x) [-4, -3] ... ");
    gsl_test_root_falsepos(sin, -4.0, -3.0, -M_PI, &num_passed, &num_failed); 
    printf("  sin(x) [-1/3, 1] ... ");
    gsl_test_root_falsepos(sin, -1.0/3.0, 1.0, 0.0, &num_passed, &num_failed); 
    printf("  cos(x) [0, 3] ... ");
    gsl_test_root_falsepos(cos, 0.0, 3.0, M_PI/2.0, &num_passed, &num_failed); 
    printf("  cos(x) [-3, 0] ... ");
    gsl_test_root_falsepos(cos, -3.0, 0.0, -M_PI/2.0, &num_passed,
                           &num_failed); 
    printf("  x^{20} - 1 [0.1, 2] ... ");
    gsl_test_root_falsepos(gsl_test_root_hairy_1, 0.1, 2.0, 1.0, &num_passed,
                           &num_failed); 
    printf("  sqrt(abs(x)) * sgn(x) [-1/3, 1] ... ");
    gsl_test_root_falsepos(gsl_test_root_hairy_2, -1.0/3.0, 1.0, 0.0,
                           &num_passed, &num_failed);
    printf("  x^2 - 1e-8 [0, 1] ... ");
    gsl_test_root_falsepos(gsl_test_root_hairy_3, 0.0, 1.0, sqrt(1e-8),
                           &num_passed, &num_failed);
    printf("  x exp(-x) [-1/3, 2] ... ");
    gsl_test_root_falsepos(gsl_test_root_hairy_4, -1.0/3.0, 2.0, 0.0,
                           &num_passed, &num_failed);
    printf("  (x - 1)^7 [pi/10, 2] ... ");
    gsl_test_root_falsepos(gsl_test_root_hairy_6, M_PI/10.0, 2.0, 1.0,
                           &num_passed, &num_failed);
  }

  /* Test secant method if so instructed. */
  if (strchr(argv[1], 's')) {
    printf("Testing `gsl_root_secant':\n");
    printf("  sin(x) {3.3, 3.4} ... ");
    gsl_test_root_secant(sin, 3.3, 3.4, M_PI, &num_passed, &num_failed, SUCC); 
    printf("  sin(x) {-3.3, -3.4} ... ");
    gsl_test_root_secant(sin, -3.3, -3.4, -M_PI, &num_passed, &num_failed,
                         SUCC); 
    printf("  sin(x) {0.4, 0.5} ... ");
    gsl_test_root_secant(sin, 0.4, 0.5, 0.0, &num_passed, &num_failed, SUCC); 
    printf("  cos(x) {0.5, 0.6} ... ");
    gsl_test_root_secant(cos, 0.5, 0.6, M_PI/2.0, &num_passed, &num_failed,
                         SUCC); 
    printf("  cos(x) {-2.5, -3.0} ... ");
    gsl_test_root_secant(cos, -2.5, -3.0, -M_PI/2.0, &num_passed, &num_failed,
                         SUCC);
    printf("  x^20 - 1 {0.9, 0.91} ... ");
    gsl_test_root_secant(gsl_test_root_hairy_1, 0.9, 0.91, 1.0, &num_passed,
                         &num_failed, SUCC); 
    printf("  x^20 - 1 {1.1, 1.11} ... ");
    gsl_test_root_secant(gsl_test_root_hairy_1, 1.1, 1.11, 1.0, &num_passed,
                         &num_failed, SUCC); 
    printf("  sqrt(abs(x)) * sgn(x) {1, 1.01} ... ");
    gsl_test_root_secant(gsl_test_root_hairy_2, 1.0, 1.01, 0.0,
                         &num_passed, &num_failed, SUCC);
    printf("  x^2 - 1e-8 {1, 1.01} ... ");
    gsl_test_root_secant(gsl_test_root_hairy_3, 1.0, 1.01, sqrt(1e-8),
                         &num_passed, &num_failed, SUCC);
    printf("  x exp(-x) {2, 2.01} ... ");
    gsl_test_root_secant(gsl_test_root_hairy_4, 2.0, 2.01, 0.0, &num_passed,
                         &num_failed, SUCC);
    printf("  1 / (1 + exp(-x)) {0, 0.01} ... ");
    gsl_test_root_secant(gsl_test_root_hairy_5, 0.0, 0.01, 0.0, &num_passed,
                         &num_failed, FAIL);
    printf("  (x - 1)^7 {0, 0.01} ... ");
    gsl_test_root_secant(gsl_test_root_hairy_6, 0.0, 0.01, 1.0, &num_passed,
                         &num_failed, SUCC);
  }

  /* Test Newton's Method if so instructed. */
  if (strchr(argv[1], 'n')) {
    printf("Testing `gsl_root_newton':\n");
    printf("  sin(x) {3.4} ... ");
    gsl_test_root_newton(sin_fdf, 3.4, M_PI, &num_passed, &num_failed, SUCC); 
    printf("  sin(x) {-3.3} ... ");
    gsl_test_root_newton(sin_fdf, -3.3, -M_PI, &num_passed, &num_failed, SUCC);
    printf("  sin(x) {0.5} ... ");
    gsl_test_root_newton(sin_fdf, 0.5, 0.0, &num_passed, &num_failed, SUCC); 
    printf("  cos(x) {0.6} ... ");
    gsl_test_root_newton(cos_fdf, 0.6, M_PI/2.0, &num_passed, &num_failed,
                         SUCC); 
    printf("  cos(x) {-2.5} ... ");
    gsl_test_root_newton(cos_fdf, -2.5, -M_PI/2.0, &num_passed, &num_failed,
                         SUCC);
    printf("  x^{20} - 1 {0.9} ... ");
    gsl_test_root_newton(gsl_test_root_hairy_1_fdf, 0.9, 1.0, &num_passed,
                         &num_failed, SUCC); 
    printf("  x^{20} - 1 {1.1} ... ");
    gsl_test_root_newton(gsl_test_root_hairy_1_fdf, 1.1, 1.0, &num_passed,
                         &num_failed, SUCC); 
    printf("  sqrt(abs(x)) * sgn(x) {1.001} ... ");
    gsl_test_root_newton(gsl_test_root_hairy_2_fdf, 1.001, 0.0, &num_passed,
                         &num_failed, SUCC);
    printf("  x^2 - 1e-8 {1} ... ");
    gsl_test_root_newton(gsl_test_root_hairy_3_fdf, 1.0, sqrt(1e-8),
                         &num_passed, &num_failed, SUCC);
    printf("  x exp(-x) {2} ... ");
    gsl_test_root_newton(gsl_test_root_hairy_4_fdf, 2.0, 0.0, &num_passed,
                         &num_failed, FAIL);
    printf("  x exp(-x) {-2} ... ");
    gsl_test_root_newton(gsl_test_root_hairy_4_fdf, -2.0, 0.0, &num_passed,
                         &num_failed, SUCC);
    printf("  1 / (1 + exp(-x)) {0} ... ");
    gsl_test_root_newton(gsl_test_root_hairy_5_fdf, 0.0, 0.0, &num_passed,
                         &num_failed, FAIL);
    printf("  (x - 1)^7 {0} ... ");
    gsl_test_root_newton(gsl_test_root_hairy_6_fdf, 0.0, 1.0, &num_passed,
                         &num_failed, SUCC);
  }

  /* now summarize the results */
  if (num_failed == 0) {
    printf("All %d tests passed.\n", num_passed);
    /* `make check' expects 0 if the test succeeded. */
    return 0;
  }
  else {
    printf("%d of %d tests failed.\n", num_failed, num_passed + num_failed);
    /* `make check' expects non-zero if the test failed. */
    return 1;
  }
}

/* Print usage instructions. */
void
print_usage(void)
{
  printf(
"Usage:

  test <tests>

where <tests> is a string indicating which tests to run. It can contain the
following characters:

  m -- test macros
  b -- test gsl_root_bisection
  f -- test gsl_root_falsepos
  s -- test gsl_root_secant
  n -- test gsl_root_newton

Example:

  test mb

tests macros and gsl_root_bisection.
");
}

/* An error handler which sets some global variables and returns. */
void
error_handler(const char * reason, const char * file, int line)
{
  /* Memory for g_reason and g_file was allocated in the last call to
     error_handler (or the are NULL because error_handler hasn't been called
     yet); we don't want a leak. */
  free(g_reason);
  free(g_file);
  /* strdup allocates memory internally. */
  g_reason = strdup(reason);
  g_file = strdup(file);
  if (g_reason == NULL || g_file == NULL) {
    printf("malloc failed\n");
    abort();
  }
  g_line = line;
}

/* Test certain macros. */
void
gsl_test_root_macros(int * num_passed, int * num_failed)
{
  /* GSL_ISREAL */
  /* 1.0 is real */
  printf("  GSL_ISREAL(1.0) ... ");
  if (GSL_ISREAL(1.0)) {
    printf("ok\n");
    (void)(*num_passed)++;
  }
  else {
    printf("failed\n");
    (void)(*num_failed)++;
  }
  /* 1.0/0.0 == Inf is not real */
  printf("  GSL_ISREAL(1.0/0.0) ... ");
  if (!GSL_ISREAL(1.0/0.0)) {
    printf("ok\n");
    (void)(*num_passed)++;
  }
  else {
    printf("failed\n");
    (void)(*num_failed)++;
  }
  /* 0.0/0.0 == NaN is not real */
  printf("  GSL_ISREAL(0.0/0.0) ... ");
  if (!GSL_ISREAL(0.0/0.0)) {
    printf("ok\n");
    (void)(*num_passed)++;
  }
  else {
    printf("failed\n");
    (void)(*num_failed)++;
  }
}

/* Using gsl_root_bisection, find the root of the function pointed to by f,
   using the interval [lower_bound, upper_bound]. Check if f succeeded and
   that it was accurate enough. */
void
gsl_test_root_bisection(double (* f)(double), double lower_bound,
                        double upper_bound, double cor_root, int * num_passed,
                        int * num_failed)
{
  int err;
  double root;
  
  err = gsl_root_bisection(&root, f, &lower_bound, &upper_bound,
                           TEST_REL_EPSILON, TEST_ABS_EPSILON,
                           TEST_MAX_ITERATIONS, TEST_MAX_DELTAY);
  /* if there was an error */
  if (err != GSL_SUCCESS) {
    printf("failed (error %d, %s)\n", gsl_errno, g_reason);
    (void)(*num_failed)++;
  }
  /* if it wasn't accurate enough */
  else if (!_WITHIN_TOL(root, cor_root, TEST_REL_EPSILON, TEST_ABS_EPSILON)) {
    printf("failed (inaccurate)\n");
    (void)(*num_failed)++;
  }
  /* if the test passed */
  else {
    printf("ok\n");
    (void)(*num_passed)++;
  }
}

/* Using gsl_root_falsepos, find the root of the function pointed to by f,
   using the interval [lower_bound, upper_bound]. Check if f succeeded and
   that it was accurate enough. */
void
gsl_test_root_falsepos(double (* f)(double), double lower_bound,
                        double upper_bound, double cor_root, int * num_passed,
                        int * num_failed)
{
  int err;
  double root;
  
  err = gsl_root_falsepos(&root, f, &lower_bound, &upper_bound,
                          TEST_REL_EPSILON, TEST_ABS_EPSILON,
                          TEST_MAX_ITERATIONS, TEST_MAX_DELTAY);
  /* if there was an error */
  if (err != GSL_SUCCESS) {
    printf("failed (error %d, %s)\n", gsl_errno, g_reason);
    (void)(*num_failed)++;
  }
  /* if it wasn't accurate enough */
  else if (!_WITHIN_TOL(root, cor_root, TEST_REL_EPSILON, TEST_ABS_EPSILON)) {
    printf("failed (inaccurate)\n");
    (void)(*num_failed)++;
  }
  /* if the test passed */
  else {
    printf("ok\n");
    (void)(*num_passed)++;
  }
}

/* Using gsl_root_secant, find the root of the function pointed to by f, with
   guesses guess1 and guess2. Check if f succeeded and that it was accurate
   enough. */
void
gsl_test_root_secant(double (* f)(double), double guess1, double guess2,
                     double cor_root, int * num_passed, int * num_failed,
                     int failure_desired)
{
  int err;
  double root;
  
  err = gsl_root_secant(&root, f, &guess1, &guess2, TEST_REL_EPSILON,
                        TEST_ABS_EPSILON, TEST_MAX_ITERATIONS,
                        TEST_MAX_STEP_SIZE);
  /* Were we supposed to fail? */
  if (failure_desired) {
    /* Did it not fail? */
    if (err == GSL_SUCCESS) {
      printf("failed (succeeded)\n");
      (void)(*num_failed)++;
    }
    /* It failed. */
    else {
      printf("ok (error %d)\n", gsl_errno);
      (void)(*num_passed)++;
    }
  }
  /* We were supposed to succeed. */
  else {
    /* Was there was an error? */
    if (err != GSL_SUCCESS) {
      printf("failed (error %d, %s)\n", gsl_errno, g_reason);
      (void)(*num_failed)++;
    }
    /* Was it not accurate enough? */
    else if (!_WITHIN_TOL(root, cor_root, TEST_REL_EPSILON,
                          TEST_ABS_EPSILON)) {
      printf("failed (inaccurate)\n");
      (void)(*num_failed)++;
    }
    /* The test passed. */
    else {
      printf("ok\n");
      (void)(*num_passed)++;
    }
  }
}

/* Using gsl_root_newton, find the root of the function pointed to by fdf,
   with guess guess. Check if f succeeded and that it was accurate enough. */
void
gsl_test_root_newton(void (* fdf)(double *, double *, double, int, int),
                     double guess, double cor_root, int * num_passed,
                     int * num_failed, int failure_desired)
{
  int err;
  double root;
  
  err = gsl_root_newton(&root, fdf, &guess, TEST_REL_EPSILON, TEST_ABS_EPSILON,
                        TEST_MAX_ITERATIONS, TEST_MAX_STEP_SIZE);
  /* Were we supposed to fail? */
  if (failure_desired) {
    /* Did it not fail? */
    if (err == GSL_SUCCESS) {
      printf("failed (succeeded)\n");
      (void)(*num_failed)++;
    }
    /* It failed. */
    else {
      printf("ok (error %d)\n", gsl_errno);
      (void)(*num_passed)++;
    }
  }
  /* We were supposed to succeed. */
  else {
    /* Was there was an error? */
    if (err != GSL_SUCCESS) {
      printf("failed (error %d, %s)\n", gsl_errno, g_reason);
      (void)(*num_failed)++;
    }
    /* Was it not accurate enough? */
    else if (!_WITHIN_TOL(root, cor_root, TEST_REL_EPSILON,
                          TEST_ABS_EPSILON)) {
      printf("failed (inaccurate)\n");
      (void)(*num_failed)++;
    }
    /* The test passed. */
    else {
      printf("ok\n");
      (void)(*num_passed)++;
    }
  }
}

/* f(x) = x^{20} - 1 */
/* f'(x) = 20x^{19} */
/* zero at x = 1 or -1 */
double
gsl_test_root_hairy_1(double x)
{
  return pow(x, 20.0) - 1;
}

void
gsl_test_root_hairy_1_fdf(double * y, double * yprime, double x, int y_wanted,
                          int yprime_wanted)
{
  *y = gsl_test_root_hairy_1(x);
  *yprime = 20.0 * pow(x, 19.0);
}

/* f(x) = sqrt(abs(x))*sgn(x) */
/* f'(x) = 1 / sqrt(abs(x) */
/* zero at x = 0 */
double
gsl_test_root_hairy_2(double x)
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
gsl_test_root_hairy_2_fdf(double * y, double * yprime, double x, int y_wanted,
                          int yprime_wanted)
{
  *y = gsl_test_root_hairy_2(x);
  *yprime = 1 / sqrt(fabs(x));
}


/* f(x) = x^2 - 1e-8 */
/* f'(x) = 2x */
/* zero at x = sqrt(1e-8) or -sqrt(1e-8) */
double
gsl_test_root_hairy_3(double x)
{
  return pow(x, 2.0) - 1e-8;
}

void
gsl_test_root_hairy_3_fdf(double * y, double * yprime, double x, int y_wanted,
                          int yprime_wanted)
{
  *y = gsl_test_root_hairy_3(x);
  *yprime = 2 * x;
}

/* f(x) = x exp(-x) */
/* f'(x) = exp(-x) - x exp(-x) */
/* zero at x = 0 */
double
gsl_test_root_hairy_4(double x)
{
  return x * exp(-x);
}

void
gsl_test_root_hairy_4_fdf(double * y, double * yprime, double x, int y_wanted,
                          int yprime_wanted)
{
  *y = gsl_test_root_hairy_4(x);
  *yprime = exp(-x) - x * exp(-x);
}

/* f(x) = 1/(1+exp(x)) */
/* f'(x) = -exp(x) / (1 + exp(x))^2 */
/* no roots! */
double
gsl_test_root_hairy_5(double x)
{
  return 1 / (1 + exp(x));
}

void
gsl_test_root_hairy_5_fdf(double * y, double * yprime, double x, int y_wanted,
                          int yprime_wanted)
{
  *y = gsl_test_root_hairy_5(x);
  *yprime = -exp(x) / pow(1 + exp(x), 2.0);
}

/* f(x) = (x - 1)^7 */
/* f'(x) = 7 * (x - 1)^6 */
/* zero at x = 1 */
double
gsl_test_root_hairy_6(double x)
{
  return pow(x - 1, 7.0);
}

void
gsl_test_root_hairy_6_fdf(double * y, double * yprime, double x, int y_wanted,
                          int yprime_wanted)
{
  *y = gsl_test_root_hairy_6(x);
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

