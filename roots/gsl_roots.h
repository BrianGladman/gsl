/* $Id$ */

#ifndef GSL_ROOTS_H
#define GSL_ROOTS_H


/* Macro Constants */

/* Requested epsilon must be this many times greater than DBL_EPSILON to
   protect against roundoff problems. */
#define GSL_ROOT_EPSILON_BUFFER 10.0

/* The minimum allowed value of max_deltay. Somewhat arbitrary. */
#define GSL_ROOT_MIN_MAX_DELTAY 1.0


/* Macros */

#if (HAVE_ISINF == 0)
#define isinf(x) (x == HUGE_VAL)
#endif /* HAVE_ISINF == 1 */

/* Return nonzero if x is a real number, i.e. non NaN or infinite. */
/* FIXME: Is this correct way to check if something is real? */
#define GSL_ISREAL(x) (((x) == (x)) && !isinf (x))

/* Function Prototypes */

int
gsl_root_bisection(double * root, double (* f)(double), double * lower_bound,
                   double * upper_bound, double rel_epsilon,
                   double abs_epsilon, unsigned int max_iterations, 
                   double max_deltay);
int
gsl_root_falsepos(double * root, double (* f)(double), double * lower_bound,
                  double * upper_bound, double rel_epsilon, double abs_epsilon,
                  unsigned int max_iterations, double max_deltay);

int
gsl_root_secant(double * root, double (* f)(double), double * guess1, 
                double * guess2, double rel_epsilon, double abs_epsilon,
                unsigned int max_iterations, double max_step_size);

int
gsl_root_newton(double * root,
                double (* fdf)(double *, double *, double, int, int),
                double * guess, double rel_epsilon, double abs_epsilon,
                unsigned int max_iterations, double max_step_size);

/* Solve for real roots of the standard quadratic equation,
 * returning the number of real roots.
 * x[] is assumed big enough.
 * Roots are returned ordered.
 */
int gsl_root_solve_quadratic(double a, double b, double c, double x[]);


/* Solve for real roots of the cubic equation
 * x^3 + a x^2 + b x + c = 0, returning the
 * number of real roots.
 * x[] is assumed big enough.
 * Roots are returned ordered.
 */
int gsl_root_solve_cubic(double a, double b, double c, double x[]);


#endif /* GSL_ROOTS_H */
