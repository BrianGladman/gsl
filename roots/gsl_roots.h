/* $Id$ */

#ifndef __GSL_ROOTS_H__
#define __GSL_ROOTS_H__


/* Macro Constants */

/* Requested epsilon must be this many times greater than DBL_EPSILON to
   protect against roundoff problems. */
#define GSL_ROOT_EPSILON_BUFFER 10.0

/* The minimum allowed value of max_deltay. Somewhat arbitrary. */
#define GSL_ROOT_MIN_MAX_DELTAY 1.0


/* Macros */

#ifndef HAVE_ISINF
#define isinf(x) (x == HUGE_VAL)
#endif /* HAVE_ISINF */

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
                void (* fdf)(double *, double *, double, int, int),
                double * guess, double rel_epsilon, double abs_epsilon,
                unsigned int max_iterations, double max_step_size);

#endif /* __GSL_ROOTS_H__ */
