/* $Id$ */

#ifndef __GSL_ROOTS_H__
#define __GSL_ROOTS_H__


/* Macro Constants */

/* Requested epsilon must be this many times greater than DBL_EPSILON to
   protect against roundoff problems. */
#define GSL_ROOT_EPSILON_BUFFER 10.0


/* Macros */

/* Return nonzero if x is a real number, i.e. non NaN or infinite. */
/* FIXME: Is this correct way to check if something is real? */
#define GSL_ISREAL(x) (((x) == (x)) && !isinf (x))

  
/* Function Prototypes */

int
gsl_root_bisection(double * root, double (* f)(double), double * lower_bound,
                   double * upper_bound, double rel_epsilon,
                   double abs_epsilon, unsigned int max_iterations);


#endif /* __GSL_ROOTS_H__ */
