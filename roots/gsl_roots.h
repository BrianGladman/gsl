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

/* Call the pointed-to function with argument x, put its result in y, and
   check if it returned something icky. */
#define GSL_ROOT_FPCALL(f, x, y) \
do { \
  y = (*f)(x); \
  if (!GSL_ISREAL(y)) \
    GSL_ERROR("function under search is not continous", GSL_EBADFUNC); \
} while (0)

/* Return the lesser of the absolute values of its two arguments. */
#define GSL_ROOT_MIN_ABS(x, y) ((fabs(x) < fabs(y)) ? fabs(x) : fabs(y))


/* Function Prototypes */

int
gsl_root_bisection(double * root, double (* f)(double), double * lower_bound,
                   double * upper_bound, double epsilon,
                   unsigned int max_iterations);


#endif /* __GSL_ROOTS_H__ */
