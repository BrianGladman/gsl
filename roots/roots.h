/* roots.h -- declarations for internal root finding and RF support stuff. */ 
/* $Id# */

#ifndef __ROOTS_H__
#define __ROOTS_H__


/* Constant Macros */

/* Something to indicate clearly that we want to ignore maximum delta-y
   restrictions. */
#define _IGNORE_DELTAY -1.0


/* Macros */

/* Call the pointed-to function with argument x, put its result in y, and barf
   if it returned something icky. */
#define _BARF_FPCALL(f, x, y) \
do { \
  y = (*f)(x); \
  if (!GSL_ISREAL(y)) \
    GSL_ERROR("function under search not continous", GSL_EBADFUNC); \
} while (0)

/* Return the minumum absolute value of its two arguments. */
#define _MINA(a, b) ((fabs(a) < fabs(b)) ? fabs(a) : fabs(b))

/* Return the maximum absolute value of its two arguments. */
#define _MAXA(a, b) ((fabs(a) > fabs(b)) ? fabs(a) : fabs(b))

/* Barf if neither rel_epsilon nor abs_epsilon is meaningful in the context of
   a and b being the bounds of the region of interest. */
#define _BARF_TOLS(a, b, rel_epsilon, abs_epsilon) \
do { \
  if (rel_epsilon * _MINA(a, b) + abs_epsilon \
      < _MAXA(a, b) * DBL_EPSILON * GSL_ROOT_EPSILON_BUFFER) \
    GSL_ERROR("tolerances too small for this context", GSL_ETOL); \
} while (0)

/* Return nonzero if a and b are within tolerance of each other. */
#define _WITHIN_TOL(a, b, rel_epsilon, abs_epsilon) \
     (fabs((a) - (b)) < rel_epsilon * _MINA(a, b) + abs_epsilon)

/* Barf if a and b are not within delta of each other; used to protect against
   finding a discontinouity. FIXME.3: This has problems if delta is too small,
   but GSL's current use of it does not run into this situation. Subtracting b
   from a is not a problem because the signs of a and b differ, making the
   subtraction an addition. */
#define _BARF_DELTAY(a, b, delta) \
do { \
  if (fabs((a) - (b)) > delta) \
    GSL_ERROR("function is probably not continuous", GSL_EBADFUNC); \
} while (0)

/* Barf if a and b are not within delta of each other or if delta is too
   small; used to protect against taking a step size which is too large.
   FIXME.3: This has problems if a and b are nearly equal, but GSL's current
   use of it does not run into this situation. */
#define _BARF_DELTAX(a, b, delta) \
do { \
  if (delta < _MINA((a), (b)) * DBL_EPSILON * GSL_ROOT_EPSILON_BUFFER) \
    GSL_ERROR("maximum step size too small for this context", GSL_ETOL); \
  if (fabs((a) - (b)) > delta) \
    GSL_ERROR("next extrapolation step too large", GSL_ERUNAWAY); \
} while (0)

/* Barf if a is zero; used to protect against dividing by zero. */
#define _BARF_ZERO(a) \
do { \
  if ((a) == 0.0) \
    GSL_ERROR("almost divided by zero", GSL_EZERODIV); \
} while (0)


/* Function Prototypes */

int
_gsl_root_validate_bfp_args(void * root, void * f, double * lower_bound,
                            double * upper_bound, double rel_epsilon,
                            double abs_epsilon, unsigned int max_iterations,
                            double max_deltay);

int
_gsl_root_validate_sm_args(void * root, void * f, double * where1,
                           double * where2, double rel_epsilon,
                           double abs_epsilon, unsigned int max_iterations,
                           double max_step_size);

int
_gsl_root_validate_args(void * root, void * f, double * lower_bound,
                        double * upper_bound, double rel_epsilon,
                        double abs_epsilon, unsigned int max_iterations);

int
_gsl_root_ivt_guar(double (* f)(double), double lower_bound,
                   double upper_bound);

int
_gsl_root_silly_user(double * root, double (* f)(double), double lower_bound,
                     double upper_bound, double rel_epsilon,
                     double abs_epsilon, double max_deltay);

#endif /* __ROOTS_H__ */

