/* roots.h -- declarations for internal root finding and RF support stuff. */

#ifndef __ROOTS_H__
#define __ROOTS_H__

/* Macros */

#ifndef HAVE_ISINF
#define isinf(x) (x == HUGE_VAL)
#endif /* HAVE_ISINF */

/* Call the pointed-to function with argument x, put its result in y, and barf
   if it returned something icky. */

#define SAFE_FUNC_CALL(f, x, y) \
do { \
  y = GSL_FN_EVAL(f,x); \
  if (!GSL_IS_REAL(y)) \
    GSL_ERROR("function not continuous", GSL_EBADFUNC); \
} while (0)

/* Return nonzero if a and b are within tolerance of each other. */

#define WITHIN_TOL(a, b, rel_epsilon, abs_epsilon) \
 (fabs((a) - (b)) < (rel_epsilon) * GSL_MIN(fabs(a), fabs(b)) + (abs_epsilon))

#define CHECK_TOL(a, b, rel_epsilon, abs_epsilon) \
do { \
  if (rel_epsilon * GSL_MIN(fabs(a), fabs(b)) + abs_epsilon \
      < GSL_MAX(fabs(a), fabs(b)) * GSL_DBL_EPSILON * GSL_ROOT_EPSILON_BUFFER)\
    GSL_ERROR("tolerances too small for this context", GSL_EBADTOL); \
} while (0)


#endif /* __ROOTS_H__ */


