/* $Id# */

#ifndef __ROOTS_H__
#define __ROOTS_H__


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

#endif /* __ROOTS_H__ */

