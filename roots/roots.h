/* $Id# */

#ifndef __ROOTS_H__
#define __ROOTS_H__


/* Macros */

/* Call the pointed-to function with argument x, put its result in y, and
   check if it returned something icky. */
#define _GSL_ROOT_FPCALL(f, x, y) \
do { \
  y = (*f)(x); \
  if (!GSL_ISREAL(y)) \
    GSL_ERROR("function under search is not continous", GSL_EBADFUNC); \
} while (0)

/* Calculate the error, taking into account that we switch from relative to
   absolute error if any part of the region of interest is within [-1,1]. */
#define _GSL_ROOT_ERR(a, b) ((fabs(a) < 1.0 || fabs(b) < 1.0) \
     ? fabs((a) - (b)) : fabs((a) - (b)) / _GSL_ROOT_MINA(fabs(a), fabs(b)))

/* Return the minumum absolute value of its two arguments. */
#define _GSL_ROOT_MINA(a, b) ((fabs(a) < fabs(b)) ? fabs(a) : fabs(b))

#endif /* __ROOTS_H__ */
