/* roots.h -- declarations for internal root finding and RF support stuff. */

#ifndef __ROOTS_H__
#define __ROOTS_H__

/* Macros */

#ifndef HAVE_ISINF
#define isinf(x) (x == HUGE_VAL)
#endif /* HAVE_ISINF */

/* Call the pointed-to function with argument x, put its result in y, and barf
   if it returned something icky. */

#define SAFE_FUNC_CALL(f, x, yp) \
do { \
  *yp = GSL_FN_EVAL(f,x); \
  if (!GSL_IS_REAL(*yp)) \
    GSL_ERROR("function not continuous", GSL_EBADFUNC); \
} while (0)

#endif /* __ROOTS_H__ */


