#define SAFE_FUNC_CALL(f, x, yp) \
do { \
  *yp = GSL_FN_EVAL(f,x); \
  if (!GSL_IS_REAL(*yp)) \
    GSL_ERROR("function not continuous", GSL_EBADFUNC); \
} while (0)
