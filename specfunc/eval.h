/* evaluate a function discarding the status value in a modifiable way */

#define EVAL_RESULT(fn) \
   gsl_sf_result result; \
   int status = fn; \
   if (status != GSL_SUCCESS) { \
     GSL_ERROR(#fn, status); \
   } ; \
   return result.val;
