/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   e.g. gsl_int_name()  and use BASE as the base datatype      */

#ifndef BASE  /* default to double */
#define BASE double
#define CONCAT(a,b) a ## _ ## b 
#define FUNCTION(dir,name) CONCAT(dir,name)
#define TYPE(dir) dir
#else
#define CONCAT(a,b) a ## _ ## b 
#define CONCAT(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT1(a,b) CONCAT(a,b)
#define CONCAT1(a,b,c) CONCAT(a,b,c)
#define TYPE(dir) CONCAT1(dir,BASE)
#endif
