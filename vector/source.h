/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   e.g. gsl_int_name()  and use BASE as the base datatype      */

#ifndef BASE  /* default to double */
#define BASE double
#define CONCAT2(a,b) a ## _ ## b 
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#else
#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define FUNCTION(a,c) CONCAT3(a,BASE,c)
#define TYPE(dir) CONCAT2(dir,BASE)
#endif
