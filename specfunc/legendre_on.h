
#define CONCAT2x(a,b)   a ## _ ## b 
#define CONCAT2(a,b)    CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c)  CONCAT3x(a,b,c)

#if defined(LEGENDRE_ARRAY)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define OUTPUT result_array
#define OUTPUT_ARG double result_array[]

#elif defined(LEGENDRE_DERIV_ARRAY)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define OUTPUT result_array, result_deriv_array
#define OUTPUT_ARG double result_array[], double result_deriv_array[]
#define QUALIFIER deriv
#define LEGENDRE_DERIV

#elif defined(LEGENDRE_DERIV_ALT_ARRAY)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define OUTPUT result_array, result_deriv_array
#define OUTPUT_ARG double result_array[], double result_deriv_array[]
#define QUALIFIER deriv_alt
#define LEGENDRE_DERIV
#define LEGENDRE_DERIV_ALT

#elif defined(LEGENDRE_DERIV2_ARRAY)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define OUTPUT result_array, result_deriv_array, result_deriv2_array
#define OUTPUT_ARG double result_array[], double result_deriv_array[], double result_deriv2_array[]
#define QUALIFIER deriv2
#define LEGENDRE_DERIV
#define LEGENDRE_DERIV2

#elif defined(LEGENDRE_DERIV2_ALT_ARRAY)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define OUTPUT result_array, result_deriv_array, result_deriv2_array
#define OUTPUT_ARG double result_array[], double result_deriv_array[], double result_deriv2_array[]
#define QUALIFIER deriv2_alt
#define LEGENDRE_DERIV
#define LEGENDRE_DERIV2
#define LEGENDRE_DERIV_ALT

#endif
