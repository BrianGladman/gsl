#ifndef GSL_IEEE_UTILS_H
#define GSL_IEEE_UTILS_H

enum {
  GSL_IEEE_TYPE_NAN = 1,
  GSL_IEEE_TYPE_INF = 2,
  GSL_IEEE_TYPE_NORMAL = 3,
  GSL_IEEE_TYPE_DENORMAL = 4,
  GSL_IEEE_TYPE_ZERO = 5
} ;

typedef struct  {
  int sign ;
  char mantissa[24] ; /* Actual bits are 0..22, element 23 is a NULL */
  int exponent ;
  int type ;
} gsl_ieee_float_rep ;

typedef struct  {
  int sign ;
  char mantissa[53] ; /* Actual bits are 0..51, element 52 is a NULL */
  int exponent ;
  int type ;
} gsl_ieee_double_rep ;


void gsl_ieee_printf_float (const float * x) ;
void gsl_ieee_printf_double (const double * x) ;

void gsl_ieee_float_to_rep (const float * x, gsl_ieee_float_rep * r) ;
void gsl_ieee_double_to_rep (const double * x, gsl_ieee_double_rep * r) ;

#endif /* GSL_IEEE_UTILS_H */

