#ifndef GSL_IEEE_UTILS_H
#define GSL_IEEE_UTILS_H

typedef struct  {
  int sign ;
  char bits[23] ;
  int exponent ;
} gsl_ieee_float_rep ;

typedef struct  {
  int sign ;
  char bits[52] ;
  int exponent ;
} gsl_ieee_double_rep ;


void gsl_ieee_printf_float (const float * x) ;
void gsl_ieee_printf_double (const double * x) ;

void gsl_ieee_float_to_rep (const float * x, gsl_ieee_float_rep * r) ;
void gsl_ieee_double_to_rep (const double * x, gsl_ieee_double_rep * r) ;

#endif /* GSL_IEEE_UTILS_H */

