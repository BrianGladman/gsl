#include <config.h>
#include <stdio.h>
#include <math.h>
#include <gsl_ieee_utils.h>

#include "ieee_utils.h"

/* A table of sign characters, 0=positive, 1=negative. We print a space
   instead of a unary + sign for compatibility with bc */

static char signs[2]={' ','-'} ;  

void 
gsl_ieee_printf_float (const float * x) {
  gsl_ieee_float_rep r ;
  gsl_ieee_float_to_rep(x, &r) ;

  /* output is compatible with bc (with ibase=2), mant*2^expb */  

  switch (r.type)
    {
    case GSL_IEEE_TYPE_NAN:
      printf("NaN") ;
      break ;
    case GSL_IEEE_TYPE_INF:
      printf("%cInf", signs[r.sign]) ;
      break ;
    case GSL_IEEE_TYPE_NORMAL:
      printf("%c1.%s*2^%d", signs[r.sign], r.mantissa, r.exponent) ;
      break ;
    case GSL_IEEE_TYPE_DENORMAL:
      printf("%c0.%s*2^%d", signs[r.sign], r.mantissa, r.exponent + 1) ;
      break ;
    case GSL_IEEE_TYPE_ZERO:
      printf("%c0", signs[r.sign]) ;
      break ;
    default:
      printf("[non-standard IEEE float]") ;
    }
}

void
gsl_ieee_printf_double (const double * x) {
  gsl_ieee_double_rep r ;
  gsl_ieee_double_to_rep (x, &r) ;

  /* output is compatible with bc (with ibase=2), mant*2^expb */  

  switch (r.type)
    {
    case GSL_IEEE_TYPE_NAN:
      printf("NaN") ;
      break ;
    case GSL_IEEE_TYPE_INF:
      printf("%cInf", signs[r.sign]) ;
      break ;
    case GSL_IEEE_TYPE_NORMAL:
      printf("%c1.%s*2^%d", signs[r.sign], r.mantissa, r.exponent) ;
      break ;
    case GSL_IEEE_TYPE_DENORMAL:
      printf("%c0.%s*2^%d", signs[r.sign], r.mantissa, r.exponent + 1) ;
      break ;
    case GSL_IEEE_TYPE_ZERO:
      printf("%c0", signs[r.sign]) ;
      break ;
    default:
      printf("[non-standard IEEE double]") ;
    }
}




