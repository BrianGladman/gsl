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
  char mantissa[24] ; /* Actual bits are 0..22, element 23 is \0 */
  int exponent ;
  int type ;
} gsl_ieee_float_rep ;

typedef struct  {
  int sign ;
  char mantissa[53] ; /* Actual bits are 0..51, element 52 is \0 */
  int exponent ;
  int type ;
} gsl_ieee_double_rep ;


void gsl_ieee_printf_float (const float * x) ;
void gsl_ieee_printf_double (const double * x) ;

void gsl_ieee_float_to_rep (const float * x, gsl_ieee_float_rep * r) ;
void gsl_ieee_double_to_rep (const double * x, gsl_ieee_double_rep * r) ;

enum {
  GSL_IEEE_SINGLE_PRECISION = 1,
  GSL_IEEE_DOUBLE_PRECISION = 2,
  GSL_IEEE_EXTENDED_PRECISION = 3
} ;

enum {
  GSL_IEEE_ROUND_TO_NEAREST = 1,
  GSL_IEEE_ROUND_DOWN = 2,
  GSL_IEEE_ROUND_UP = 3,
  GSL_IEEE_ROUND_TO_ZERO = 4
} ;

enum {
  GSL_IEEE_MASK_INVALID = 1,
  GSL_IEEE_MASK_DENORMALIZED = 2,
  GSL_IEEE_MASK_DIVISION_BY_ZERO = 4,
  GSL_IEEE_MASK_OVERFLOW = 8,
  GSL_IEEE_MASK_UNDERFLOW = 16,
  GSL_IEEE_TRAP_INEXACT = 32
} ;

void gsl_ieee_env_setup (void) ;
int gsl_ieee_read_mode_string (const char * description, int * precision,
			       int * rounding, int * exception_mask) ;
int gsl_ieee_set_mode (int precision, int rounding, int exception_mask) ;

#endif /* GSL_IEEE_UTILS_H */

