#include <config.h>
#include <string.h>
#include <stdlib.h>
#include <gsl_ieee_utils.h>
#include <gsl_errno.h>

void
gsl_ieee_env_setup (void)
{
  const char * p = getenv("GSL_IEEE_MODE") ;

  int precision = 0, rounding = 0, exception_mask = 0 ;

  if (p == 0)  /* GSL_IEEE_MODE environment variable is not set */
    return ;

  if (strlen(p) == 0) /* GSL_IEEE_MODE environment variable is empty */
    return ;

  gsl_ieee_read_mode_string (p, &precision, &rounding, &exception_mask) ;

  gsl_ieee_set_mode (precision, rounding, exception_mask) ;
  
  printf("GSL_IEEE_MODE=\"") ;
  
  switch (precision) 
    {
    case GSL_IEEE_SINGLE_PRECISION:
      printf("single-precision;") ;
      break ;
    case GSL_IEEE_DOUBLE_PRECISION:
      printf("double-precision;") ;
      break ;
    case GSL_IEEE_EXTENDED_PRECISION:
      printf("extended-precision;") ;
      break ;
    }

  switch (rounding) 
    {
    case GSL_IEEE_ROUND_TO_NEAREST:
      printf("round-to-nearest;") ;
      break ;
    case GSL_IEEE_ROUND_DOWN:
      printf("round-down;") ;
      break ;
    case GSL_IEEE_ROUND_UP:
      printf("round-up;") ;
      break ;
    case GSL_IEEE_ROUND_TO_ZERO:
      printf("round-to-zero;") ;
      break ;
    }

  if ((exception_mask & GSL_IEEE_MASK_ALL) == GSL_IEEE_MASK_ALL)
    {
      printf("mask-all;") ;
    }
  else 
    {
      if (exception_mask & GSL_IEEE_MASK_INVALID)
	printf("mask-invalid;") ;
      
      if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
	printf("mask-denormalized;") ;
      
      if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
	printf("mask-division-by-zero;") ;
      
      if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
	printf("mask-overflow;") ;
      
      if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
	printf("mask-underflow;") ;
    }

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    printf("trap-inexact;") ;
  
  printf("\"\n") ;
}





