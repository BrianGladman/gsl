#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_errno.h>

void
gsl_ieee_env_setup (void)
{
  const char * p = getenv("GSL_IEEE_MODE") ;

  int precision = 0, rounding = 0, exception_mask = 0 ;

  int comma = 0 ;

  if (p == 0)  /* GSL_IEEE_MODE environment variable is not set */
    return ;

  if (*p == '\0') /* GSL_IEEE_MODE environment variable is empty */
    return ;

  gsl_ieee_read_mode_string (p, &precision, &rounding, &exception_mask) ;

  gsl_ieee_set_mode (precision, rounding, exception_mask) ;
  
  printf("GSL_IEEE_MODE=\"") ;

  /* Print string with a preceeding comma if the list has already begun */

#define PRINTC(x) do {if(comma) printf(","); printf(x); comma++ ;} while(0)
  
  switch (precision) 
    {
    case GSL_IEEE_SINGLE_PRECISION:
      PRINTC("single-precision") ;
      break ;
    case GSL_IEEE_DOUBLE_PRECISION:
      PRINTC("double-precision") ;
      break ;
    case GSL_IEEE_EXTENDED_PRECISION:
      PRINTC("extended-precision") ;
      break ;
    }

  switch (rounding) 
    {
    case GSL_IEEE_ROUND_TO_NEAREST:
      PRINTC("round-to-nearest") ;
      break ;
    case GSL_IEEE_ROUND_DOWN:
      PRINTC("round-down") ;
      break ;
    case GSL_IEEE_ROUND_UP:
      PRINTC("round-up") ;
      break ;
    case GSL_IEEE_ROUND_TO_ZERO:
      PRINTC("round-to-zero") ;
      break ;
    }

  if ((exception_mask & GSL_IEEE_MASK_ALL) == GSL_IEEE_MASK_ALL)
    {
      PRINTC("mask-all") ;
    }
  else if ((exception_mask & GSL_IEEE_MASK_ALL) == 0)
    {
      PRINTC("trap-common") ;
    }
  else 
    {
      if (exception_mask & GSL_IEEE_MASK_INVALID)
	PRINTC("mask-invalid") ;
      
      if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
	PRINTC("mask-denormalized") ;
      
      if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
	PRINTC("mask-division-by-zero") ;
      
      if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
	PRINTC("mask-overflow") ;
      
      if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
	PRINTC("mask-underflow") ;
    }

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    PRINTC("trap-inexact") ;
  
  printf("\"\n") ;
}





