#include <stdio.h>
#include <fpu_control.h>
#include <gsl_ieee_utils.h>

int
gsl_ieee_set_mode (int precision, int rounding, int exception_mask)
{
  int mode = 0 ;

  switch (precision)
    {
    case GSL_IEEE_SINGLE_PRECISION:
      mode |= _FPU_SINGLE ;
      break ;
    case GSL_IEEE_DOUBLE_PRECISION:
      mode |= _FPU_DOUBLE ;
      break ;
    case GSL_IEEE_EXTENDED_PRECISION:
      mode |= _FPU_EXTENDED ;
      break ;
    default:
      mode |= _FPU_EXTENDED ;
    }

  switch (rounding)
    {
    case GSL_IEEE_ROUND_TO_NEAREST:
      mode |= _FPU_RC_NEAREST ;
      break ;
    case GSL_IEEE_ROUND_DOWN:
      mode |= _FPU_RC_DOWN ;
      break ;
    case GSL_IEEE_ROUND_UP:
      mode |= _FPU_RC_UP ;
      break ;
    case GSL_IEEE_ROUND_TO_ZERO:
      mode |= _FPU_RC_ZERO ;
      break ;
    default:
      mode |= _FPU_RC_NEAREST ;
    }

  if (exception_mask & GSL_IEEE_MASK_INVALID)
    mode |= _FPU_MASK_IM ;

  if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
    mode |= _FPU_MASK_DM ;

  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    mode |= _FPU_MASK_ZM ;

  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    mode |= _FPU_MASK_OM ;

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    mode |= _FPU_MASK_UM ;

  if (exception_mask & GSL_IEEE_CATCH_INEXACT)
    {
      mode &= ~ _FPU_MASK_PM ;
    }
  else
    {
      mode |= _FPU_MASK_PM ;
    }

  printf("p = %d, r = %d, m = %d, mode = %x\n",precision, rounding,
	 exception_mask, mode) ;

  __setfpucw (mode) ;

}
