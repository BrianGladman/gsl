#include <stdio.h>
#include <fpu_control.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_ieee_utils.h>

int
gsl_ieee_set_mode (int precision, int rounding, int exception_mask)
{
  unsigned short mode = 0 ;

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

  /* FIXME: I don't have documentation for the M68K so I'm not sure
     about the mapping of the exceptions below. Maybe someone who does
     know could correct this. */

  if (exception_mask & GSL_IEEE_MASK_INVALID)
    mode |= _FPU_MASK_OPERR ;
  
  if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
    {
      GSL_ERROR ("the denormalized operand exception has not been implemented for m68k linux yet. Use 'mask-denormalized' to work around this.", GSL_EUNSUP) ;
      /*mode |= _FPU_MASK_DM ; ???? */ 
    }
  
  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    mode |= _FPU_MASK_DZ ;

  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    mode |= _FPU_MASK_OVFL ;

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    mode |= _FPU_MASK_UNFL ;

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    {
      mode &= ~ (_FPU_MASK_INEX1 | _FPU_MASK_INEX2) ;
    }
  else
    {
      mode |= (_FPU_MASK_INEX1 | _FPU_MASK_INEX2) ;
    }

  _FPU_SETCW(mode) ;

  return GSL_SUCCESS ;
}
