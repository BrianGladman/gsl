#include <sys/ieeefp.h>
#include <floatingpoint.h>
#include <signal.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_errno.h>

int
gsl_ieee_set_mode (int precision, int rounding, int exception_mask)
{
  char * out ;

  switch (precision)
    {
    case GSL_IEEE_SINGLE_PRECISION:
      ieee_flags ("set", "precision", "single", out) ;
      break ;
    case GSL_IEEE_DOUBLE_PRECISION:
      ieee_flags ("set", "precision", "double", out) ;
      break ;
    case GSL_IEEE_EXTENDED_PRECISION:
      ieee_flags ("set", "precision", "extended", out) ;
      break ;
    default:
      ieee_flags ("set", "precision", "extended", out) ;
    }

  switch (rounding)
    {
    case GSL_IEEE_ROUND_TO_NEAREST:
      ieee_flags ("set", "direction", "nearest", out) ;
      break ;
    case GSL_IEEE_ROUND_DOWN:
      ieee_flags ("set", "direction", "negative", out) ;
      break ;
    case GSL_IEEE_ROUND_UP:
      ieee_flags ("set", "direction", "positive", out) ;
      break ;
    case GSL_IEEE_ROUND_TO_ZERO:
      ieee_flags ("set", "direction", "tozero", out) ;
      break ;
    default:
      ieee_flags ("set", "direction", "nearest", out) ;
    }

  if (exception_mask & GSL_IEEE_MASK_INVALID)
    {
      ieee_handler ("set", "invalid", SIGFPE_IGNORE) ;
    }
  else 
    {
      ieee_handler ("set", "invalid", SIGFPE_ABORT) ;
    }

  if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
    {
      ieee_handler ("set", "denormalized", SIGFPE_IGNORE) ;
    }
  else
    {
      GSL_ERROR ("sunos4 does not support the denormalized operand exception. "
		 "Use 'mask-denormalized' to work around this.",
		 GSL_EUNSUP) ;
    }


  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    {
      ieee_handler ("set", "division", SIGFPE_IGNORE) ;
    } 
  else
    {
      ieee_handler ("set", "division", SIGFPE_ABORT) ;
    }
  
  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    {
      ieee_handler ("set", "overflow", SIGFPE_IGNORE) ;
    }
  else 
    {
      ieee_handler ("set", "overflow", SIGFPE_ABORT) ;
    }

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    {
      ieee_handler ("set", "underflow", SIGFPE_IGNORE) ;
    }
  else
    {
      ieee_handler ("set", "underflow", SIGFPE_ABORT) ;
    }

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    {
      ieee_handler ("set", "inexact", SIGFPE_ABORT) ;
    }
  else
    {
      ieee_handler ("set", "inexact", SIGFPE_IGNORE) ;
    }

  return GSL_SUCCESS ;
}
