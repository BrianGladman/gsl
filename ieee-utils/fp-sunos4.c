#include <stdio.h>
#include <sys/ieeefp.h>
#include <gsl_ieee_utils.h>
#include <gsl_errno.h>

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
    ieee_flags ("set", "exception", "invalid", out) ;

  if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
    GSL_ERROR ("sunos4 does not support the denormalized operand exception",
	       GSL_EUNSUP) ;

  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    ieee_flags ("set", "exception", "division", out) ;

  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    ieee_flags ("set", "exception", "overflow", out) ;

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    ieee_flags ("set", "exception", "underflow", out) ;

  if (exception_mask & GSL_IEEE_CATCH_INEXACT)
    {
      ieee_flags ("set", "exception", "inexact", out) ;
    }

}
