#include <math.h>
#include <ieeefp.h>
#include <gsl_ieee_utils.h>
#include <gsl_errno.h>

int
gsl_ieee_set_mode (int precision, int rounding, int exception_mask)
{
  fp_except mode = 0 ;

  switch (precision)
    {
      /* FIXME: it may be that it actually only supports DOUBLE, not
         EXTENDED */

    case GSL_IEEE_SINGLE_PRECISION:
      GSL_ERROR ("solaris only supports extended precision rounding (??) "
		 "(single precision is not supported)",
		 GSL_EUNSUP) ;
      break ;
    case GSL_IEEE_DOUBLE_PRECISION:
      GSL_ERROR ("solaris only supports extended precision rounding (??) "
		 "(double precision is not supported)",
		 GSL_EUNSUP) ;
      break ;
    case GSL_IEEE_EXTENDED_PRECISION:
      break ;
    }

  switch (rounding)
    {
    case GSL_IEEE_ROUND_TO_NEAREST:
      fpsetround (FP_RN) ;
      break ;
    case GSL_IEEE_ROUND_DOWN:
      fpsetround (FP_RM) ;
      break ;
    case GSL_IEEE_ROUND_UP:
      fpsetround (FP_RP) ;
      break ;
    case GSL_IEEE_ROUND_TO_ZERO:
      fpsetround (FP_RZ) ;
      break ;
    default:
      fpsetround (FP_RN) ;
    }

  /* Turn on all the exceptions apart from 'inexact' */

  mode = FP_X_INV | FP_X_DZ | FP_X_OFL | FP_X_UFL ;

  if (exception_mask & GSL_IEEE_MASK_INVALID)
    mode &= ~ FP_X_INV ;

  if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
    GSL_ERROR ("solaris does not support the denormalized operand exception",
	       GSL_EUNSUP) ;

  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    mode &= ~ FP_X_DZ ;

  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    mode &= ~ FP_X_OFL ;

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    mode &=  ~ FP_X_UFL ;

  if (exception_mask & GSL_IEEE_CATCH_INEXACT)
    {
      mode |= FP_X_IMP ;
    }
  else
    {
      mode &= ~ FP_X_IMP ;
    }

  fpsetmask (mode) ;

  return GSL_SUCCESS ;

}
