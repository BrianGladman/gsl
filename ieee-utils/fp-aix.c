#include <math.h>
#include <fptrap.h>
#include <float.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_errno.h>

int
gsl_ieee_set_mode (int precision, int rounding, int exception_mask)
{
  fptrap_t   mode = 0 ;
  fprnd_t    rnd  = 0 ;

  switch (precision)
    {

    /* I'm not positive about AIX only supporting default precision rounding,
     * but this is the best assumption until it's proven otherwise. */

    case GSL_IEEE_SINGLE_PRECISION:
      GSL_ERROR ("AIX only supports default precision rounding",
		 GSL_EUNSUP) ;
      break ;
    case GSL_IEEE_DOUBLE_PRECISION:
      GSL_ERROR ("AIX only supports default precision rounding",
		 GSL_EUNSUP) ;
      break ;
    case GSL_IEEE_EXTENDED_PRECISION:
      GSL_ERROR ("AIX only supports default precision rounding",
		 GSL_EUNSUP) ;
      break ;
    }

  switch (rounding)
    {
    case GSL_IEEE_ROUND_TO_NEAREST:
      rnd = FP_RND_RN ;
      fp_swap_rnd (rnd) ;
      break ;
    case GSL_IEEE_ROUND_DOWN:
      rnd = FP_RND_RM ;
      fp_swap_rnd (rnd) ;
      break ;
    case GSL_IEEE_ROUND_UP:
      rnd = FP_RND_RP ;
      fp_swap_rnd (rnd) ;
      break ;
    case GSL_IEEE_ROUND_TO_ZERO:
      rnd = FP_RND_RZ ;
      fp_swap_rnd (rnd) ;
      break ;
    default:
      rnd = FP_RND_RN ;
      fp_swap_rnd (rnd) ;
    }

  /* Turn on all the exceptions apart from 'inexact' */

  mode = TRAP_INVALID | TRAP_DIV_BY_ZERO | TRAP_OVERFLOW | TRAP_UNDERFLOW ;

  if (exception_mask & GSL_IEEE_MASK_INVALID)
    mode &= ~ TRAP_INVALID ;

  if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
    GSL_ERROR ("AIX does not support the denormalized operand exception. "
	       "Use 'mask-denormalized' to work around this.",
	       GSL_EUNSUP) ;

  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    mode &= ~ TRAP_DIV_BY_ZERO ;

  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    mode &= ~ TRAP_OVERFLOW ;

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    mode &=  ~ TRAP_UNDERFLOW ;

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    {
      mode |= TRAP_INEXACT ;
    }
  else
    {
      mode &= ~ TRAP_INEXACT ;
    }

  /* AIX appears to require two steps -- first enable floating point traps
   * in general... */
  fp_trap(FP_TRAP_SYNC);

  /* next, enable the traps we're interested in */
  fp_enable(mode);

  return GSL_SUCCESS ;

}
