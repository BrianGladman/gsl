
/*
 * Under Compaq's Unix with the silly name, read the man pages for read_rnd,
 * write_rnd, and ieee(3) for more information on the functions used here.
 *
 * Note that enabling control of dynamic rounding mode (via write_rnd) requires
 * that you pass a special flag to your C compiler.  For Compaq's C compiler
 * the flag is `-fprm d', for gcc it's `-mfp-rounding-mode=d'.
 *
 * Enabling the trap control (via ieee_set_fp_control) also requires a flag be
 * passed to the C compiler.  If you don't need the `inexact' stuff the flag
 * for Compaq's C compiler is `-ieee' and for gcc it's `-mieee'.  If you *do*
 * need the `inexact' stuff, the flag for Compaq's compiler is
 * `-ieee_with_inexact', and the flag for gcc is `-mieee-with-inexact'.
 */

#include <float.h>
#include <machine/fpu.h>
#include <stdio.h>
#include <gsl_ieee_utils.h>
#include <gsl_errno.h>

int
gsl_ieee_set_mode (int precision, int rounding, int exception_mask)
{
  unsigned long int mode = 0 ;
  unsigned int    rnd  = 0 ;

/* I'm actually not completely sure that the alpha only supports default
 * precisions rounding, but I couldn't find any information regarding this, so
 * it seems safe to assume this for now until it's proven otherwise.
 */

  switch (precision)
    {
    case GSL_IEEE_SINGLE_PRECISION:
      GSL_ERROR ("Tru64 Unix on the alpha only supports default precision rounding",
		 GSL_EUNSUP) ;
      break ;
    case GSL_IEEE_DOUBLE_PRECISION:
      GSL_ERROR ("Tru64 Unix on the alpha only supports default precision rounding",
		 GSL_EUNSUP) ;
      break ;
    case GSL_IEEE_EXTENDED_PRECISION:
      GSL_ERROR ("Tru64 Unix on the alpha only supports default precision rounding",
		 GSL_EUNSUP) ;
      break ;
    }


  switch (rounding)
    {
    case GSL_IEEE_ROUND_TO_NEAREST:
      rnd = FP_RND_RN ;
      write_rnd (rnd) ;
      break ;
    case GSL_IEEE_ROUND_DOWN:
      rnd = FP_RND_RM ;
      write_rnd (rnd) ;
      break ;
    case GSL_IEEE_ROUND_UP:
      rnd = FP_RND_RP ;
      write_rnd (rnd) ;
      break ;
    case GSL_IEEE_ROUND_TO_ZERO:
      rnd = FP_RND_RZ ;
      write_rnd (rnd) ;
      break ;
    default:
      rnd = FP_RND_RN ;
      write_rnd (rnd) ;
    }

  /* Turn on all the exceptions apart from 'inexact' */

  /* from the ieee(3) man page:
   * IEEE_TRAP_ENABLE_INV	->	Invalid operation
   * IEEE_TRAP_ENABLE_DZE	->	Divide by 0
   * IEEE_TRAP_ENABLE_OVF	->	Overflow
   * IEEE_TRAP_ENABLE_UNF	->	Underflow
   * IEEE_TRAP_ENABLE_INE	->	Inexact (requires special option to C compiler)
   * IEEE_TRAP_ENABLE_DNO	->	denormal operand
   * IEEE_TRAP_ENABLE_MASK	->	mask of all the trap enables
   * IEEE_MAP_DMZ			->	map denormal inputs to zero
   * IEEE_MAP_UMZ			->	map underflow results to zero
   */

  mode = IEEE_TRAP_ENABLE_INV | IEEE_TRAP_ENABLE_DZE | IEEE_TRAP_ENABLE_OVF
		| IEEE_TRAP_ENABLE_UNF ;

  if (exception_mask & GSL_IEEE_MASK_INVALID)
    mode &= ~ IEEE_TRAP_ENABLE_INV ;

  if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
    mode &= ~ IEEE_TRAP_ENABLE_DNO ;

  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    mode &= ~ IEEE_TRAP_ENABLE_DZE ;

  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    mode &= ~ IEEE_TRAP_ENABLE_OVF ;

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    mode &=  ~ IEEE_TRAP_ENABLE_UNF ;

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    {
	/* requires special flag to C compiler */
      mode |= IEEE_TRAP_ENABLE_INE ;
    }
  else
    {
	/* requires special flag to C compiler */
      mode &= ~ IEEE_TRAP_ENABLE_INE ;
    }

  ieee_set_fp_control (mode) ;

  return GSL_SUCCESS ;
}
