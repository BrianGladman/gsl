/* ieee-utils/fp-darwin.c
 * 
 * Copyright (C) 2001 Rodney Sparapani <rsparapa@mcw.edu>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <float.h>
#include <mach/ppc/exception.h>
/* #include <architecture/ppc/fp_regs.h> */
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_errno.h>

int
gsl_ieee_set_mode (int precision, int rounding, int exception_mask)
{
  int prec = 0 ;
  int mode = 0 ;
  int rnd  = 0 ;

  switch (precision)
    {
    case GSL_IEEE_SINGLE_PRECISION:
      prec = FLT_DIG;
      fpsetprec(prec);      
      break ;
    case GSL_IEEE_DOUBLE_PRECISION:
      prec = DBL_DIG;
      fpsetprec(prec);
      break ;
    case GSL_IEEE_EXTENDED_PRECISION:
      prec = LDBL_DIG;
      fpsetprec(prec);
      break ;
    }

  switch (rounding)
    {
    case GSL_IEEE_ROUND_TO_NEAREST:
      rnd = 0; /* RN_NEAREST */;
      fpsetround (rnd) ;
      break ;
    case GSL_IEEE_ROUND_DOWN:
      rnd = 3; /* RN_TOWARD_MINUS */ ;
      fpsetround (rnd) ;
      break ;
    case GSL_IEEE_ROUND_UP:
      rnd = 2; /* RN_TOWARD_PLUS */;
      fpsetround (rnd) ;
      break ;
    case GSL_IEEE_ROUND_TO_ZERO:
      rnd = 1; /* RN_TOWARD_ZERO */;
      fpsetround (rnd) ;
      break ;
    default:
      rnd = 0; /* RN_NEAREST */;
      fpsetround (rnd) ;
    }

  /* Turn on all the exceptions apart from 'inexact' */

  mode = EXC_PPC_FLT_ZERO_DIVIDE | EXC_PPC_FLT_UNDERFLOW | EXC_PPC_FLT_OVERFLOW | EXC_PPC_FLT_NOT_A_NUMBER;

  if (exception_mask & GSL_IEEE_MASK_INVALID)
    mode &= ~ EXC_PPC_FLT_NOT_A_NUMBER ;

  /* if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
     mode &= ~ FP_X_DNML ; */

  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    mode &= ~ EXC_PPC_FLT_ZERO_DIVIDE ;

  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    mode &= ~ EXC_PPC_FLT_OVERFLOW ;

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    mode &=  ~ EXC_PPC_FLT_UNDERFLOW ;

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    {
      mode |= EXC_PPC_FLT_INEXACT ;
    }
  else
    {
      mode &= ~ EXC_PPC_FLT_INEXACT ;
    }

  fpsetmask (mode) ;

  return GSL_SUCCESS ;

}
