/* specfunc/pow_int.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
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

/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include "gsl_sf_pow_int.h"


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_pow_2(const double x) { return x*x;   }
double gsl_sf_pow_3(const double x) { return x*x*x; }
double gsl_sf_pow_4(const double x) { double x2 = x*x;   return x2*x2;    }
double gsl_sf_pow_5(const double x) { double x2 = x*x;   return x2*x2*x;  }
double gsl_sf_pow_6(const double x) { double x2 = x*x;   return x2*x2*x2; }
double gsl_sf_pow_7(const double x) { double x3 = x*x*x; return x3*x3*x;  }
double gsl_sf_pow_8(const double x) { double x2 = x*x;   double x4 = x2*x2; return x4*x4; }
double gsl_sf_pow_9(const double x) { double x3 = x*x*x; return x3*x3*x3; }


int gsl_sf_pow_int_e(double x, int n, gsl_sf_result * result)
{
  double value = 1.0;
  int count = 0;

  /* CHECK_POINTER(result) */

  if(n < 0) {
    if(x == 0.0) return 0.0; /* FIXME: should be Inf */
    x = 1.0/x;
    n = -n;
  }

  /* repeated squaring method 
   * returns 0.0^0 = 1.0, so continuous in x
   */
  do {
     if(GSL_IS_ODD(n)) value *= x;
     n >>= 1;
     x *= x;
     ++count;
  } while (n);

  result->val = value;
  result->err = 2.0 * GSL_DBL_EPSILON * (count + 1.0) * fabs(value); 

  return GSL_SUCCESS;
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

double gsl_sf_pow_int(const double x, const int n)
{
  EVAL_RESULT(gsl_sf_pow_int_e(x, n, &result));
}
