/* roots/test_macros.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Reid Priedhorsky, Brian Gough
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

#include <gsl/gsl_test.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "test.h"

/* Test certain macros. */
void
test_macros (void)
{
  int result;
  double inf, nan ;

  /* 1.0 is real */
  result = GSL_IS_REAL (1.0);
  gsl_test (result != 1, "GSL_IS_REAL(1.0) is 1");

  inf = 1.0 / (sqrt(1.0) - 1) ;

  /* 1.0/0.0 == Inf is not real */
  result = GSL_IS_REAL (inf);
  gsl_test (result != 0, "GSL_IS_REAL(Inf) is 0");

  nan = inf - inf ;

  /* 0.0/0.0 == NaN is not real */
  result = GSL_IS_REAL (nan);
  gsl_test (result != 0, "GSL_IS_REAL(NaN) is 0");
}
