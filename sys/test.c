/* sys/test.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
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

#include <config.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>

#include <gsl/gsl_ieee_utils.h>

int
main (void)
{
  double y, y_expected;

  gsl_ieee_env_setup ();

  /* Test for expm1 */

  y = gsl_expm1 (0.0); y_expected = 0.0;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_expm1(0.0)");

  y = gsl_expm1 (1e-10); y_expected = 1.000000000050000000002e-10;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_expm1(1e-10)");

  y = gsl_expm1 (-1e-10); y_expected = -9.999999999500000000017e-11;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_expm1(-1e-10)");

  y = gsl_expm1 (0.1); y_expected = 0.1051709180756476248117078264902;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_expm1(0.1)");

  y = gsl_expm1 (-0.1); y_expected = -0.09516258196404042683575094055356;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_expm1(-0.1)");

  y = gsl_expm1 (10.0); y_expected = 22025.465794806716516957900645284;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_expm1(10.0)");

  y = gsl_expm1 (-10.0); y_expected = -0.99995460007023751514846440848444;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_expm1(-10.0)");
   
  /* Test for log1p */

  y = gsl_log1p (0.0); y_expected = 0.0;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_log1p(0.0)");

  y = gsl_log1p (1e-10); y_expected = 9.9999999995000000000333333333308e-11;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_log1p(1e-10)");

  y = gsl_log1p (0.1); y_expected = 0.095310179804324860043952123280765;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_log1p(0.1)");

  y = gsl_log1p (10.0); y_expected = 2.3978952727983705440619435779651;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_log1p(10.0)");

  /* Test for gsl_hypot */

  y = gsl_hypot (0.0, 0.0) ; y_expected = 0.0;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(0.0, 0.0)");

  y = gsl_hypot (1e-10, 1e-10) ; y_expected = 1.414213562373095048801688e-10;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(1e-10, 1e-10)");

  y = gsl_hypot (1e-10, -1.0) ; y_expected = 1.000000000000000000005;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(1e-10, -1)");

  y = gsl_hypot (-1.0, 1e-10) ; y_expected = 1.000000000000000000005;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(-1, 1e-10)");

  y = gsl_hypot (1e307, 1e301) ; y_expected = 1.000000000000499999999999e307;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(1e307, 1e301)");

  y = gsl_hypot (1e301, 1e307) ; y_expected = 1.000000000000499999999999e307;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(1e301, 1e307)");

  y = gsl_hypot (1e307, 1e307) ; y_expected = 1.414213562373095048801688e307;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(1e307, 1e307)");

  return gsl_test_summary ();
}


