/* poly/test.c
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
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_poly.h>

int
main (void)
{
  const double eps = 100.0 * GSL_DBL_EPSILON;

  gsl_ieee_env_setup ();

  {
    double x, y;
    double c[3] = { 1.0, 0.5, 0.3 };
    x = 0.5;
    y = gsl_poly_eval (c, 3, x);
    gsl_test_rel (y, 1 + 0.5 * x + 0.3 * x * x, eps,
		  "gsl_poly_eval({1, 0.5, 0.3}, 0.5)");
  }

  {
    double x, y;
    double d[11] = { 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1 };
    x = 1.0;
    y = gsl_poly_eval (d, 11, x);
    gsl_test_rel (y, 1.0, eps,
		  "gsl_poly_eval({1,-1, 1, -1, 1, -1, 1, -1, 1, -1, 1}, 1.0)");

  }

  /* Quadratic */

  {
    double x0, x1;

    int n = gsl_poly_solve_quadratic (4.0, -20.0, 26.0, &x0, &x1);

    gsl_test (n != 0, "gsl_poly_solve_quadratic, no roots, (2x - 5)^2 = -1");
  }

  {
    double x0, x1;

    int n = gsl_poly_solve_quadratic (4.0, -20.0, 25.0, &x0, &x1);

    gsl_test (n != 2, "gsl_poly_solve_quadratic, one root, (2x - 5)^2 = 0");
    gsl_test_rel (x0, 2.5, 1e-9, "x0, (2x - 5)^2 = 0");
    gsl_test_rel (x1, 2.5, 1e-9, "x1, (2x - 5)^2 = 0");
    gsl_test (x0 != x1, "x0 == x1, (2x - 5)^2 = 0");
  }

  {
    double x0, x1;

    int n = gsl_poly_solve_quadratic (4.0, -20.0, 21.0, &x0, &x1);

    gsl_test (n != 2, "gsl_poly_solve_quadratic, two roots, (2x - 5)^2 = 4");
    gsl_test_rel (x0, 1.5, 1e-9, "x0, (2x - 5)^2 = 4");
    gsl_test_rel (x1, 3.5, 1e-9, "x1, (2x - 5)^2 = 4");
  }

  {
    double x0, x1;

    int n = gsl_poly_solve_quadratic (4.0, 7.0, 0.0, &x0, &x1);

    gsl_test (n != 2, "gsl_poly_solve_quadratic, two roots, x(4x + 7) = 0");
    gsl_test_rel (x0, -1.75, 1e-9, "x0, x(4x + 7) = 0");
    gsl_test_rel (x1, 0.0, 1e-9, "x1, x(4x + 7) = 0");
  }

  {
    double x0, x1;

    int n = gsl_poly_solve_quadratic (5.0, 0.0, -20.0, &x0, &x1);

    gsl_test (n != 2,
	      "gsl_poly_solve_quadratic, two roots b = 0, 5 x^2 = 20");
    gsl_test_rel (x0, -2.0, 1e-9, "x0, 5 x^2 = 20");
    gsl_test_rel (x1, 2.0, 1e-9, "x1, 5 x^2 = 20");
  }

  /* Cubic */

  {
    double x0, x1, x2;

    int n = gsl_poly_solve_cubic (0.0, 0.0, -27.0, &x0, &x1, &x2);

    gsl_test (n != 1, "gsl_poly_solve_cubic, one root, x^3 = 27");
    gsl_test_rel (x0, 3.0, 1e-9, "x0, x^3 = 27");
  }

  {
    double x0, x1, x2;

    int n = gsl_poly_solve_cubic (-51.0, 867.0, -4913.0, &x0, &x1, &x2);

    gsl_test (n != 3, "gsl_poly_solve_cubic, three roots, (x-17)^3=0");
    gsl_test_rel (x0, 17.0, 1e-9, "x0, (x-17)^3=0");
    gsl_test_rel (x1, 17.0, 1e-9, "x1, (x-17)^3=0");
    gsl_test_rel (x2, 17.0, 1e-9, "x2, (x-17)^3=0");
  }

  {
    double x0, x1, x2;

    int n = gsl_poly_solve_cubic (-57.0, 1071.0, -6647.0, &x0, &x1, &x2);

    gsl_test (n != 3,
	      "gsl_poly_solve_cubic, three roots, (x-17)(x-17)(x-23)=0");
    gsl_test_rel (x0, 17.0, 1e-9, "x0, (x-17)(x-17)(x-23)=0");
    gsl_test_rel (x1, 17.0, 1e-9, "x1, (x-17)(x-17)(x-23)=0");
    gsl_test_rel (x2, 23.0, 1e-9, "x2, (x-17)(x-17)(x-23)=0");
  }

  {
    double x0, x1, x2;

    int n = gsl_poly_solve_cubic (-11.0, -493.0, +6647.0, &x0, &x1, &x2);

    gsl_test (n != 3,
	      "gsl_poly_solve_cubic, three roots, (x+23)(x-17)(x-17)=0");
    gsl_test_rel (x0, -23.0, 1e-9, "x0, (x+23)(x-17)(x-17)=0");
    gsl_test_rel (x1, 17.0, 1e-9, "x1, (x+23)(x-17)(x-17)=0");
    gsl_test_rel (x2, 17.0, 1e-9, "x2, (x+23)(x-17)(x-17)=0");
  }

  {
    double x0, x1, x2;

    int n = gsl_poly_solve_cubic (-143.0, 5087.0, -50065.0, &x0, &x1, &x2);

    gsl_test (n != 3,
	      "gsl_poly_solve_cubic, three roots, (x-17)(x-31)(x-95)=0");
    gsl_test_rel (x0, 17.0, 1e-9, "x0, (x-17)(x-31)(x-95)=0");
    gsl_test_rel (x1, 31.0, 1e-9, "x1, (x-17)(x-31)(x-95)=0");
    gsl_test_rel (x2, 95.0, 1e-9, "x2, (x-17)(x-31)(x-95)=0");
  }

  {
    double x0, x1, x2;

    int n = gsl_poly_solve_cubic (-109.0, 803.0, 50065.0, &x0, &x1, &x2);

    gsl_test (n != 3,
	      "gsl_poly_solve_cubic, three roots, (x+17)(x-31)(x-95)=0");
    gsl_test_rel (x0, -17.0, 1e-9, "x0, (x+17)(x-31)(x-95)=0");
    gsl_test_rel (x1, 31.0, 1e-9, "x1, (x+17)(x-31)(x-95)=0");
    gsl_test_rel (x2, 95.0, 1e-9, "x2, (x+17)(x-31)(x-95)=0");
  }

  /* Quadratic with complex roots */

  {
    gsl_complex z0, z1;

    int n = gsl_poly_complex_solve_quadratic (4.0, -20.0, 26.0, &z0, &z1);

    gsl_test (n != 2,
	      "gsl_poly_complex_solve_quadratic, 2 roots (2x - 5)^2 = -1");
    gsl_test_rel (GSL_REAL (z0), 2.5, 1e-9, "z0.real, (2x - 5)^2 = -1");
    gsl_test_rel (GSL_IMAG (z0), -0.5, 1e-9, "z0.imag, (2x - 5)^2 = -1");

    gsl_test_rel (GSL_REAL (z1), 2.5, 1e-9, "z1.real, (2x - 5)^2 = -1");
    gsl_test_rel (GSL_IMAG (z1), 0.5, 1e-9, "z1.imag, (2x - 5)^2 = -1");
  }

  {
    gsl_complex z0, z1;

    int n = gsl_poly_complex_solve_quadratic (4.0, -20.0, 25.0, &z0, &z1);

    gsl_test (n != 2,
	      "gsl_poly_complex_solve_quadratic, one root, (2x - 5)^2 = 0");
    gsl_test_rel (GSL_REAL (z0), 2.5, 1e-9, "z0.real, (2x - 5)^2 = 0");
    gsl_test_rel (GSL_IMAG (z0), 0.0, 1e-9, "z0.imag (2x - 5)^2 = 0");
    gsl_test_rel (GSL_REAL (z1), 2.5, 1e-9, "z1.real, (2x - 5)^2 = 0");
    gsl_test_rel (GSL_IMAG (z1), 0.0, 1e-9, "z1.imag (2x - 5)^2 = 0");
    gsl_test (GSL_REAL (z0) != GSL_REAL (z1),
	      "z0.real == z1.real, (2x - 5)^2 = 0");
    gsl_test (GSL_IMAG (z0) != GSL_IMAG (z1),
	      "z0.imag == z1.imag, (2x - 5)^2 = 0");
  }

  {
    gsl_complex z0, z1;

    int n = gsl_poly_complex_solve_quadratic (4.0, -20.0, 21.0, &z0, &z1);

    gsl_test (n != 2,
	      "gsl_poly_complex_solve_quadratic, two roots, (2x - 5)^2 = 4");
    gsl_test_rel (GSL_REAL (z0), 1.5, 1e-9, "z0.real, (2x - 5)^2 = 4");
    gsl_test_rel (GSL_IMAG (z0), 0.0, 1e-9, "z0.imag, (2x - 5)^2 = 4");
    gsl_test_rel (GSL_REAL (z1), 3.5, 1e-9, "z1.real, (2x - 5)^2 = 4");
    gsl_test_rel (GSL_IMAG (z1), 0.0, 1e-9, "z1.imag, (2x - 5)^2 = 4");
  }

  {
    gsl_complex z0, z1;

    int n = gsl_poly_complex_solve_quadratic (4.0, 7.0, 0.0, &z0, &z1);

    gsl_test (n != 2,
	      "gsl_poly_complex_solve_quadratic, two roots, x(4x + 7) = 0");
    gsl_test_rel (GSL_REAL (z0), -1.75, 1e-9, "z0.real, x(4x + 7) = 0");
    gsl_test_rel (GSL_IMAG (z0), 0.0, 1e-9, "z0.imag, x(4x + 7) = 0");
    gsl_test_rel (GSL_REAL (z1), 0.0, 1e-9, "z1.real, x(4x + 7) = 0");
    gsl_test_rel (GSL_IMAG (z1), 0.0, 1e-9, "z1.imag, x(4x + 7) = 0");
  }

  {
    gsl_complex z0, z1;

    int n = gsl_poly_complex_solve_quadratic (5.0, 0.0, -20.0, &z0, &z1);

    gsl_test (n != 2,
	      "gsl_poly_complex_solve_quadratic, two roots b = 0, 5 x^2 = 20");
    gsl_test_rel (GSL_REAL (z0), -2.0, 1e-9, "z0.real, 5 x^2 = 20");
    gsl_test_rel (GSL_IMAG (z0), 0.0, 1e-9, "z0.imag, 5 x^2 = 20");
    gsl_test_rel (GSL_REAL (z1), 2.0, 1e-9, "z1.real, 5 x^2 = 20");
    gsl_test_rel (GSL_IMAG (z1), 0.0, 1e-9, "z1.imag, 5 x^2 = 20");
  }

  {
    gsl_complex z0, z1;

    int n = gsl_poly_complex_solve_quadratic (5.0, 0.0, 20.0, &z0, &z1);

    gsl_test (n != 2,
	      "gsl_poly_complex_solve_quadratic, two roots b = 0, 5 x^2 = -20");
    gsl_test_rel (GSL_REAL (z0), 0.0, 1e-9, "z0.real, 5 x^2 = -20");
    gsl_test_rel (GSL_IMAG (z0), -2.0, 1e-9, "z0.imag, 5 x^2 = -20");
    gsl_test_rel (GSL_REAL (z1), 0.0, 1e-9, "z1.real, 5 x^2 = -20");
    gsl_test_rel (GSL_IMAG (z1), 2.0, 1e-9, "z1.imag, 5 x^2 = -20");
  }

  /* Cubic with complex roots */

  {
    gsl_complex z0, z1, z2;

    int n = gsl_poly_complex_solve_cubic (0.0, 0.0, -27.0, &z0, &z1, &z2);

    gsl_test (n != 3, "gsl_poly_complex_solve_cubic, three root, x^3 = 27");
    gsl_test_rel (GSL_REAL (z0), -1.5, 1e-9, "z0.real, x^3 = 27");
    gsl_test_rel (GSL_IMAG (z0), -1.5 * sqrt (3.0), 1e-9,
		  "z0.imag, x^3 = 27");
    gsl_test_rel (GSL_REAL (z1), -1.5, 1e-9, "z1.real, x^3 = 27");
    gsl_test_rel (GSL_IMAG (z1), 1.5 * sqrt (3.0), 1e-9, "z1.imag, x^3 = 27");
    gsl_test_rel (GSL_REAL (z2), 3.0, 1e-9, "z2.real, x^3 = 27");
    gsl_test_rel (GSL_IMAG (z2), 0.0, 1e-9, "z2.imag, x^3 = 27");
  }

  {
    gsl_complex z0, z1, z2;

    int n = gsl_poly_complex_solve_cubic (-1.0, 1.0, 39.0, &z0, &z1, &z2);

    gsl_test (n != 3,
	      "gsl_poly_complex_solve_cubic, three root, (x+3)(x^2-4x+13) = 0");
    gsl_test_rel (GSL_REAL (z0), -3.0, 1e-9, "z0.real, (x+3)(x^2+1) = 0");
    gsl_test_rel (GSL_IMAG (z0), 0.0, 1e-9, "z0.imag, (x+3)(x^2+1) = 0");
    gsl_test_rel (GSL_REAL (z1), 2.0, 1e-9, "z1.real, (x+3)(x^2+1) = 0");
    gsl_test_rel (GSL_IMAG (z1), -3.0, 1e-9, "z1.imag, (x+3)(x^2+1) = 0");
    gsl_test_rel (GSL_REAL (z2), 2.0, 1e-9, "z2.real, (x+3)(x^2+1) = 0");
    gsl_test_rel (GSL_IMAG (z2), 3.0, 1e-9, "z2.imag, (x+3)(x^2+1) = 0");
  }

  {
    gsl_complex z0, z1, z2;

    int n =
      gsl_poly_complex_solve_cubic (-51.0, 867.0, -4913.0, &z0, &z1, &z2);

    gsl_test (n != 3,
	      "gsl_poly_complex_solve_cubic, three roots, (x-17)^3=0");
    gsl_test_rel (GSL_REAL (z0), 17.0, 1e-9, "z0.real, (x-17)^3=0");
    gsl_test_rel (GSL_IMAG (z0), 0.0, 1e-9, "z0.imag, (x-17)^3=0");
    gsl_test_rel (GSL_REAL (z1), 17.0, 1e-9, "z1.real, (x-17)^3=0");
    gsl_test_rel (GSL_IMAG (z1), 0.0, 1e-9, "z1.imag, (x-17)^3=0");
    gsl_test_rel (GSL_REAL (z2), 17.0, 1e-9, "z2.real, (x-17)^3=0");
    gsl_test_rel (GSL_IMAG (z2), 0.0, 1e-9, "z2.imag, (x-17)^3=0");
  }

  {
    gsl_complex z0, z1, z2;

    int n =
      gsl_poly_complex_solve_cubic (-57.0, 1071.0, -6647.0, &z0, &z1, &z2);

    gsl_test (n != 3,
	      "gsl_poly_complex_solve_cubic, three roots, (x-17)(x-17)(x-23)=0");
    gsl_test_rel (GSL_REAL (z0), 17.0, 1e-9, "z0.real, (x-17)(x-17)(x-23)=0");
    gsl_test_rel (GSL_IMAG (z0), 0.0, 1e-9, "z0.imag, (x-17)(x-17)(x-23)=0");
    gsl_test_rel (GSL_REAL (z1), 17.0, 1e-9, "z1.real, (x-17)(x-17)(x-23)=0");
    gsl_test_rel (GSL_IMAG (z1), 0.0, 1e-9, "z1.imag, (x-17)(x-17)(x-23)=0");
    gsl_test_rel (GSL_REAL (z2), 23.0, 1e-9, "z2.real, (x-17)(x-17)(x-23)=0");
    gsl_test_rel (GSL_IMAG (z2), 0.0, 1e-9, "z2.imag, (x-17)(x-17)(x-23)=0");
  }

  {
    gsl_complex z0, z1, z2;

    int n =
      gsl_poly_complex_solve_cubic (-11.0, -493.0, +6647.0, &z0, &z1, &z2);

    gsl_test (n != 3,
	      "gsl_poly_complex_solve_cubic, three roots, (x+23)(x-17)(x-17)=0");
    gsl_test_rel (GSL_REAL (z0), -23.0, 1e-9,
		  "z0.real, (x+23)(x-17)(x-17)=0");
    gsl_test_rel (GSL_IMAG (z0), 0.0, 1e-9, "z0.imag, (x+23)(x-17)(x-17)=0");
    gsl_test_rel (GSL_REAL (z1), 17.0, 1e-9, "z1.real, (x+23)(x-17)(x-17)=0");
    gsl_test_rel (GSL_IMAG (z1), 0.0, 1e-9, "z1.imag, (x+23)(x-17)(x-17)=0");
    gsl_test_rel (GSL_REAL (z2), 17.0, 1e-9, "z2.real, (x+23)(x-17)(x-17)=0");
    gsl_test_rel (GSL_IMAG (z2), 0.0, 1e-9, "z2.imag, (x+23)(x-17)(x-17)=0");
  }


  {
    gsl_complex z0, z1, z2;

    int n =
      gsl_poly_complex_solve_cubic (-143.0, 5087.0, -50065.0, &z0, &z1, &z2);

    gsl_test (n != 3,
	      "gsl_poly_complex_solve_cubic, three roots, (x-17)(x-31)(x-95)=0");
    gsl_test_rel (GSL_REAL (z0), 17.0, 1e-9, "z0.real, (x-17)(x-31)(x-95)=0");
    gsl_test_rel (GSL_IMAG (z0), 0.0, 1e-9, "z0.imag, (x-17)(x-31)(x-95)=0");
    gsl_test_rel (GSL_REAL (z1), 31.0, 1e-9, "z1.real, (x-17)(x-31)(x-95)=0");
    gsl_test_rel (GSL_IMAG (z1), 0.0, 1e-9, "z1.imag, (x-17)(x-31)(x-95)=0");
    gsl_test_rel (GSL_REAL (z2), 95.0, 1e-9, "z2.real, (x-17)(x-31)(x-95)=0");
    gsl_test_rel (GSL_IMAG (z2), 0.0, 1e-9, "z2.imag, (x-17)(x-31)(x-95)=0");
  }


  {
    /* Wilkinson polynomial: y = (x-1)(x-2)(x-3)(x-4)(x-5) */

    double a[6] = { -120, 274, -225, 85, -15, 1 };
    double z[6*2];

    gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc (6);

    int status = gsl_poly_complex_solve (a, 6, w, z);

    gsl_poly_complex_workspace_free (w);

    gsl_test (status,
	      "gsl_poly_complex_solve, 5th-order Wilkinson polynomial");
    gsl_test_rel (z[0], 1.0, 1e-9, "z0.real, 5th-order polynomial");
    gsl_test_rel (z[1], 0.0, 1e-9, "z0.imag, 5th-order polynomial");
    gsl_test_rel (z[2], 2.0, 1e-9, "z1.real, 5th-order polynomial");
    gsl_test_rel (z[3], 0.0, 1e-9, "z1.imag, 5th-order polynomial");
    gsl_test_rel (z[4], 3.0, 1e-9, "z2.real, 5th-order polynomial");
    gsl_test_rel (z[5], 0.0, 1e-9, "z2.imag, 5th-order polynomial");
    gsl_test_rel (z[6], 4.0, 1e-9, "z3.real, 5th-order polynomial");
    gsl_test_rel (z[7], 0.0, 1e-9, "z3.imag, 5th-order polynomial");
    gsl_test_rel (z[8], 5.0, 1e-9, "z4.real, 5th-order polynomial");
    gsl_test_rel (z[9], 0.0, 1e-9, "z4.imag, 5th-order polynomial");
  }

  {
    /* : 8-th order polynomial y = x^8 + x^4 + 1 */

    double a[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    double z[8*2];

    double C = 0.5;
    double S = sqrt (3.0) / 2.0;

    gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc (9);

    int status = gsl_poly_complex_solve (a, 9, w, z);

    gsl_poly_complex_workspace_free (w);

    gsl_test (status, "gsl_poly_complex_solve, 8th-order polynomial");

    gsl_test_rel (z[0], -S, 1e-9, "z0.real, 8th-order polynomial");
    gsl_test_rel (z[1], C, 1e-9, "z0.imag, 8th-order polynomial");
    gsl_test_rel (z[2], -S, 1e-9, "z1.real, 8th-order polynomial");
    gsl_test_rel (z[3], -C, 1e-9, "z1.imag, 8th-order polynomial");
    gsl_test_rel (z[4], -C, 1e-9, "z2.real, 8th-order polynomial");
    gsl_test_rel (z[5], S, 1e-9, "z2.imag, 8th-order polynomial");
    gsl_test_rel (z[6], -C, 1e-9, "z3.real, 8th-order polynomial");
    gsl_test_rel (z[7], -S, 1e-9, "z3.imag, 8th-order polynomial");
    gsl_test_rel (z[8], C, 1e-9, "z4.real, 8th-order polynomial");
    gsl_test_rel (z[9], S, 1e-9, "z4.imag, 8th-order polynomial");
    gsl_test_rel (z[10], C, 1e-9, "z5.real, 8th-order polynomial");
    gsl_test_rel (z[11], -S, 1e-9, "z5.imag, 8th-order polynomial");
    gsl_test_rel (z[12], S, 1e-9, "z6.real, 8th-order polynomial");
    gsl_test_rel (z[13], C, 1e-9, "z6.imag, 8th-order polynomial");
    gsl_test_rel (z[14], S, 1e-9, "z7.real, 8th-order polynomial");
    gsl_test_rel (z[15], -C, 1e-9, "z7.imag, 8th-order polynomial");

  }

  {
    int i;

    double xa[7] = {0.16, 0.97, 1.94, 2.74, 3.58, 3.73, 4.70 };
    double ya[7] = {0.73, 1.11, 1.49, 1.84, 2.30, 2.41, 3.07 };

    double dd_expected[7] = {  7.30000000000000e-01,
                               4.69135802469136e-01,
                              -4.34737219941284e-02,
                               2.68681098870099e-02,
                              -3.22937056934996e-03,
                               6.12763259971375e-03,
                              -6.45402453527083e-03 };

    double dd[7], coeff[7], work[7];
    
    gsl_poly_dd_init (dd, xa, ya, 7);

    for (i = 0; i < 7; i++)
      {
        gsl_test_rel (dd[i], dd_expected[i], 1e-10, "divided difference dd[%d]", i);
      }

    for (i = 0; i < 7; i++)
      {
        double y = gsl_poly_dd_eval(dd, xa, 7, xa[i]);
        gsl_test_rel (y, ya[i], 1e-10, "divided difference y[%d]", i);
      }

    gsl_poly_dd_taylor (coeff, 1.5, dd, xa, 7, work);
    
    for (i = 0; i < 7; i++)
      {
        double y = gsl_poly_eval(coeff, 7, xa[i] - 1.5);
        gsl_test_rel (y, ya[i], 1e-10, "taylor expansion about 1.5 y[%d]", i);
      }
  }

  {
    gsl_poly * p;
    size_t i,j;
    int status = 0;
    gsl_poly * p1 = gsl_poly_calloc(6);
    gsl_poly * p2 = gsl_poly_calloc(10);
    const double eps = 100.0 * GSL_DBL_EPSILON;

    

    for (i = 0; i < 5; i++)
       { 
         p = gsl_poly_calloc(i);
         gsl_test (p == 0, "gsl_poly_calloc returns valid pointer");
         gsl_test (p->c == 0, "gsl_poly_calloc returns valid coefficients pointer");
         for (j = 0; j < i; j++)
            {
              gsl_poly_set(p, j, (double) j);
            }
     
         status = 0;

         for (j = 0; j < i ; j++)
            {
              if (p->c[j] != (double) j)
                status = 1;
            }

         gsl_test (status, "gsl_poly_set writes into array correctly");

         status = 0;

         for (j = 0; j < i ; j++)
            {
              if (gsl_poly_get(p, j) != (double) j)
                status = 1;
            }

         gsl_test (status, "gsl_poly_get reads from array correctly");
          
         gsl_poly_free(p);
       }


    for (i = 0; i <= 5; i++) 
       {
         p = gsl_poly_calloc(i+1);
         gsl_poly_set(p,i,1.0);
         gsl_test (p->degree != i, "gsl_poly sets the degree correctly");
         gsl_poly_free(p);
       }

    for (i = 0; i <= 5; i++)
       {
         gsl_poly_set(p1,i,(double) i);
       }

    gsl_poly_memcpy(p2,p1);

    status = 0;

    if (p2->degree != p1->degree)
      {
        status = 1;
      }
  
    for (i = 0; i <= 5; i++)
       {
         if (p2->c[i] != (double) i)
           {
             status = 1;
           }
       }
  gsl_test(status, "gsl_poly_memcpy");
  gsl_poly_free(p1);
  gsl_poly_free(p2);
 
  gsl_ieee_env_setup ();

  {
    double x, y;
    gsl_poly * p = gsl_poly_calloc(10);

    gsl_poly_set(p,0,1.0);
    gsl_poly_set(p,1,0.5);
    gsl_poly_set(p,2,0.3);
    x = 0.5;
    y = gsl_poly_eval2 (p, x);
    gsl_test_rel (y, 1 + 0.5 * x + 0.3 * x * x, eps,
		  "gsl_poly_eval2({1, 0.5, 0.3}, 0.5)");

    gsl_poly_free(p);
  }

  {
    double x, y;
    gsl_poly * p = gsl_poly_calloc(20);

    gsl_poly_set(p,0,1.0);
    gsl_poly_set(p,1,-1.0);
    gsl_poly_set(p,2,1.0);
    gsl_poly_set(p,3,-1.0);
    gsl_poly_set(p,4,1.0);
    gsl_poly_set(p,5,-1.0);
    gsl_poly_set(p,6,1.0);
    gsl_poly_set(p,7,-1.0);
    gsl_poly_set(p,8,1.0);
    gsl_poly_set(p,9,-1.0);
    gsl_poly_set(p,10,1.0);
    x = 1.0;
    y = gsl_poly_eval2 (p, x);
    gsl_test_rel (y, 1.0, eps,
		  "gsl_poly_eval2({1,-1, 1, -1, 1, -1, 1, -1, 1, -1, 1}, 1.0)");

    gsl_poly_free(p);
  }
 
  {
    int status = 0;
    gsl_poly * p = gsl_poly_calloc(6);

    gsl_poly_set_all(p,5,1.0);
    status = gsl_poly_consistent(p, GSL_DBL_EPSILON);
    gsl_test(status, "gsl_poly_consistent GSL_SUCCESS"); 
    gsl_poly_set(p,5,0.0);
    gsl_poly_set(p,4,0.0);
    status = gsl_poly_consistent(p, GSL_DBL_EPSILON);
    gsl_test(!status, "gsl_poly_consistent GSL_FAILURE"); 
  
    gsl_poly_free(p);
  }

  {
    gsl_poly * p = gsl_poly_calloc(5);

    gsl_poly_set(p,0,1.0);
    gsl_poly_set(p,1,1.0);
    gsl_poly_set(p,2,1.0);
    gsl_poly_set(p,3,1.0);
    gsl_poly_set(p,4,1.0);
    gsl_poly_set_degree(p,GSL_DBL_EPSILON);
    gsl_test(p->degree != 4, "gsl_poly_set_degree {1,1,1,1}, GSL_DBL_EPSILON");
    gsl_poly_set(p,4,1.0e-20);
    gsl_poly_set_degree(p,GSL_DBL_EPSILON);
    gsl_test(p->degree != 3, "gsl_poly_set_degree {1,1,1,1.0e-20}, GSL_DBL_EPSILON");

    gsl_poly_free(p); 
  }

  {
    int status = 0;
    size_t i;
    gsl_poly * p = gsl_poly_calloc(6);

    gsl_poly_set(p,0,1.0); 
    gsl_poly_set(p,1,2.0); 
    gsl_poly_set(p,2,3.0); 
    gsl_poly_set(p,3,4.0); 
    gsl_poly_set(p,4,5.0); 
    gsl_poly_set(p,5,6.0); 
    gsl_poly_scale(p,2.0);
    
    gsl_test_rel(p->c[0], 2.0, 1.0e-09, "x^0, poly scale");
    gsl_test_rel(p->c[1], 4.0, 1.0e-09, "x^1, poly scale");
    gsl_test_rel(p->c[2], 6.0, 1.0e-09, "x^2, poly scale");
    gsl_test_rel(p->c[3], 8.0, 1.0e-09, "x^3, poly scale");
    gsl_test_rel(p->c[4], 10.0, 1.0e-09, "x^4, poly scale");
    gsl_test_rel(p->c[5], 12.0, 1.0e-09, "x^5, poly scale");

    gsl_poly_free(p);
    
  }

  {
    int status = 0;
    size_t i;
    gsl_poly * p1 = gsl_poly_calloc(6);
    gsl_poly * p2 = gsl_poly_calloc(6);

    gsl_poly_set(p1,0,1.0); 
    gsl_poly_set(p1,1,1.0); 
    gsl_poly_set(p1,2,1.0); 
    gsl_poly_set(p1,3,1.0); 
    gsl_poly_set(p2,0,1.0); 
    gsl_poly_set(p2,1,1.0); 
    gsl_poly_set(p2,2,1.0); 
    gsl_poly_set(p2,3,1.0); 
    gsl_poly_set(p2,4,1.0); 
    gsl_poly_set(p2,5,1.0); 
    gsl_poly_add(p1,p2);
    
    gsl_test_rel(p1->c[0], 2.0, 1.0e-09, "x^0, poly addition");
    gsl_test_rel(p1->c[1], 2.0, 1.0e-09, "x^1, poly addition");
    gsl_test_rel(p1->c[2], 2.0, 1.0e-09, "x^2, poly addition");
    gsl_test_rel(p1->c[3], 2.0, 1.0e-09, "x^3, poly addition");
    gsl_test_rel(p1->c[4], 1.0, 1.0e-09, "x^4, poly addition");
    gsl_test_rel(p1->c[5], 1.0, 1.0e-09, "x^5, poly addition");

    gsl_poly_free(p1);
    gsl_poly_free(p2);
    
  }

  {
    int status = 0;
    size_t i;
    gsl_poly * p1 = gsl_poly_calloc(6);
    gsl_poly * p2 = gsl_poly_calloc(6);

    gsl_poly_set(p1,0,1.0); 
    gsl_poly_set(p1,1,1.0); 
    gsl_poly_set(p1,2,1.0); 
    gsl_poly_set(p1,3,1.0); 
    gsl_poly_set(p1,4,1.0); 
    gsl_poly_set(p1,5,1.0); 
    gsl_poly_set(p2,0,-1.0); 
    gsl_poly_set(p2,1,1.0); 
    gsl_poly_set(p2,2,1.0); 
    gsl_poly_set(p2,3,1.0); 
    gsl_poly_set(p2,4,-1.0); 
    gsl_poly_set(p2,5,-1.0); 
    gsl_poly_add(p1,p2);

    gsl_test_abs(p1->c[0], 0.0, 1.0e-09, "x^0, poly addition");
    gsl_test_rel(p1->c[1], 2.0, 1.0e-09, "x^1, poly addition");
    gsl_test_rel(p1->c[2], 2.0, 1.0e-09, "x^2, poly addition");
    gsl_test_rel(p1->c[3], 2.0, 1.0e-09, "x^3, poly addition");
    gsl_test_abs(p1->c[4], 0.0, 1.0e-09, "x^4, poly addition");
    gsl_test_abs(p1->c[5], 0.0, 1.0e-09, "x^5, poly addition");

    gsl_poly_set_degree(p1,GSL_DBL_EPSILON);
    if (p1->degree != 3)
      {
        status = 1;
      }
    gsl_test(status, "gsl_poly_add with lower degree for result");

    gsl_poly_free(p1);
    gsl_poly_free(p2);
    
  }

  {
    int status = 0;
    size_t i;
    gsl_poly * p1 = gsl_poly_calloc(6);
    gsl_poly * p2 = gsl_poly_calloc(6);

    gsl_poly_set(p1,0,1.0); 
    gsl_poly_set(p1,1,1.0); 
    gsl_poly_set(p1,2,1.0); 
    gsl_poly_set(p1,3,1.0); 
    gsl_poly_set(p2,0,2.0); 
    gsl_poly_set(p2,1,2.0); 
    gsl_poly_set(p2,2,2.0); 
    gsl_poly_set(p2,3,2.0); 
    gsl_poly_set(p2,4,2.0); 
    gsl_poly_set(p2,5,2.0); 
    gsl_poly_sub(p1,p2);
    gsl_test_rel(p1->c[0], -1.0, 1.0e-09, "x^0, poly substracttion");
    gsl_test_rel(p1->c[1], -1.0, 1.0e-09, "x^1, poly substraction");
    gsl_test_rel(p1->c[2], -1.0, 1.0e-09, "x^2, poly substraction");
    gsl_test_rel(p1->c[3], -1.0, 1.0e-09, "x^3, poly substraction");
    gsl_test_rel(p1->c[4], -2.0, 1.0e-09, "x^4, poly substraction");
    gsl_test_rel(p1->c[5], -2.0, 1.0e-09, "x^5, poly substraction");

    gsl_poly_free(p1);
    gsl_poly_free(p2);
    
  }

  {
    int status = 0;
    size_t i;
    gsl_poly * p1 = gsl_poly_calloc(6);
    gsl_poly * p2 = gsl_poly_calloc(6);

    gsl_poly_set(p1,0,1.0); 
    gsl_poly_set(p1,1,1.0); 
    gsl_poly_set(p1,2,1.0); 
    gsl_poly_set(p1,3,1.0); 
    gsl_poly_set(p1,4,1.0); 
    gsl_poly_set(p1,5,1.0); 
    gsl_poly_set(p2,0,1.0); 
    gsl_poly_set(p2,1,-1.0); 
    gsl_poly_set(p2,2,-1.0); 
    gsl_poly_set(p2,3,-1.0); 
    gsl_poly_set(p2,4,1.0); 
    gsl_poly_set(p2,5,1.0); 
    gsl_poly_sub(p1,p2);
    gsl_test_abs(p1->c[0], 0.0, 1.0e-09, "x^0, poly substracttion");
    gsl_test_rel(p1->c[1], 2.0, 1.0e-09, "x^1, poly substraction");
    gsl_test_rel(p1->c[2], 2.0, 1.0e-09, "x^2, poly substraction");
    gsl_test_rel(p1->c[3], 2.0, 1.0e-09, "x^3, poly substraction");
    gsl_test_abs(p1->c[4], 0.0, 1.0e-09, "x^4, poly substraction");
    gsl_test_abs(p1->c[5], 0.0, 1.0e-09, "x^5, poly substraction");
    gsl_poly_set_degree(p1,GSL_DBL_EPSILON);
    if (p1->degree != 3)
      {
        status = 1;
      }
    gsl_test(status, "gsl_poly_sub with lower degree for result");
 
    gsl_poly_free(p1);
    gsl_poly_free(p2);
    
  }
 

  {
    gsl_poly * p1 = gsl_poly_calloc(8);
    gsl_poly * p2 = gsl_poly_calloc(4);
    gsl_poly *  w = gsl_poly_calloc(8);

    double p[8] = {1,2,3,4,4,3,2,1};

    gsl_poly_set_all(p1,4,1.0);
    gsl_poly_set_all(p2,3,1.0);
    gsl_poly_mul(p1,p2,w);
    gsl_test_rel(p1->c[0], p[0], 1.0e-09, "x^0, poly multiplication");
    gsl_test_rel(p1->c[1], p[1], 1.0e-09, "x^1, poly multiplication");
    gsl_test_rel(p1->c[2], p[2], 1.0e-09, "x^2, poly multiplication");
    gsl_test_rel(p1->c[3], p[3], 1.0e-09, "x^3, poly multiplication");
    gsl_test_rel(p1->c[4], p[4], 1.0e-09, "x^4, poly multiplication");
    gsl_test_rel(p1->c[5], p[5], 1.0e-09, "x^5, poly multiplication");
    gsl_test_rel(p1->c[6], p[6], 1.0e-09, "x^6, poly multiplication");
    gsl_test_rel(p1->c[7], p[7], 1.0e-09, "x^7, poly multiplication");

    gsl_poly_free(p1);
    gsl_poly_free(p2);
    gsl_poly_free(w);
  }

  {
    double qok[4] = {1,0,0,1};
    double rok[2] = {0,0};
    int status = 0;
    size_t i;
    gsl_poly * p1 = gsl_poly_calloc(6);
    gsl_poly * p2 = gsl_poly_calloc(3);
    gsl_poly *  q = gsl_poly_calloc(6);
    gsl_poly *  r = gsl_poly_calloc(6);

    gsl_poly_set_all(p1,5,1.0);
    gsl_poly_set_all(p2,2,1.0);
    gsl_poly_div(p1,p2,q,r);
    
    gsl_test_rel(q->c[0],qok[0],1.0e-9,"quotient x^0, poly division");
    gsl_test_abs(q->c[1],qok[1],1.0e-9,"quotient x^1, poly division");
    gsl_test_abs(q->c[2],qok[2],1.0e-9,"quotient x^2, poly division");
    gsl_test_rel(q->c[3],qok[3],1.0e-9,"quotient x^3, poly division");
    gsl_test_abs(r->c[0],rok[0],1.0e-9,"remainder x^0, poly division");
    gsl_test_abs(r->c[1],rok[1],1.0e-9,"remainder x^1, poly division");

    gsl_poly_free(p1);
    gsl_poly_free(p2);
    gsl_poly_free(q);
    gsl_poly_free(r);
  }

  {
    double qok[4] = {5.0/16.0,1.0/8.0,1.0/4.0,1.0/2.0};
    double rok[2] = {11.0/16.0,9.0/16.0};
    gsl_poly * p1 = gsl_poly_calloc(6);
    gsl_poly * p2 = gsl_poly_calloc(3);
    gsl_poly *  q = gsl_poly_calloc(6);
    gsl_poly *  r = gsl_poly_calloc(6);

    gsl_poly_set_all(p1,5,1.0);
    gsl_poly_set_all(p2,2,1.0);
    gsl_poly_set(p2,2,2.0);
    gsl_poly_div(p1,p2,q,r);
    gsl_test_rel(q->c[0],qok[0],1.0e-9,"quotient x^0, poly division");
    gsl_test_rel(q->c[1],qok[1],1.0e-9,"quotient x^1, poly division");
    gsl_test_rel(q->c[2],qok[2],1.0e-9,"quotient x^2, poly division");
    gsl_test_rel(q->c[3],qok[3],1.0e-9,"quotient x^3, poly division");
    gsl_test_rel(r->c[0],rok[0],1.0e-9,"remainder x^0, poly division");
    gsl_test_rel(r->c[1],rok[1],1.0e-9,"remainder x^1, poly division");
    gsl_poly_free(p1);
    gsl_poly_free(p2);
    gsl_poly_free(q);
    gsl_poly_free(r);

  }

  {
    double dpok[5] = {1.0,2.0,3.0,4.0,5.0};
    gsl_poly * p  = gsl_poly_calloc(6);
    gsl_poly * dp = gsl_poly_calloc(6);

    gsl_poly_set_all(p,5,1.0);
    gsl_poly_diff(p,dp);
    gsl_test_rel(dp->c[0],dpok[0],1.0e-9,"x^0, poly diff");
    gsl_test_rel(dp->c[1],dpok[1],1.0e-9,"x^1, poly diff");
    gsl_test_rel(dp->c[2],dpok[2],1.0e-9,"x^2, poly diff");
    gsl_test_rel(dp->c[3],dpok[3],1.0e-9,"x^3, poly diff");
    gsl_test_rel(dp->c[4],dpok[4],1.0e-9,"x^4, poly diff");

    gsl_poly_free(p);
    gsl_poly_free(dp);
  }


  {
     gsl_poly * p  = gsl_poly_calloc(6);
     int c1,c2,c3,nr;
     gsl_poly_sturm * pss = gsl_poly_sturm_calloc(6);
     gsl_poly * w = gsl_poly_calloc(6);
     gsl_poly_set(p,0,-1);
     gsl_poly_set(p,1,-3);
     gsl_poly_set(p,5,1);
     gsl_test (pss == 0, "gsl_poly_sturm_calloc returns valid pointer");
     gsl_poly_sturm_build(pss,p,w);
     gsl_test(pss->sturmseq[2]->degree!=1,"Sturm function 2, degree = 1");
     gsl_test(pss->index!=3,"Sturm function index = 3");
     gsl_test_rel(pss->sturmseq[2]->c[0],10.0/24.0,1.0e-9,"Sturm function 2, x^0");
     gsl_test_rel(pss->sturmseq[2]->c[1],1.0,1.0e-9,"Sturm function 2, x^1");
     gsl_test_rel(pss->sturmseq[3]->c[0],5.69859182098765e-01,1.0e-9,"Sturm function 3, x^0");
     gsl_test(pss->sturmseq[3]->degree!=0,"Sturm function 3, degree = 0");
     gsl_test(pss->sturmseq[4]->degree!=0,"Sturm function 4, degree = 0");
     gsl_test(pss->sturmseq[5]->degree!=0,"Sturm function 5, degree = 0");
     c1 = gsl_poly_sturm_changes(pss,-2.0);
     c2 = gsl_poly_sturm_changes(pss, 0.0);
     c3 = gsl_poly_sturm_changes(pss,+2.0);
     gsl_test_rel(c1,3.0,1.0e-9,"Sturm changes in x = -2.0");
     gsl_test_rel(c2,1.0,1.0e-9,"Sturm changes in x =  0.0");
     gsl_test_abs(c3,0.0,1.0e-9,"Sturm changes in x = +2.0");
     nr = gsl_poly_sturm_numroots(pss,-2.0,2.0);
     gsl_test_rel(nr,3,1.0e-9,"Number of roots in [-2,2]");

     free(p); 
     free(w); 
     free(pss);
  }

  {
    gsl_poly_sturm * pss = gsl_poly_sturm_calloc(11);
    gsl_poly       *   p = gsl_poly_calloc(11);
    gsl_poly       *   w = gsl_poly_calloc(11);
    double a = 0.0, b = 0.0;
    int i = 0, nr = 0;
    char buf[200];

    /* Wilkinson poly*/
    gsl_poly_set(p,0,3628800.0);
    gsl_poly_set(p,1,-10628640.0);
    gsl_poly_set(p,2,12753576.0);
    gsl_poly_set(p,3,-8409500.0);
    gsl_poly_set(p,4,3416930.0);
    gsl_poly_set(p,5,-902055.0);
    gsl_poly_set(p,6,157773.0);
    gsl_poly_set(p,7,-18150.0);
    gsl_poly_set(p,8,1320.0);
    gsl_poly_set(p,9,-55.0);
    gsl_poly_set(p,10,1.0);
    gsl_poly_sturm_build(pss,p,w);
    for (i = 0; i < 11; i++) {
      b = ((double) i) + 1.0/2.0;
      nr = gsl_poly_sturm_numroots(pss,a,b);
      sprintf(buf,"Number of roots in [0,%f] Wilkinson p10",b);
      gsl_test_rel(nr,i,1.0e-12,buf);
    }

    gsl_poly_sturm_free(pss);
    gsl_poly_free(p);
    gsl_poly_free(w);
  }

  {

    gsl_poly_sturm * pss = gsl_poly_sturm_alloc(3);
    gsl_poly       * p   = gsl_poly_calloc(3);
    gsl_poly       * w   = gsl_poly_calloc(3);
    int              nr  = 0;

    gsl_poly_set(p,0,1.0);
    gsl_poly_set(p,2,0.825);

    gsl_poly_sturm_build(pss,p,w);
    nr = gsl_poly_sturm_numroots(pss,0.0,1.0);
    gsl_test_abs(nr,0,1.0e-12,"Number of roots of 1+0.825x^2 in [0,1]");   

    gsl_poly_sturm_free(pss);
    gsl_poly_free(p);
    gsl_poly_free(w);
  }

    
 
  {
  /* Test Legendre sets */
    int i = 0;
    double a = 0.0;
    char buf[200];
    gsl_poly * p = gsl_poly_calloc(16);
    gsl_poly * w1= gsl_poly_calloc(16);
    gsl_poly * w2 = gsl_poly_calloc(16);
    gsl_poly * w3 = gsl_poly_calloc(16);

    for (i = 0; i <= 15; i++) {
      gsl_poly_set_Legendre(p,i,w1,w2,w3);
      a = gsl_poly_eval2(p,1.0);
      sprintf(buf,"Legendre sets: value of P%d in 1.0",i);
      gsl_test_rel(a,1.0,1.0e-12,buf);
    }

    gsl_poly_free(p);
    gsl_poly_free(w1);
    gsl_poly_free(w2);
    gsl_poly_free(w3);
  }


  }
  /* now summarize the results */

  exit (gsl_test_summary ());
}
