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

    gsl_test (n != 3, "gsl_poly_complex_solve_cubic, three roots, x^3 = 27");
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
	      "gsl_poly_complex_solve_cubic, three roots, (x+3)(x^2-4x+13) = 0");
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
    double z[6 * 2];

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
    double z[8 * 2];

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

    double xa[7] = { 0.16, 0.97, 1.94, 2.74, 3.58, 3.73, 4.70 };
    double ya[7] = { 0.73, 1.11, 1.49, 1.84, 2.30, 2.41, 3.07 };

    double dd_expected[7] = { 7.30000000000000e-01,
      4.69135802469136e-01,
      -4.34737219941284e-02,
      2.68681098870099e-02,
      -3.22937056934996e-03,
      6.12763259971375e-03,
      -6.45402453527083e-03
    };

    double dd[7], coeff[7], work[7];

    gsl_poly_dd_init (dd, xa, ya, 7);

    for (i = 0; i < 7; i++)
      {
	gsl_test_rel (dd[i], dd_expected[i], 1e-10,
		      "divided difference dd[%d]", i);
      }

    for (i = 0; i < 7; i++)
      {
	double y = gsl_poly_dd_eval (dd, xa, 7, xa[i]);
	gsl_test_rel (y, ya[i], 1e-10, "divided difference y[%d]", i);
      }

    gsl_poly_dd_taylor (coeff, 1.5, dd, xa, 7, work);

    for (i = 0; i < 7; i++)
      {
	double y = gsl_poly_eval (coeff, 7, xa[i] - 1.5);
	gsl_test_rel (y, ya[i], 1e-10, "taylor expansion about 1.5 y[%d]", i);
      }
  }

  /* Test polynomial allocation */

  {
    size_t i, N = 6;

    gsl_poly *p = gsl_poly_calloc (N);

    gsl_test (p == 0, "gsl_poly_calloc returns valid pointer");
    gsl_test (p->c == 0,
	      "gsl_poly_calloc returns valid coefficients pointer");
    gsl_test (p->size != N, "gsl_poly_calloc sets the size correctly");

    for (i = 0; i < N; i++)
      {
	gsl_poly_set (p, i, 1.25 + (double) i);
      }

    {
      int status = 0;

      for (i = 0; i < N; i++)
	{
	  if (p->c[i] != 1.25 + (double) i)
	    status = 1;
	}

      gsl_test (status, "gsl_poly_set writes into array correctly");
    }

    {
      int status = 0;

      for (i = 0; i < N; i++)
	{
	  if (gsl_poly_get (p, i) != 1.25 + (double) i)
	    status = 1;
	}

      gsl_test (status, "gsl_poly_get reads from array correctly");
    }

    {
      int status = 0;

      gsl_poly *q = gsl_poly_alloc (N);

      gsl_poly_memcpy (q, p);

      for (i = 0; i <N; i++)
	{
	  if (q->c[i] != p->c[i])
	    status = 1;
	}

      gsl_test (status, "gsl_poly_memcpy copies coefficients correctly");

      gsl_poly_free (q);
    }

    gsl_poly_free (p);
  }

  {
    double x, y;
    gsl_poly *p = gsl_poly_calloc (10);

    gsl_poly_set (p, 0, 1.0);
    gsl_poly_set (p, 1, 0.5);
    gsl_poly_set (p, 2, 0.3);

    x = 0.5;
    y = gsl_poly_eval2 (p, x);

    gsl_test_rel (y, 1 + 0.5 * x + 0.3 * x * x, 10 * GSL_DBL_EPSILON,
		  "gsl_poly_eval2({1, 0.5, 0.3}, 0.5)");

    gsl_poly_free (p);
  }

  {
    double x, y;
    gsl_poly *p = gsl_poly_calloc (7);

    double a[10] = { 1, -1, 1, -1, 1, -1, 1};

    gsl_poly_set_from_array (p, a, 7);

    x = 1.0;
    y = gsl_poly_eval2 (p, x);

    gsl_test_rel (y, 1.0, 10 * GSL_DBL_EPSILON,
		  "gsl_poly_eval2({1,-1, 1, -1, 1, -1, 1, -1, 1}, 1.0)");

    gsl_poly_free (p);
  }

#ifdef DEGREE
  {
    int status = 0;
    gsl_poly *p = gsl_poly_calloc (6);

    gsl_poly_set_all (p, 5, 1.0);
    status = gsl_poly_consistent (p, GSL_DBL_EPSILON);
    gsl_test (status, "gsl_poly_consistent GSL_SUCCESS");
    gsl_poly_set (p, 5, 0.0);
    gsl_poly_set (p, 4, 0.0);
    status = gsl_poly_consistent (p, GSL_DBL_EPSILON);
    gsl_test (!status, "gsl_poly_consistent GSL_FAILURE");

    gsl_poly_free (p);
  }
#endif

#ifdef DEGREE
  {
    gsl_poly *p = gsl_poly_calloc (5);

    gsl_poly_set (p, 0, 1.0);
    gsl_poly_set (p, 1, 1.0);
    gsl_poly_set (p, 2, 1.0);
    gsl_poly_set (p, 3, 1.0);
    gsl_poly_set (p, 4, 1.0);
    gsl_poly_set_degree (p, GSL_DBL_EPSILON);
    gsl_test (p->degree != 4,
	      "gsl_poly_set_degree {1,1,1,1}, GSL_DBL_EPSILON");
    gsl_poly_set (p, 4, 1.0e-20);
    gsl_poly_set_degree (p, GSL_DBL_EPSILON);
    gsl_test (p->degree != 3,
	      "gsl_poly_set_degree {1,1,1,1.0e-20}, GSL_DBL_EPSILON");

    gsl_poly_free (p);
  }
#endif

  {
    size_t i;
    gsl_poly *p = gsl_poly_calloc (6);

    double a[6] = {1, 2, 3, 4, 5, 6};

    gsl_poly_set_from_array (p, a, 6);
    gsl_poly_scale (p, 2.1354);

    for (i = 0; i < 6; i++)
      gsl_test_rel (p->c[i], 2.1354*a[i], 1.0e-09, "x^%u, poly scale", i);

    gsl_poly_free (p);
  }

  {
    size_t i;
    gsl_poly *p1 = gsl_poly_calloc (6);
    gsl_poly *p2 = gsl_poly_calloc (6);

    double a[6] = { 2.18, -1.01, 7.39, 1.69, 5.11, 2.16 };
    double b[6] = { 9.97, 8.13, -1.34, -9.39, 1.37, -8.24};

    gsl_poly_set_from_array (p1, a, 6);
    gsl_poly_set_from_array (p2, b, 6);

    gsl_poly_add (p1, p2);

    for (i = 0; i < 6; i++)
      gsl_test_rel (p1->c[i], a[i]+b[i], 1.0e-09, "x^%u, poly addition", i);

    gsl_poly_free (p1);
    gsl_poly_free (p2);
  }

  {
    size_t i;
    gsl_poly *p1 = gsl_poly_calloc (6);
    gsl_poly *p2 = gsl_poly_calloc (6);

    double a[6] = { 2.18, -1.01, 7.39, 1.69, 5.11, 2.16 };
    double b[6] = { 9.97, 8.13, -1.34, -9.39, 1.37, -8.24};

    gsl_poly_set_from_array (p1, a, 6);
    gsl_poly_set_from_array (p2, b, 6);

    gsl_poly_sub (p1, p2);

    for (i = 0; i < 6; i++)
      gsl_test_rel (p1->c[i], a[i]-b[i], 1.0e-09, "x^%u, poly subtraction", i);

    gsl_poly_free (p1);
    gsl_poly_free (p2);
  }


#ifdef DEGREE
  {
    int status = 0;
    size_t i;
    gsl_poly *p1 = gsl_poly_calloc (6);
    gsl_poly *p2 = gsl_poly_calloc (6);

    double a[6] = { 2.18, -1.01, 7.39, 1.69, 5.11, 2.16 };
    double b[6] = { 9.97, 8.13, -1.34, -9.39, -5.11, -2.16};

    gsl_poly_set_from_array (p1, a, 6);
    gsl_poly_set_from_array (p2, b, 6);

    gsl_poly_add (p1, p2);

    for (i = 0; i < 6; i++)
      gsl_test_rel (p1->c[i], a[i]+b[i], 1.0e-09, "x^%u, poly addition", i);

    gsl_test_abs (p1->c[0], 0.0, 1.0e-09, "x^0, poly addition");
    gsl_test_rel (p1->c[1], 2.0, 1.0e-09, "x^1, poly addition");
    gsl_test_rel (p1->c[2], 2.0, 1.0e-09, "x^2, poly addition");
    gsl_test_rel (p1->c[3], 2.0, 1.0e-09, "x^3, poly addition");
    gsl_test_abs (p1->c[4], 0.0, 1.0e-09, "x^4, poly addition");
    gsl_test_abs (p1->c[5], 0.0, 1.0e-09, "x^5, poly addition");


    gsl_poly_set_degree (p1, GSL_DBL_EPSILON);
    if (p1->degree != 3)
      {
	status = 1;
      }
    gsl_test (status, "gsl_poly_add with lower degree for result");


    gsl_poly_free (p1);
    gsl_poly_free (p2);
  }
#endif

  {
    size_t i;

    gsl_poly *p1 = gsl_poly_calloc (4);
    gsl_poly *p2 = gsl_poly_calloc (3);
    gsl_poly *w = gsl_poly_calloc (6);

    double u[8] = { 1, 1, 1, 1 };
    double v[8] = { 1, 1, 1 };
    double p[8] = { 1, 2, 3, 3, 2, 1 };

    gsl_poly_set_from_array (p1, u, 4);
    gsl_poly_set_from_array (p2, v, 3);

    gsl_poly_mul (w, p1, p2);

    for (i = 0; i < 6; i++)
      gsl_test_rel (w->c[i], p[i], 1.0e-09, "x^%u, poly multiplication", i);

    gsl_poly_free (p1);
    gsl_poly_free (p2);
    gsl_poly_free (w);
  }

  {
    size_t i;

    gsl_poly *p1 = gsl_poly_calloc (6);
    gsl_poly *p2 = gsl_poly_calloc (3);
    gsl_poly *q = gsl_poly_calloc (6);
    gsl_poly *r = gsl_poly_calloc (6);

    double u[6] = { 5.74, -2.11, -3.33, 9.06, 2.40, 4.63};
    double v[3] = { 9.54, -3.28, -6.81 };

    double qok[6] = { 1.54774201200558894e0, -2.27080863126619150e0,
                      -2.49611329973844286e-2, -6.79882525697503671e-1,
                      0.0, 0.0};
    double rok[6] = { -9.02545879453331846e0, 2.46301081416577986e1,
                      0.0, 0.0, 0.0, 0.0};

    gsl_poly_set_from_array (p1, u, 6);
    gsl_poly_set_from_array (p2, v, 3);

    gsl_poly_div (q, r, p1, p2);

    for (i = 0; i < 6; i++)
      gsl_test_rel (q->c[i], qok[i], 1.0e-9, "quotient x^%u, poly division", i);
    for (i = 0; i < 6; i++)    
      gsl_test_rel (r->c[i], rok[i], 1.0e-9, "remainder x^%u, poly division", i);

    gsl_poly_free (p1);
    gsl_poly_free (p2);
    gsl_poly_free (q);
    gsl_poly_free (r);
  }

  {
    size_t i;
    gsl_poly *p1 = gsl_poly_calloc (6);
    gsl_poly *p2 = gsl_poly_calloc (3);
    gsl_poly *q = gsl_poly_calloc (6);
    gsl_poly *r = gsl_poly_calloc (6);

    double u[6] = { 1, 1, 1, 1, 1, 0 };
    double v[3] = { 1, 1, 0 };

    double qok[6] = { 1, 0, 1, 0, 1, 0};
    double rok[6] = { 0, 0, 0, 0, 0, 0};

    gsl_poly_set_from_array (p1, u, 6);
    gsl_poly_set_from_array (p2, v, 3);

    gsl_poly_div (q, r, p1, p2);

    for (i = 0; i < 6; i++)
      gsl_test_rel (q->c[i], qok[i], 1.0e-9, "quotient x^%u, poly division", i);
    for (i = 0; i < 6; i++)    
      gsl_test_rel (r->c[i], rok[i], 1.0e-9, "remainder x^%u, poly division", i);

    gsl_poly_free (p1);
    gsl_poly_free (p2);
    gsl_poly_free (q);
    gsl_poly_free (r);

  }

  {
    size_t i;
    gsl_poly *p = gsl_poly_calloc (6);
    gsl_poly *dp = gsl_poly_calloc (5);

    double c[6] = { 1, 1, 1, 1, 1, 1 };
    double dpok[5] = { 1.0, 2.0, 3.0, 4.0, 5.0 };

    gsl_poly_set_from_array (p, c, 5);

    gsl_poly_diff (dp, p);

    for (i = 0; i < 6; i++)
      gsl_test_rel (dp->c[i], dpok[i], 1.0e-9, "x^%u, poly diff", i);

    gsl_poly_free (p);
    gsl_poly_free (dp);
  }

  /* now summarize the results */

  exit (gsl_test_summary ());
}
