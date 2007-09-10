/* roots/demo.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Reid Priedhorsky, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "demof.h"
#include "demof.c"

int
main ()
{
  int status;
  int iterations = 0, max_iterations = 100;
  gsl_root_fsolver *s;
  double r = 0, r_expected = sqrt (5.0);
  double x_lower x = 0.0, x_upper = 5.0;
  gsl_function F;
  struct quadratic_params params =
  {1.0, 0.0, -5.0};

  F.function = &quadratic;
  F.params = &params;

  s = gsl_root_fsolver_alloc (gsl_root_fsolver_bisection);
  gsl_root_fsolver_set (s, &F, x);

  printf ("using %s method\n", gsl_root_fsolver_name (s));

  printf ("%5s [%9s, %9s] %9s %9s %10s %9s\n",
          "iter", "lower", "upper", "root", "actual", "err", "err(est)");

  do
    {
      iterations++;
      status = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      x_lower = gsl_root_fsolver_x_lower (s);
      x_upper = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x, 0, 0.001);

      if (status == GSL_SUCCESS)
        printf ("Converged:\n");

      printf ("%5d [%.7f, %.7f] %.7f %.7f %+.7f %.7f\n",
              iterations, x_lower, x_upper,
              r, r_expected, r - r_expected, x_upper - x_lower);
    }
  while (status == GSL_CONTINUE && iterations < max_iterations);

}
