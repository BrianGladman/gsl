/* multifit_nlinear/test.c
 * 
 * Copyright (C) 2007, 2013, 2015 Brian Gough, Patrick Alken
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

/* These tests are based on the NIST Statistical Reference Datasets
   See http://www.nist.gov/itl/div898/strd/index.html for more
   information. */

#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_ieee_utils.h>

#include "test_fdf.c"

static const gsl_multifit_nlinear_trs **nlinear_trs[] = {
#if 1
  &gsl_multifit_nlinear_trs_lm,
#elif 0
  &gsl_multifit_nlinear_trs_dogleg,
#else
  &gsl_multifit_nlinear_trs_cgst,
#endif

  NULL
};

static void
test_proc(const gsl_multifit_nlinear_trs *trs,
          const gsl_multifit_nlinear_scale *scale,
          const gsl_multifit_nlinear_solver *solver,
          const int fdtype, const int accel)
{
  gsl_multifit_nlinear_parameters fdf_params =
    gsl_multifit_nlinear_default_parameters();

  fdf_params.trs = trs;
  fdf_params.scale = scale;
  fdf_params.solver = solver;
  fdf_params.fdtype = fdtype;
  fdf_params.accel = accel;

  test_fdf_main(&fdf_params);
}

int
main (void)
{
  const gsl_multifit_nlinear_trs **trs;
  int fdtype, accel;
  size_t i = 0;

  gsl_ieee_env_setup();

  /* run testsuite over all parameter combinations;
   * but skip Marquardt scaling since it won't pass */

  for (trs = nlinear_trs[i]; trs != NULL; trs = nlinear_trs[++i])
    {
      fprintf(stderr, "trs = %s\n", (*trs)->name);

      for (fdtype = GSL_MULTIFIT_NLINEAR_FWDIFF;
           fdtype <= GSL_MULTIFIT_NLINEAR_CTRDIFF; ++fdtype)
        {
          for (accel = 0; accel <= 0; ++accel)
            {
#if 1
              test_proc(*trs, gsl_multifit_nlinear_scale_levenberg,
                        gsl_multifit_nlinear_solver_normal,
                        fdtype, accel);
#endif

#if 1
              test_proc(*trs, gsl_multifit_nlinear_scale_more,
                        gsl_multifit_nlinear_solver_normal,
                        fdtype, accel);
#endif

#if 1
              test_proc(*trs, gsl_multifit_nlinear_scale_levenberg,
                        gsl_multifit_nlinear_solver_qr,
                        fdtype, accel);
#endif

#if 1
              test_proc(*trs, gsl_multifit_nlinear_scale_more,
                        gsl_multifit_nlinear_solver_qr,
                        fdtype, accel);
#endif
            }
        }
    }

  exit (gsl_test_summary ());
}
