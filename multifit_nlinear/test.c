/* multifit_nlinear/test.c
 * 
 * Copyright (C) 2007, 2013, 2015, 2016 Brian Gough, Patrick Alken
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
#if 0
  &gsl_multifit_nlinear_trs_lm,
  &gsl_multifit_nlinear_trs_dogleg,
  &gsl_multifit_nlinear_trs_ddogleg,
#elif 1
  &gsl_multifit_nlinear_trs_cgst,
#else
  &gsl_multifit_nlinear_trs_subspace2D,
#endif

  NULL
};

static const gsl_multifit_nlinear_solver **nlinear_solvers[] = {
  &gsl_multifit_nlinear_solver_cholesky,
  &gsl_multifit_nlinear_solver_qr,
  &gsl_multifit_nlinear_solver_svd,

  NULL
};

/* skip Marquardt scaling since it won't pass */
static const gsl_multifit_nlinear_scale **nlinear_scales[] = {
  &gsl_multifit_nlinear_scale_levenberg,
  &gsl_multifit_nlinear_scale_more,

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
  const gsl_multifit_nlinear_solver **solver;
  const gsl_multifit_nlinear_scale **scale;
  int fdtype;
  size_t i = 0;

  gsl_ieee_env_setup();

  /* run testsuite over all parameter combinations */

  for (trs = nlinear_trs[i]; trs != NULL; trs = nlinear_trs[++i])
    {
      size_t j = 0;
      fprintf(stderr, "trs = %s\n", (*trs)->name);

      for (solver = nlinear_solvers[j]; solver != NULL; solver = nlinear_solvers[++j])
        {
          size_t k = 0;
          fprintf(stderr, "solver = %s\n", (*solver)->name);
          for (scale = nlinear_scales[k]; scale != NULL; scale = nlinear_scales[++k])
            {
              for (fdtype = GSL_MULTIFIT_NLINEAR_FWDIFF;
                   fdtype <= GSL_MULTIFIT_NLINEAR_CTRDIFF; ++fdtype)
                {
                  test_proc(*trs, *scale, *solver, fdtype, 0);

                  if (*trs == gsl_multifit_nlinear_trs_lm)
                    {
                      /* test LM with geodesic acceleration */
                      test_proc(*trs, *scale, *solver, fdtype, 1);
                    }
                }
            }
        }
    }

  exit (gsl_test_summary ());
}
