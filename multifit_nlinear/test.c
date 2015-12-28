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
#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_ieee_utils.h>

#include "test_fdf.c"

int
main (void)
{
  gsl_multifit_nlinear_parameters fdf_params;

  gsl_ieee_env_setup();

  /* test normal equations method */
  fdf_params.solver = GSL_MULTIFIT_NLINEAR_SOLVER_NORMAL;

  fdf_params.scale = GSL_MULTIFIT_NLINEAR_SCALE_LEVENBERG;
  test_fdf_main(&fdf_params);

  fdf_params.scale = GSL_MULTIFIT_NLINEAR_SCALE_MORE;
  test_fdf_main(&fdf_params);

  /* test QR method */
  fdf_params.solver = GSL_MULTIFIT_NLINEAR_SOLVER_QR;

  fdf_params.scale = GSL_MULTIFIT_NLINEAR_SCALE_LEVENBERG;
  test_fdf_main(&fdf_params);

  fdf_params.scale = GSL_MULTIFIT_NLINEAR_SCALE_MORE;
  test_fdf_main(&fdf_params);

  exit (gsl_test_summary ());
}
