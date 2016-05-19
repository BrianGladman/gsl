/* multilarge_nlinear/test.c
 * 
 * Copyright (C) 2015, 2016 Patrick Alken
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
#include <gsl/gsl_multilarge_nlinear.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_ieee_utils.h>

#include "test_fdf.c"

static const gsl_multilarge_nlinear_trs **nlinear_trs[] = {
  &gsl_multilarge_nlinear_trs_cgst,

  NULL
};

static void
test_proc(const gsl_multilarge_nlinear_trs *trs,
          const int fdtype, const int accel)
{
  gsl_multilarge_nlinear_parameters fdf_params =
    gsl_multilarge_nlinear_default_parameters();

  fdf_params.trs = trs;
  fdf_params.fdtype = fdtype;
  fdf_params.accel = accel;

  test_fdf_main(&fdf_params);
}

int
main (void)
{
  const gsl_multilarge_nlinear_trs **trs;
  int fdtype;
  size_t i = 0;

  gsl_ieee_env_setup();

  /* run testsuite over all parameter combinations */

  for (trs = nlinear_trs[i]; trs != NULL; trs = nlinear_trs[++i])
    {
      fprintf(stderr, "trs = %s\n", (*trs)->name);

      for (fdtype = GSL_MULTILARGE_NLINEAR_FWDIFF;
           fdtype <= GSL_MULTILARGE_NLINEAR_CTRDIFF; ++fdtype)
        {
          test_proc(*trs, fdtype, 0);
        }
    }

  exit (gsl_test_summary ());
}
