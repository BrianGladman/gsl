/* fft/errs.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
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

#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_dft_complex.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

#include <gsl/gsl_test.h>
#include "compare.h"
#include "complex_internal.h"

#include "urand.c"

int verbose = 0;

size_t tests = 0;
size_t passed = 0;
size_t failed = 0;

int
main (int argc, char *argv[])
{
  int status, factor_sum;
  size_t i, start, end, n;
  double *complex_data, *complex_tmp;
  double rms, total;

  gsl_fft_complex_wavetable * cw;

  if (argc == 2)
    {
      start = strtol (argv[1], NULL, 0);
      end = start + 1;
    }
  else
    {
      start = 1 ;
      end = 1000 ;
    }

  for (n = start; n < end; n++)
    {

      complex_data = (double *) malloc (n * 2 * sizeof (double));
      complex_tmp = (double *) malloc (n * 2 * sizeof (double));

      cw = gsl_fft_complex_wavetable_alloc (n);
      status = gsl_fft_complex_init (n, cw);
      status = gsl_fft_complex_generate (n, cw);

      for (i = 0; i < n; i++)
        {
          REAL(complex_data,1,i) = urand();
          IMAG(complex_data,1,i) = urand();
        }

      memcpy (complex_tmp, complex_data, n * 2 * sizeof (double));
      gsl_fft_complex_forward (complex_data, 1, n, cw);
      gsl_fft_complex_inverse (complex_data, 1, n, cw);

      total = 0.0;
      for (i = 0; i < n; i++)
        {
          double dr = REAL(complex_data,1,i) - REAL(complex_tmp,1,i);
          double di = IMAG(complex_data,1,i) - IMAG(complex_tmp,1,i);
          total += dr * dr + di * di;
        }

      rms = sqrt (total / n);

      factor_sum = 0;
      for (i = 0; i < cw->nf; i++)
        {
          int j = cw->factor[i];
          factor_sum += j;
        }

      printf ("n = %d factor_sum = %d rms = %e\n", n, factor_sum, rms);

      free (complex_data);
      free (complex_tmp);

    }

  return 0;
}


