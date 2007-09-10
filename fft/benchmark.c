/* fft/benchmark.c
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
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <time.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_fft_complex.h>

#include <gsl/gsl_errno.h>

#include "complex_internal.h"
#include "bitreverse.h"
#include "urand.c"

void my_error_handler (const char *reason, const char *file,
                       int line, int err);

int
main (int argc, char *argv[])
{
  double *data, *fft_data;
  gsl_fft_complex_wavetable * cw;
  unsigned int i, logn;
  int result;
  int status;
  clock_t start, end;
  int resolution = CLOCKS_PER_SEC;
  unsigned int n = 1;

  gsl_set_error_handler (&my_error_handler);

  if (argc == 2)
    {
      n = strtol (argv[1], NULL, 0);
    }
  else
    {
      printf ("test n\n");
      exit (EXIT_FAILURE);
    }

  cw = gsl_fft_complex_wavetable_alloc (n);
  cwork = gsl_fft_complex_wavetable_workspace_alloc (n);
  status = gsl_fft_complex_init (n, cw, cwork);
  status = gsl_fft_complex_generate (n, cw, cwork);

  data = (double *) malloc (n * 2 * sizeof (double));
  fft_data = (double *) malloc (n * 2 * sizeof (double));

  for (i = 0; i < n; i++)
    {
      REAL(data,1,i) = urand ();
      IMAG(data,1,i) = urand ();
    }


  /* compute the fft */

  memcpy (fft_data, data, n * 2 * sizeof (double));

  start = clock ();
  i = 0;
  do
    {
      status = gsl_fft_complex_forward (fft_data, 1, n, cw, cwork);
      i++;
      end = clock ();
    }
  while (end < start + resolution && status == 0);

  if (status == 0)
    {
      printf ("n = %d gsl_fft_complex_forward %f seconds\n", n, (end - start) / ((double) i) / ((double) CLOCKS_PER_SEC));
    }
  else
    {
      printf ("MR fft failed\n");
    }


  /* compute the fft with radix2 */
  memcpy (fft_data, data, n * 2 * sizeof (double));

  result = fft_binary_logn(n) ;

  if (result == -1) {
    exit(EXIT_SUCCESS) ;
  } else {
    logn = result ;
  }

  start = clock ();
  i = 0;
  do
    {
      status = gsl_fft_complex_radix2_forward (fft_data, 1, n);
      i++;
      end = clock ();
    }
  while (end < start + resolution && status == 0);

  if (status == 0)
    {
      printf ("n = %d gsl_fft_complex_radix2_forward %f seconds\n", n, (end - start) / ((double) i) / ((double) CLOCKS_PER_SEC));
    }
  else
    {
      printf ("fft_radix2: not a power of 2\n");
    }

  start = clock ();
  i = 0;
  do
    {
      status = fft_complex_bitreverse_order (fft_data, 1, n, logn);
      i++;
      end = clock ();
    }
  while (end < start + resolution && status == 0);

  printf ("n = %d gsl_fft_complex_bitreverse_order %f seconds\n", n, (end - start) / ((double) i) / ((double) CLOCKS_PER_SEC));


  return 0;
}


void
my_error_handler (const char *reason, const char *file, int line, int err)
{
  printf ("error: %s in %s at %d (gsl_errno=%d)\n", reason, file, line, err);
}
