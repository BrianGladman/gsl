#include <config.h>

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <gsl_complex.h>
#include <gsl_errno.h>
#include <gsl_dft_complex.h>
#include <gsl_fft_complex.h>
#include <gsl_fft_real.h>
#include <gsl_fft_halfcomplex.h>
#include <gsl_fft_signals.h>
#include <gsl_test.h>

#include "compare.h"

void check_complex (size_t n) ;

/* Usage: test [n]
   Exercise the fft routines for length n. By default n runs from 1 to 100.
   The exit status indicates success or failure. */

int
main (int argc, char *argv[])
{
  size_t i;
  size_t n = 0;

  if (argc == 2) 
    n = strtol (argv[1], NULL, 0);
 
  if (n)
    {
      check_complex (n) ;
    }
  else
    {
      for (i = 1 ; i < 100 ; i++) 
	{
	  check_complex (i) ;
	}
    }

  return gsl_test_summary ();

}

void check_complex (size_t n) 
{
  size_t i ;
  int status ;

  double *real_data, *real_tmp;
  double *fft_real_data, *fft_real_tmp;
  gsl_complex *complex_data, *complex_tmp;
  gsl_complex *fft_complex_data, *fft_complex_tmp;

  gsl_fft_complex_wavetable * cw ;
  gsl_fft_real_wavetable * rw ;
  gsl_fft_halfcomplex_wavetable * hcw ;

  real_data = malloc (n * sizeof (double));
  complex_data = malloc (n * sizeof (gsl_complex));
  complex_tmp = malloc (n * sizeof (gsl_complex));
  fft_complex_data = malloc (n * sizeof (gsl_complex));
  fft_complex_tmp = malloc (n * sizeof (gsl_complex));

  gsl_set_error_handler (NULL);	/* abort on any errors */

  cw = gsl_fft_complex_wavetable_alloc (n);
  gsl_test (cw == 0, "gsl_fft_complex_wavetable_alloc, n = %d", n);
  
  status = gsl_fft_complex_init (n, cw);
  gsl_test (status, "gsl_fft_complex_wavetable_init, n = %d", n);
  
  /* the wavetable generation is also performed by gsl_fft_complex_init
     but test it here too */
  
  status = gsl_fft_complex_generate_wavetable (n, cw);
  gsl_test (status, "gsl_fft_complex_generate_wavetable, n = %d", n);
  
  complex_data = malloc (n * sizeof (gsl_complex));
  complex_tmp = malloc (n * sizeof (gsl_complex));
  fft_complex_data = malloc (n * sizeof (gsl_complex));
  fft_complex_tmp = malloc (n * sizeof (gsl_complex));
  
  real_data = malloc (n * sizeof (double));
  
  /* mixed radix fft */
  gsl_fft_signal_complex_noise (n, complex_data, fft_complex_data);
  memcpy (complex_tmp, complex_data, n * sizeof (gsl_complex));
  gsl_fft_complex_forward (complex_data, n, cw);
  memcpy (fft_complex_tmp, complex_data, n * sizeof (gsl_complex));
  status = compare_complex_results ("dft", fft_complex_data,
				    "fft of noise", complex_data,
				    n, 1e6);
  gsl_test (status, "gsl_fft_complex_forward with signal_noise, n = %d", n);
  
  /* compute the inverse fft */
  status = gsl_fft_complex_inverse (complex_data, n, cw);
  status = compare_complex_results ("orig", complex_tmp,
				    "fft inverse", complex_data,
				    n, 1e6);
  gsl_test (status, "gsl_fft_complex_inverse with signal_noise, n = %d", n);
  
  /* compute the backward fft */

  status = gsl_fft_complex_backward (fft_complex_tmp, n, cw);

  for (i = 0; i < n; i++)
    {
      complex_tmp[i].real *= n;
      complex_tmp[i].imag *= n;
    }
  status = compare_complex_results ("orig", complex_tmp,
				    "fft backward", fft_complex_tmp,
				    n, 1e6);
  gsl_test (status, "gsl_fft_complex_backward with signal_noise, n = %d", n);
  
  /* mixed radix real fft */
  
  rw = gsl_fft_real_wavetable_alloc (n);
  gsl_test (rw == 0, "gsl_fft_real_wavetable_alloc, n = %d", n);
  
  status = gsl_fft_real_init (n, rw);
  gsl_test (status, "gsl_fft_real_wavetable_init, n = %d", n);
  
  status = gsl_fft_real_generate_wavetable (n, rw);
  gsl_test (status, "gsl_fft_real_generate_wavetable, n = %d", n);
  
  real_data = malloc (n * sizeof (double));
  real_tmp = malloc (n * sizeof (double));
  fft_real_data = malloc (n * sizeof (double));
  fft_real_tmp = malloc (n * sizeof (double));
  
  gsl_fft_signal_real_noise (n, complex_data, fft_complex_data);
  memcpy (complex_tmp, complex_data, n * sizeof (gsl_complex));

  for (i = 0; i < n; i++)
    {
      real_data[i] = complex_data[i].real;
    }
  
  gsl_fft_real (real_data, n, rw);
  gsl_fft_halfcomplex_unpack (real_data, complex_data, n);
  
  status = compare_complex_results ("dft", fft_complex_data,
				    "fft of noise", complex_data,
				    n, 1e6);
  gsl_test (status, "gsl_fft_real with signal_real_noise, n = %d", n);
  
  /* compute the inverse fft */
  

  hcw = gsl_fft_halfcomplex_wavetable_alloc (n);
  gsl_test (hcw == 0, "gsl_fft_halfcomplex_wavetable_alloc, n = %d", n);
  
  status = gsl_fft_halfcomplex_init (n, hcw);
  gsl_test (status, "gsl_fft_halfcomplex_wavetable_init, n = %d", n);
  
  status = gsl_fft_halfcomplex_generate_wavetable (n,
						   hcw);
  gsl_test (status, "gsl_fft_halfcomplex_generate_wavetable, n = %d", n);
  
  status = gsl_fft_halfcomplex (real_data, n, hcw);
  
  for (i = 0; i < n; i++)
    {
      real_data[i] /= n;
    }
  
  gsl_fft_real_unpack (real_data, complex_data, n);
  
  status = compare_complex_results ("orig", complex_tmp,
				    "fft inverse", complex_data,
				    n, 1e6);
  gsl_test (status, "gsl_fft_halfcomplex with data from signal_noise, n = %d", n);
  
  /* pulse */
  gsl_fft_signal_complex_pulse (1, n, 1.0, 0.0, complex_data,
				     fft_complex_data);
  gsl_fft_complex_forward (complex_data, n, cw);
  status = compare_complex_results ("analytic", fft_complex_data,
				    "fft of pulse", complex_data, n, 1e6);
  gsl_test (status, "gsl_fft_complex_forward with signal_pulse, n = %d", n);
  
  gsl_fft_signal_complex_constant (n, 1.0, 0.0, complex_data,
					fft_complex_data);
  gsl_fft_complex_forward (complex_data, n, cw);
  status = compare_complex_results ("analytic", fft_complex_data,
				    "fft of constant", complex_data,
				    n, 1e6);
  gsl_test (status, "gsl_fft_complex_forward with signal_constant, n = %d", n);
  
  status = 0;
  for (i = 0; i < n; i++)
    {
      gsl_fft_signal_complex_exp ((int)i, n, 1.0, 0.0, complex_data,
				       fft_complex_data);
      gsl_fft_complex_forward (complex_data, n, cw);
      status |= compare_complex_results ("analytic", fft_complex_data,
					 "fft of exp", complex_data,
					 n, 1e6);
    };
  gsl_test (status, "gsl_fft_complex_forward with signal_exp, n = %d", n);
  
  gsl_fft_complex_wavetable_free (cw);
  gsl_test (status, "gsl_fft_complex_wavetable_free, n = %d", n);
  
  /* check for memory leaks here if mstats available */
  
  free (complex_data);
  free (complex_tmp);
  free (fft_complex_data);
  free (fft_complex_tmp);
}
