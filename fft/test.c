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

#include <autotest.h>
#include <compare.h>

int verbose = 0;

unsigned int tests = 0;
unsigned int passed = 0;
unsigned int failed = 0;

void check_complex (unsigned int n) ;

int
main (int argc, char *argv[])
{
  int status;
  unsigned int i, n ;

  if (argc == 2)
    {
      n = strtol (argv[1], NULL, 0);
      check_complex (n) ;
    }
  else
    {
      for (i = 1 ; i < 100 ; i++) 
	{
	  check_complex (i) ;
	}
    }

  status = msg_summary (tests, passed, failed);
  
  return status;
}

void check_complex (unsigned int n) 
{
  unsigned int i ;
  int status ;

  double *real_data, *real_tmp;
  double *fft_real_data, *fft_real_tmp;
  complex *complex_data, *complex_tmp;
  complex *fft_complex_data, *fft_complex_tmp;
  char length[256];

  gsl_fft_complex_wavetable complex_wavetable;
  gsl_fft_real_wavetable real_wavetable;
  gsl_fft_halfcomplex_wavetable halfcomplex_wavetable;

  real_data = malloc (n * sizeof (double));
  complex_data = malloc (n * sizeof (complex));
  complex_tmp = malloc (n * sizeof (complex));
  fft_complex_data = malloc (n * sizeof (complex));
  fft_complex_tmp = malloc (n * sizeof (complex));

  gsl_set_error_handler (NULL);	/* abort on any errors */

  sprintf (length, "n = %d", n);
  
  msg_checking_params (length, "gsl_fft_complex_wavetable_alloc");
  status = gsl_fft_complex_wavetable_alloc (n, &complex_wavetable);
  msg_result_status (status);
  
  msg_checking_params (length, "gsl_fft_complex_wavetable_init");
  status = gsl_fft_complex_init (n, &complex_wavetable);
  msg_result_status (status);
  
  /* the wavetable generation is also performed by gsl_fft_complex_init
     but test it here too */
  
  msg_checking_params (length, "gsl_fft_complex_generate_wavetable");
  status = gsl_fft_complex_generate_wavetable (n, &complex_wavetable);
  msg_result_status (status);
  
  complex_data = malloc (n * sizeof (complex));
  complex_tmp = malloc (n * sizeof (complex));
  fft_complex_data = malloc (n * sizeof (complex));
  fft_complex_tmp = malloc (n * sizeof (complex));
  
  real_data = malloc (n * sizeof (double));
  
  /* mixed radix fft */
  msg_checking_params (length, 
		       "gsl_fft_complex_forward with signal_noise");
  gsl_fft_signal_complex_noise (n, complex_data, fft_complex_data);
  memcpy (complex_tmp, complex_data, n * sizeof (complex));
  gsl_fft_complex_forward (complex_data, n, &complex_wavetable);
  memcpy (fft_complex_tmp, complex_data, n * sizeof (complex));
  status = compare_complex_results ("dft", fft_complex_data,
				    "fft of noise", complex_data,
				    n, 1e6);
  msg_result_status (status);
  
  /* compute the inverse fft */
  msg_checking_params (length, 
		       "gsl_fft_complex_inverse with signal_noise");
  status = gsl_fft_complex_inverse (complex_data, n, &complex_wavetable);
  status = compare_complex_results ("orig", complex_tmp,
				    "fft inverse", complex_data,
				    n, 1e6);
  msg_result_status (status);
  
  /* compute the backward fft */
  msg_checking_params (length, 
		       "gsl_fft_complex_backward with signal_noise");
  status = gsl_fft_complex_backward (fft_complex_tmp, n, &complex_wavetable);

  for (i = 0; i < n; i++)
    {
      complex_tmp[i].real *= n;
      complex_tmp[i].imag *= n;
    }
  status = compare_complex_results ("orig", complex_tmp,
				    "fft backward", fft_complex_tmp,
				    n, 1e6);
  msg_result_status (status);
  
  
  /* mixed radix real fft */
  
  msg_checking_params (length, "gsl_fft_real_wavetable_alloc");
  status = gsl_fft_real_wavetable_alloc (n, &real_wavetable);
  msg_result_status (status);
  
  msg_checking_params (length, "gsl_fft_real_wavetable_init");
  status = gsl_fft_real_init (n, &real_wavetable);
  msg_result_status (status);
  
  msg_checking_params (length, "gsl_fft_real_generate_wavetable");
  status = gsl_fft_real_generate_wavetable (n, &real_wavetable);
  msg_result_status (status);
  
  real_data = malloc (n * sizeof (double));
  real_tmp = malloc (n * sizeof (double));
  fft_real_data = malloc (n * sizeof (double));
  fft_real_tmp = malloc (n * sizeof (double));
  
  msg_checking_params (length, "gsl_fft_real with signal_real_noise");
  gsl_fft_signal_real_noise (n, complex_data, fft_complex_data);
  memcpy (complex_tmp, complex_data, n * sizeof (complex));

  for (i = 0; i < n; i++)
    {
      real_data[i] = complex_data[i].real;
    }
  
  gsl_fft_real (real_data, n, &real_wavetable);
  gsl_fft_halfcomplex_unpack (real_data, complex_data, n);
  
  status = compare_complex_results ("dft", fft_complex_data,
				    "fft of noise", complex_data,
				    n, 1e6);
  msg_result_status (status);
  
  /* compute the inverse fft */
  
  msg_checking_params (length, "gsl_fft_halfcomplex_wavetable_alloc");
  status = gsl_fft_halfcomplex_wavetable_alloc (n, &halfcomplex_wavetable);
  msg_result_status (status);
  
  msg_checking_params (length, "gsl_fft_halfcomplex_wavetable_init");
  status = gsl_fft_halfcomplex_init (n, &halfcomplex_wavetable);
  msg_result_status (status);
  
  msg_checking_params (length, "gsl_fft_halfcomplex_generate_wavetable");
  status = gsl_fft_halfcomplex_generate_wavetable (n,
						   &halfcomplex_wavetable);
  msg_result_status (status);
  
  msg_checking_params (length, "gsl_fft_halfcomplex with data from signal_noise");
  status = gsl_fft_halfcomplex (real_data, n, &halfcomplex_wavetable);
  
  for (i = 0; i < n; i++)
    {
      real_data[i] /= n;
    }
  
  gsl_fft_real_unpack (real_data, complex_data, n);
  
  status = compare_complex_results ("orig", complex_tmp,
				    "fft inverse", complex_data,
				    n, 1e6);
  msg_result_status (status);
  
  /* pulse */
  msg_checking_params (length, "gsl_fft_complex_forward with signal_pulse");
  gsl_fft_signal_complex_pulse (1, n, 1.0, 0.0, complex_data,
				     fft_complex_data);
  gsl_fft_complex_forward (complex_data, n, &complex_wavetable);
  status = compare_complex_results ("analytic", fft_complex_data,
				    "fft of pulse", complex_data, n, 1e6);
  msg_result_status (status);
  
  msg_checking_params (length, "gsl_fft_complex_forward with signal_constant");
  gsl_fft_signal_complex_constant (n, 1.0, 0.0, complex_data,
					fft_complex_data);
  gsl_fft_complex_forward (complex_data, n, &complex_wavetable);
  status = compare_complex_results ("analytic", fft_complex_data,
				    "fft of constant", complex_data,
				    n, 1e6);
  msg_result_status (status);
  
  msg_checking_params (length, "gsl_fft_complex_forward with signal_exp");
  status = 0;
  for (i = 0; i < n; i++)
    {
      gsl_fft_signal_complex_exp ((int)i, n, 1.0, 0.0, complex_data,
				       fft_complex_data);
      gsl_fft_complex_forward (complex_data, n, &complex_wavetable);
      status |= compare_complex_results ("analytic", fft_complex_data,
					 "fft of exp", complex_data,
					 n, 1e6);
    };
  msg_result_status (status);
  
  msg_checking_params (length, "gsl_fft_complex_wavetable_free");
  status = gsl_fft_complex_wavetable_free (&complex_wavetable);
  msg_result_status (status);
  
  /* check for memory leaks here if mstats available */
  
  free (complex_data);
  free (complex_tmp);
  free (fft_complex_data);
  free (fft_complex_tmp);
}
