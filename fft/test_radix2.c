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
#include <gsl_fft_test.h>

#include <autotest.h>
#include <compare.h>

int verbose = 0;

unsigned int tests = 0;
unsigned int passed = 0;
unsigned int failed = 0;

int
main (int argc, char *argv[])
{

  double *real_data, *real_tmp;
  double *fft_real_data, *fft_real_tmp;
  complex *complex_data, *complex_tmp;
  complex *fft_complex_data, *fft_complex_tmp;
  char length[256];

  gsl_fft_complex_wavetable complex_wavetable;
  gsl_fft_real_wavetable real_wavetable;
  gsl_fft_halfcomplex_wavetable halfcomplex_wavetable;
  int i, status;
  unsigned int n, n_min, n_max;

  if (argc == 2)
    {
      n = strtol (argv[1], NULL, 0);
    }
  else
    {
      printf ("test n\n");
      exit (EXIT_FAILURE);
    }
  
  real_data = malloc (n * sizeof (double));
  complex_data = malloc (n * sizeof (complex));
  complex_tmp = malloc (n * sizeof (complex));
  fft_complex_data = malloc (n * sizeof (complex));
  fft_complex_tmp = malloc (n * sizeof (complex));

  sprintf (length, "n = %d", n);
  
  msg_checking_params (length, 
		       "gsl_fft_complex_radix2_dif with test_signal_noise");
  gsl_fft_test_signal_complex_noise (n, complex_data, fft_complex_data);
  memcpy (complex_tmp, complex_data, n * sizeof (complex));
  gsl_fft_complex_radix2_dif_forward (complex_data, n);
  status = compare_complex_results ("dft", fft_complex_data,
				    "fft of noise", complex_data,
				    n, 1000.0);
  msg_result_status (status);
  
  
  msg_checking_params (length, 
		       "gsl_fft_complex_radix2_forward with test_signal_noise");
  gsl_fft_test_signal_complex_noise (n, complex_data, fft_complex_data);
  memcpy (complex_tmp, complex_data, n * sizeof (complex));
  gsl_fft_complex_radix2_forward (complex_data, n);
  status = compare_complex_results ("dft", fft_complex_data,
				    "fft of noise", complex_data,
				    n, 1000.0);
  msg_result_status (status);
  
  /* compute the inverse fft */
  msg_checking_params (length, 
		       "gsl_fft_complex_radix2_inverse with data from test_signal_noise");
  status = gsl_fft_complex_radix2_inverse (complex_data, n);
  status = compare_complex_results ("orig", complex_tmp,
				    "fft_real", complex_data,
				    n, 1000.0);
  msg_result_status (status);
  
  
}




