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
#include <gsl_fft_test.h>

#include <autotest.h>
#include <compare.h>

int verbose = 0;

unsigned int tests = 0;
unsigned int passed = 0;
unsigned int failed = 0;

void check_complex_radix2 (unsigned int n) ;
void check_real_radix2 (unsigned int n) ;

int
main (int argc, char *argv[])
{
  unsigned int i, n;
  int status ;

  if (argc == 2)
    {
      n = strtol (argv[1], NULL, 0);
      check_complex_radix2 (n) ;
      check_real_radix2 (n) ;
    }
  else
    {
      for (i = 1 ; i <= 1024 ; i *= 2) 
	{
	  check_complex_radix2 (i) ;
	  check_real_radix2 (i) ;
	}
    }

  status = msg_summary (tests, passed, failed);

  return status ;
}


void check_complex_radix2 (unsigned int n) 
{
  int status ;

  complex *complex_data, *complex_tmp;
  complex *fft_complex_data, *fft_complex_tmp;

  char length[256];

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
				    n, 1e6);
  msg_result_status (status);
  
  
  msg_checking_params (length, 
		       "gsl_fft_complex_radix2_forward with test_signal_noise");
  gsl_fft_test_signal_complex_noise (n, complex_data, fft_complex_data);
  memcpy (complex_tmp, complex_data, n * sizeof (complex));
  gsl_fft_complex_radix2_forward (complex_data, n);
  status = compare_complex_results ("dft", fft_complex_data,
				    "fft of noise", complex_data,
				    n, 1e6);
  msg_result_status (status);
  
  /* compute the inverse fft */
  msg_checking_params (length, 
		       "gsl_fft_complex_radix2_inverse with data from test_signal_noise");
  status = gsl_fft_complex_radix2_inverse (complex_data, n);
  status = compare_complex_results ("orig", complex_tmp,
				    "fft_real", complex_data,
				    n, 1e6);
  msg_result_status (status);

  free (complex_data) ;
  free (complex_tmp) ;
  free (fft_complex_data) ;
  free (fft_complex_tmp) ;
  
}



void check_real_radix2 (unsigned int n) 
{

  unsigned int i ;
  int status ;

  double *real_data, *real_tmp;
  double *fft_real_data, *fft_real_tmp;

  complex *complex_data, *complex_tmp;
  complex *fft_complex_data, *fft_complex_tmp;

  char length[256];

  real_data = malloc (n * sizeof (double));
  real_tmp = malloc (n * sizeof (double));

  fft_real_data = malloc (n * sizeof (double));
  fft_real_tmp = malloc (n * sizeof (double));

  complex_data = malloc (n * sizeof (complex));
  complex_tmp = malloc (n * sizeof (complex));

  fft_complex_data = malloc (n * sizeof (complex));
  fft_complex_tmp = malloc (n * sizeof (complex));
  
  sprintf (length, "n = %d", n);
  
  msg_checking_params (length, "gsl_fft_real_radix2 with test_signal_noise");

  gsl_fft_test_signal_real_noise (n, complex_data, fft_complex_data);
  memcpy (complex_tmp, complex_data, n * sizeof (complex));

  for (i = 0; i < n; i++)
    {
      real_data[i] = complex_data[i].real;
    }

  gsl_fft_real_radix2 (real_data, n);

  complex_data[0].real = real_data[0] ;
  complex_data[0].imag = 0.0 ;
  for (i = 1 ; i < n/2 ; i++) {
    complex_data[i].real = real_data[i] ;
    complex_data[i].imag = real_data[n-i] ;
    complex_data[n-i].real = real_data[i] ;
    complex_data[n-i].imag = -real_data[n-i] ;
  }
  complex_data[n/2].real = real_data[n/2] ;
  complex_data[n/2].imag = 0.0 ;

  status = compare_complex_results ("dft", fft_complex_data,
				    "fft of noise", complex_data,
				    n, 1e6);
  msg_result_status (status);

  msg_checking_params (length, "gsl_fft_halfcomplex_radix2 with data from test_signal_noise");
  status = gsl_fft_halfcomplex_radix2 (real_data, n) ;

   for (i = 0; i < n ; i++) { 
     real_data[i] /= n ; 
   } 

  gsl_fft_real_unpack (real_data, complex_data, n) ;

  status = compare_complex_results ("orig", complex_tmp,
				    "fft inverse", complex_data,
				    n, 1e6);
  msg_result_status (status);

  free (complex_data) ;
  free (complex_tmp) ;
  free (fft_complex_data) ;
  free (fft_complex_tmp) ;
  
}






