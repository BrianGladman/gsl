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

#include <getopt.h>
#include <compare.h>

void usage (void);
void check_complex_bitreverse_order (unsigned int n) ;
void check_complex_radix2 (unsigned int n) ;
void check_real_radix2 (unsigned int n) ;

void
usage (void)
{
  printf("Usage: test_radix2 [OPTION]\n"
"Exercise the radix-2 fft routines for length n. By default n runs\n"
"through small powers of two: 1, 2, 4, 8, ... , 1024.\n"
"\n"
"  -n, --number=NUM       tests on length n\n"
"  -v, --verbose          verbosely list tests\n"
"\n"
"Without the -v option the test is quiet. The exit status indicates\n"
"success or failure.\n"
) ; 
  exit(0) ;
}


int
main (int argc, char *argv[])
{
  unsigned int i ;
  unsigned int n = 0;

  while (1) {

    static struct option long_options[] = 
    {
      {"verbose", 0, 0, 'v'},
      {"number", 1, 0, 'n'},
      {"help", 0, 0, 'h'},
      {0, 0, 0, 0}
    } ;

    int option_index = 0 ;

    int c = getopt_long (argc, argv, "hn:v",
			 long_options, &option_index) ;
   
    if (c == -1)   /* end of options */
      break ;   

    if (c == 0 && long_options[option_index].flag == 0)
      c = long_options[option_index].val;

    switch (c) 
      {
      case 'v':
	/* gsl_test_verbose () ; */
	break ;
      case 'n':
	if (optarg) 
	  n = strtol (optarg, NULL, 0);
	else 
	  usage () ;
	break ;
      case 'h':
      default:
	usage () ;
      }
  }

  if (n)
    {
      check_complex_bitreverse_order (n) ;
      check_complex_radix2 (n) ;
      check_real_radix2 (n) ;
    }
  else
    {
      for (i = 1 ; i <= 1024 ; i *= 2) 
	{
	  check_complex_bitreverse_order (i) ;
	  check_complex_radix2 (i) ;
	  check_real_radix2 (i) ;
	}
    }

  return gsl_test_summary ();

}


void check_complex_bitreverse_order (unsigned int n) 
{
  int status ;
  int result ;
  unsigned int logn, i ;
  complex *complex_data, *complex_tmp, *complex_reversed_data;

  complex_tmp = malloc (n * sizeof (complex));
  complex_data = malloc (n * sizeof (complex));
  complex_reversed_data = malloc (n * sizeof (complex));
  
  for (i = 0; i < n; i++) 
    {
      complex_data[i].real = i + 1;
      complex_data[i].real = n + i + 1 ;
    }

  memcpy (complex_tmp, complex_data, n * sizeof(complex)) ;

  result = gsl_fft_binary_logn(n) ;
  
  if (result == -1) {
    abort() ;
  } else {
    logn = result ;
  }

  /* do a naive bit reversal as a baseline for testing the other routines */

  for (i = 0; i < n; i++) 
    {
      
      unsigned int i_tmp = i ;
      unsigned int j = 0 ;
      unsigned int bit ;

      for (bit = 0; bit < logn; bit++)
	{
	  j <<= 1;		/* reverse shift i into j */
	  j |= i_tmp & 1;
	  i_tmp >>= 1;
	}

      complex_reversed_data[j] = complex_data[i] ;
    }

  gsl_fft_complex_bitreverse_order (complex_data, n, logn);
  status = compare_complex_results ("naive bit reverse", 
				    complex_reversed_data,
				    "gsl_fft_complex_bitreverse_order", 
				    complex_data,
				    n, 1e6);

  gsl_test (status, "gsl_fft_complex_bitreverse_order, n = %d", n);

#ifdef UNUSED
  memcpy (complex_data, complex_tmp, n * sizeof(complex)) ;
  gsl_fft_complex_goldrader_bitreverse_order (complex_data, n);
  status = compare_complex_results ("naive bit reverse", 
				    complex_reversed_data,
				    "gsl_fft_complex_goldrader_bitreverse_order", 
				    complex_data,
				    n, 1e6);
  gsl_test (status, "gsl_fft_complex_goldrader_bitreverse_order, n = %d", n);

  memcpy (complex_data, complex_tmp, n * sizeof(complex)) ;
  gsl_fft_complex_rodriguez_bitreverse_order (complex_data, n, logn);
  status = compare_complex_results ("naive bit reverse", 
				    complex_reversed_data,
				    "gsl_fft_complex_rodriguez_bit_reverse", 
				    complex_data, 
				    n, 1e6);
  gsl_test (status, "gsl_fft_complex_rodriguez_bitreverse_order, n = %d", n);
#endif

  free (complex_data) ;
  free (complex_tmp) ;
}

void check_complex_radix2 (unsigned int n) 
{
  int status ;

  complex *complex_data, *complex_tmp;
  complex *fft_complex_data, *fft_complex_tmp;

  complex_data = malloc (n * sizeof (complex));
  complex_tmp = malloc (n * sizeof (complex));
  fft_complex_data = malloc (n * sizeof (complex));
  fft_complex_tmp = malloc (n * sizeof (complex));
  
  gsl_fft_signal_complex_noise (n, complex_data, fft_complex_data);
  memcpy (complex_tmp, complex_data, n * sizeof (complex));
  gsl_fft_complex_radix2_dif_forward (complex_data, n);
  status = compare_complex_results ("dft", fft_complex_data,
				    "fft of noise", complex_data,
				    n, 1e6);
  gsl_test (status, "gsl_fft_complex_radix2_dif with signal_noise, n = %d", n);


  gsl_fft_signal_complex_noise (n, complex_data, fft_complex_data);
  memcpy (complex_tmp, complex_data, n * sizeof (complex));
  gsl_fft_complex_radix2_forward (complex_data, n);
  status = compare_complex_results ("dft", fft_complex_data,
				    "fft of noise", complex_data,
				    n, 1e6);
  gsl_test (status, "gsl_fft_complex_radix2_forward with signal_noise, n = %d", n);
  
  /* compute the inverse fft */
  status = gsl_fft_complex_radix2_inverse (complex_data, n);
  status = compare_complex_results ("orig", complex_tmp,
				    "fft_real", complex_data,
				    n, 1e6);
  gsl_test (status, "gsl_fft_complex_radix2_inverse with signal_noise, n = %d", n);

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
  
  gsl_fft_signal_real_noise (n, complex_data, fft_complex_data);
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
  gsl_test (status, "gsl_fft_real_radix2 with signal_noise, n = %d", n);

  status = gsl_fft_halfcomplex_radix2 (real_data, n) ;

   for (i = 0; i < n ; i++) { 
     real_data[i] /= n ; 
   } 

  gsl_fft_real_unpack (real_data, complex_data, n) ;

  status = compare_complex_results ("orig", complex_tmp,
				    "fft inverse", complex_data,
				    n, 1e6);
  gsl_test (status, "gsl_fft_halfcomplex_radix2 with signal_noise, n = %d", n);

  free (complex_data) ;
  free (complex_tmp) ;
  free (fft_complex_data) ;
  free (fft_complex_tmp) ;
  
}
