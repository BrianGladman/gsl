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
void check_complex (unsigned int n) ;

void
usage (void)
{
  printf("Usage: test [OPTION]\n"
"Exercise the fft routines for length n. By default n runs from 1 to 100.\n"
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
  unsigned int i;
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

void check_complex (unsigned int n) 
{
  unsigned int i ;
  int status ;

  double *real_data, *real_tmp;
  double *fft_real_data, *fft_real_tmp;
  complex *complex_data, *complex_tmp;
  complex *fft_complex_data, *fft_complex_tmp;

  gsl_fft_complex_wavetable complex_wavetable;
  gsl_fft_real_wavetable real_wavetable;
  gsl_fft_halfcomplex_wavetable halfcomplex_wavetable;

  real_data = malloc (n * sizeof (double));
  complex_data = malloc (n * sizeof (complex));
  complex_tmp = malloc (n * sizeof (complex));
  fft_complex_data = malloc (n * sizeof (complex));
  fft_complex_tmp = malloc (n * sizeof (complex));

  gsl_set_error_handler (NULL);	/* abort on any errors */

  status = gsl_fft_complex_wavetable_alloc (n, &complex_wavetable);
  gsl_test (status, "gsl_fft_complex_wavetable_alloc, n = %d", n);
  
  status = gsl_fft_complex_init (n, &complex_wavetable);
  gsl_test (status, "gsl_fft_complex_wavetable_init, n = %d", n);
  
  /* the wavetable generation is also performed by gsl_fft_complex_init
     but test it here too */
  
  status = gsl_fft_complex_generate_wavetable (n, &complex_wavetable);
  gsl_test (status, "gsl_fft_complex_generate_wavetable, n = %d", n);
  
  complex_data = malloc (n * sizeof (complex));
  complex_tmp = malloc (n * sizeof (complex));
  fft_complex_data = malloc (n * sizeof (complex));
  fft_complex_tmp = malloc (n * sizeof (complex));
  
  real_data = malloc (n * sizeof (double));
  
  /* mixed radix fft */
  gsl_fft_signal_complex_noise (n, complex_data, fft_complex_data);
  memcpy (complex_tmp, complex_data, n * sizeof (complex));
  gsl_fft_complex_forward (complex_data, n, &complex_wavetable);
  memcpy (fft_complex_tmp, complex_data, n * sizeof (complex));
  status = compare_complex_results ("dft", fft_complex_data,
				    "fft of noise", complex_data,
				    n, 1e6);
  gsl_test (status, "gsl_fft_complex_forward with signal_noise, n = %d", n);
  
  /* compute the inverse fft */
  status = gsl_fft_complex_inverse (complex_data, n, &complex_wavetable);
  status = compare_complex_results ("orig", complex_tmp,
				    "fft inverse", complex_data,
				    n, 1e6);
  gsl_test (status, "gsl_fft_complex_inverse with signal_noise, n = %d", n);
  
  /* compute the backward fft */

  status = gsl_fft_complex_backward (fft_complex_tmp, n, &complex_wavetable);

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
  
  status = gsl_fft_real_wavetable_alloc (n, &real_wavetable);
  gsl_test (status, "gsl_fft_real_wavetable_alloc, n = %d", n);
  
  status = gsl_fft_real_init (n, &real_wavetable);
  gsl_test (status, "gsl_fft_real_wavetable_init, n = %d", n);
  
  status = gsl_fft_real_generate_wavetable (n, &real_wavetable);
  gsl_test (status, "gsl_fft_real_generate_wavetable, n = %d", n);
  
  real_data = malloc (n * sizeof (double));
  real_tmp = malloc (n * sizeof (double));
  fft_real_data = malloc (n * sizeof (double));
  fft_real_tmp = malloc (n * sizeof (double));
  
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
  gsl_test (status, "gsl_fft_real with signal_real_noise, n = %d", n);
  
  /* compute the inverse fft */
  

  status = gsl_fft_halfcomplex_wavetable_alloc (n, &halfcomplex_wavetable);
  gsl_test (status, "gsl_fft_halfcomplex_wavetable_alloc, n = %d", n);
  
  status = gsl_fft_halfcomplex_init (n, &halfcomplex_wavetable);
  gsl_test (status, "gsl_fft_halfcomplex_wavetable_init, n = %d", n);
  
  status = gsl_fft_halfcomplex_generate_wavetable (n,
						   &halfcomplex_wavetable);
  gsl_test (status, "gsl_fft_halfcomplex_generate_wavetable, n = %d", n);
  
  status = gsl_fft_halfcomplex (real_data, n, &halfcomplex_wavetable);
  
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
  gsl_fft_complex_forward (complex_data, n, &complex_wavetable);
  status = compare_complex_results ("analytic", fft_complex_data,
				    "fft of pulse", complex_data, n, 1e6);
  gsl_test (status, "gsl_fft_complex_forward with signal_pulse, n = %d", n);
  
  gsl_fft_signal_complex_constant (n, 1.0, 0.0, complex_data,
					fft_complex_data);
  gsl_fft_complex_forward (complex_data, n, &complex_wavetable);
  status = compare_complex_results ("analytic", fft_complex_data,
				    "fft of constant", complex_data,
				    n, 1e6);
  gsl_test (status, "gsl_fft_complex_forward with signal_constant, n = %d", n);
  
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
  gsl_test (status, "gsl_fft_complex_forward with signal_exp, n = %d", n);
  
  status = gsl_fft_complex_wavetable_free (&complex_wavetable);
  gsl_test (status, "gsl_fft_complex_wavetable_free, n = %d", n);
  
  /* check for memory leaks here if mstats available */
  
  free (complex_data);
  free (complex_tmp);
  free (fft_complex_data);
  free (fft_complex_tmp);
}
