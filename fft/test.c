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

void my_error_handler (const char *reason, const char *file, int line);

int
main (int argc, char *argv[])
{
  double *real_data, *real_tmp;
  double *fft_real_data, *fft_real_tmp;
  complex *complex_data, *complex_tmp;
  complex *fft_complex_data, *fft_complex_tmp;

  gsl_fft_complex_wavetable complex_wavetable;
  gsl_fft_real_wavetable real_wavetable;
  gsl_fft_halfcomplex_wavetable halfcomplex_wavetable;
  int i, status;
  unsigned int n, n_min, n_max;

  if (argc == 3)
    {
      n_min = strtol (argv[1], NULL, 0);
      n_max = strtol (argv[2], NULL, 0);
    }
  else if (argc == 2)
    {
      n_min = strtol (argv[1], NULL, 0);
      n_max = n_min;
    }
  else
    {
      printf ("test n_min n_max, test n\n");
      exit (0);
    }

  n = 1;

  real_data = malloc (n * sizeof (double));
  complex_data = malloc (n * sizeof (complex));
  complex_tmp = malloc (n * sizeof (complex));
  fft_complex_data = malloc (n * sizeof (complex));
  fft_complex_tmp = malloc (n * sizeof (complex));

  gsl_error_set_handler (&my_error_handler);


  /* n = 0 in alloc */

  msg_checking ("trap for n = 0 in gsl_fft_complex_wavetable_alloc ");
  status = gsl_fft_complex_wavetable_alloc (0, &complex_wavetable);
  msg_result_status (!status);

  msg_checking ("trap for n = 0 in gsl_fft_real_wavetable_alloc ");
  status = gsl_fft_real_wavetable_alloc (0, &real_wavetable);
  msg_result_status (!status);

  msg_checking ("trap for n = 0 in gsl_fft_halfcomplex_wavetable_alloc ");
  status = gsl_fft_halfcomplex_wavetable_alloc (0, &halfcomplex_wavetable);
  msg_result_status (!status);


  /* n = 0 in wavetable_init */

  msg_checking ("trap for n = 0 in gsl_fft_complex_wavetable_init");
  status = gsl_fft_complex_init (0, &complex_wavetable);
  msg_result_status (!status);

  msg_checking ("trap for n = 0 in gsl_fft_real_wavetable_init");
  status = gsl_fft_real_init (0, &real_wavetable);
  msg_result_status (!status);

  msg_checking ("trap for n = 0 in gsl_fft_halfcomplex_wavetable_init");
  status = gsl_fft_halfcomplex_init (0, &halfcomplex_wavetable);
  msg_result_status (!status);


  /* n = 0 in generate_wavetable */

  msg_checking ("trap for n = 0 in gsl_fft_complex_generate_wavetable");
  status = gsl_fft_complex_generate_wavetable (0, &complex_wavetable);
  msg_result_status (!status);

  msg_checking ("trap for n = 0 in gsl_fft_real_generate_wavetable");
  status = gsl_fft_real_generate_wavetable (0, &real_wavetable);
  msg_result_status (!status);

  msg_checking ("trap for n = 0 in gsl_fft_halfcomplex_generate_wavetable");
  status = gsl_fft_halfcomplex_generate_wavetable (0, &halfcomplex_wavetable);
  msg_result_status (!status);


  /* n = 0 in fft forward */

  msg_checking ("trap for n = 0 in gsl_fft_complex_forward");
  status = gsl_fft_complex_forward (complex_data, 0, &complex_wavetable);
  msg_result_status (!status);

  msg_checking ("trap for n = 0 in gsl_fft_real");
  status = gsl_fft_real (real_data, 0, &real_wavetable);
  msg_result_status (!status);

  msg_checking ("trap for n = 0 in gsl_fft_halfcomplex");
  status = gsl_fft_halfcomplex (real_data, 0, &halfcomplex_wavetable);
  msg_result_status (!status);

  msg_checking ("trap for n = 0 in gsl_fft_complex_radix2_forward");
  status = gsl_fft_complex_radix2_forward (complex_data, 0);
  msg_result_status (!status);


  /* n = 0 in fft backward */

  msg_checking ("trap for n = 0 in gsl_fft_complex_backward");
  status = gsl_fft_complex_backward (complex_data, 0, &complex_wavetable);
  msg_result_status (!status);

  msg_checking ("trap for n = 0 in gsl_fft_complex_radix2_backward");
  status = gsl_fft_complex_radix2_backward (complex_data, 0);
  msg_result_status (!status);

  /* n = 0 in fft inverse */

  msg_checking ("trap for n = 0 in gsl_fft_complex_inverse");
  status = gsl_fft_complex_inverse (complex_data, 0, &complex_wavetable);
  msg_result_status (!status);

  msg_checking ("trap for n = 0 in gsl_fft_complex_radix2_inverse");
  status = gsl_fft_complex_radix2_inverse (complex_data, 0);
  msg_result_status (!status);


  /* n != 2^k in power of 2 routines */

  msg_checking ("trap for n != 2^k in gsl_fft_complex_radix2_forward");
  status = gsl_fft_complex_radix2_forward (complex_data, 17);
  msg_result_status (!status);

  msg_checking ("trap for n != 2^k in gsl_fft_complex_radix2_backward");
  status = gsl_fft_complex_radix2_backward (complex_data, 17);
  msg_result_status (!status);

  msg_checking ("trap for n != 2^k in gsl_fft_complex_radix2_inverse");
  status = gsl_fft_complex_radix2_inverse (complex_data, 17);
  msg_result_status (!status);

  /* n != wavetable.n in mixed radix routines */

  msg_checking ("trap for n != nw in gsl_fft_complex_forward");
  complex_wavetable.n = 3;
  status = gsl_fft_complex_forward (complex_data, 4, &complex_wavetable);
  msg_result_status (!status);

  msg_checking ("trap for n != nw in gsl_fft_complex_backward");
  complex_wavetable.n = 3;
  status = gsl_fft_complex_backward (complex_data, 4, &complex_wavetable);
  msg_result_status (!status);

  msg_checking ("trap for n != nw in gsl_fft_complex_inverse");
  complex_wavetable.n = 3;
  status = gsl_fft_complex_inverse (complex_data, 4, &complex_wavetable);
  msg_result_status (!status);

  msg_checking ("trap for n != nw in gsl_fft_real");
  real_wavetable.n = 3;
  status = gsl_fft_real (real_data, 4, &real_wavetable);
  msg_result_status (!status);

  msg_checking ("trap for n != nw in gsl_fft_halfcomplex");
  halfcomplex_wavetable.n = 3;
  status = gsl_fft_halfcomplex (real_data, 4, &halfcomplex_wavetable);
  msg_result_status (!status);

  gsl_error_set_handler (NULL);	/* abort on any errors */

  for (n = n_min; n <= n_max; n++)
    {
      char length[256];
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
      msg_checking_params (length, "gsl_fft_complex_forward with test_signal_noise");
      gsl_fft_test_signal_complex_noise (n, complex_data, fft_complex_data);
      memcpy (complex_tmp, complex_data, n * sizeof (complex));
      gsl_fft_complex_forward (complex_data, n, &complex_wavetable);
      memcpy (fft_complex_tmp, complex_data, n * sizeof (complex));
      status = compare_complex_results ("dft", fft_complex_data,
					"fft of noise", complex_data,
					n, 10 + n * n / 100.0);
      msg_result_status (status);

      /* compute the inverse fft */
      msg_checking_params (length, "gsl_fft_complex_inverse with data from test_signal_noise");
      status = gsl_fft_complex_inverse (complex_data, n, &complex_wavetable);
      status = compare_complex_results ("orig", complex_tmp,
					"fft inverse", complex_data,
					n, 10 + n * n / 100.0);
      msg_result_status (status);

      /* compute the backward fft */
      msg_checking_params (length, "gsl_fft_complex_backward with data from test_signal_noise");
      status = gsl_fft_complex_backward (fft_complex_tmp, n, &complex_wavetable);
      for (i = 0; i < n; i++)
	{
	  complex_tmp[i].real *= n;
	  complex_tmp[i].imag *= n;
	}
      status = compare_complex_results ("orig", complex_tmp,
					"fft backward", fft_complex_tmp,
					n, 1000.0);
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

      msg_checking_params (length, "gsl_fft_real with test_signal_real_noise");
      gsl_fft_test_signal_real_noise (n, complex_data, fft_complex_data);
      memcpy (complex_tmp, complex_data, n * sizeof (complex));
      for (i = 0; i < n; i++)
	{
	  real_data[i] = complex_data[i].real;
	  fft_complex_data[i].imag *= -1;
	  /* fft_real is actually a backward tranform */
	}

      gsl_fft_real (real_data, n, &real_wavetable);
      gsl_fft_halfcomplex_unpack (real_data, complex_data, n);

      status = compare_complex_results ("dft", fft_complex_data,
					"fft of noise", complex_data,
					n, sqrt ((double) n) * 10.0);
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

      msg_checking_params (length, "gsl_fft_halfcomplex with data from test_signal_noise");
      status = gsl_fft_halfcomplex (real_data, n, &halfcomplex_wavetable);

      for (i = 0; i < n; i++)
	{
	  real_data[i] /= n;
	}

      gsl_fft_real_unpack (real_data, complex_data, n);

      status = compare_complex_results ("orig", complex_tmp,
					"fft inverse", complex_data,
					n, 1000.0);
      msg_result_status (status);

      /* pulse */
      msg_checking_params (length, "gsl_fft_complex_forward with test_signal_pulse");
      gsl_fft_test_signal_complex_pulse (1, n, 1.0, 0.0, complex_data,
					 fft_complex_data);
      gsl_fft_complex_forward (complex_data, n, &complex_wavetable);
      status = compare_complex_results ("analytic", fft_complex_data,
				   "fft of pulse", complex_data, n, 1000.0);
      msg_result_status (status);

      msg_checking_params (length, "gsl_fft_complex_forward with test_signal_constant");
      gsl_fft_test_signal_complex_constant (n, 1.0, 0.0, complex_data,
					    fft_complex_data);
      gsl_fft_complex_forward (complex_data, n, &complex_wavetable);
      status = compare_complex_results ("analytic", fft_complex_data,
					"fft of constant", complex_data,
					n, 1000.0);
      msg_result_status (status);

      msg_checking_params (length, "gsl_fft_complex_forward with test_signal_exp");
      status = 0;
      for (i = 0; i < n; i++)
	{
	  gsl_fft_test_signal_complex_exp (i, n, 1.0, 0.0, complex_data,
					   fft_complex_data);
	  gsl_fft_complex_forward (complex_data, n, &complex_wavetable);
	  status |= compare_complex_results ("analytic", fft_complex_data,
					     "fft of exp", complex_data,
					     n, 1000.0);
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

  for (n = 1; n <= n_max; n *= 2)
    {
      char length[256];
      sprintf (length, "n = %d", n);

      msg_checking_params (length, "gsl_fft_complex_radix2_dif with test_signal_noise");
      gsl_fft_test_signal_complex_noise (n, complex_data, fft_complex_data);
      memcpy (complex_tmp, complex_data, n * sizeof (complex));
      gsl_fft_complex_radix2_dif (complex_data, n, 1);
      status = compare_complex_results ("dft", fft_complex_data,
					"fft of noise", complex_data,
					n, 1000.0);
      msg_result_status (status);


      msg_checking_params (length, "gsl_fft_complex_radix2_forward with test_signal_noise");
      gsl_fft_test_signal_complex_noise (n, complex_data, fft_complex_data);
      memcpy (complex_tmp, complex_data, n * sizeof (complex));
      gsl_fft_complex_radix2_forward (complex_data, n);
      status = compare_complex_results ("dft", fft_complex_data,
					"fft of noise", complex_data,
					n, 1000.0);
      msg_result_status (status);

      /* compute the inverse fft */
      msg_checking_params (length, "gsl_fft_complex_radix2_inverse with data from test_signal_noise");
      status = gsl_fft_complex_radix2_inverse (complex_data, n);
      status = compare_complex_results ("orig", complex_tmp,
					"fft_real", complex_data,
					n, 1000.0);
      msg_result_status (status);


    }

  msg_summary (tests, passed, failed);

  return 0;
}


void
my_error_handler (const char *reason, const char *file, int line)
{
  printf ("(caught) ");
}
