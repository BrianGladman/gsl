#include "factorize.h"
#include "signals.h"
#include "compare.h"

void FUNCTION(test_complex,func) (size_t n, size_t stride);

void FUNCTION(test_complex,func) (size_t n, size_t stride) 
{
  size_t i ;
  int status ;

  TYPE(gsl_fft_wavetable_complex) * cw ;

  BASE * complex_data = (BASE *) malloc (2 * n * stride * sizeof (BASE));
  BASE * complex_tmp = (BASE *) malloc (2 * n * stride * sizeof (BASE));
  BASE * fft_complex_data = (BASE *) malloc (2 * n * stride * sizeof (BASE));
  BASE * fft_complex_tmp = (BASE *) malloc (2 * n * stride * sizeof (BASE));

  for (i = 0 ; i < 2 * n * stride ; i++)
    {
      complex_data[i] = i ;
      complex_tmp[i] = i ;
      fft_complex_data[i] = i ;
      fft_complex_tmp[i] = i ;
    }

  gsl_set_error_handler (NULL);	/* abort on any errors */

  /* Test allocation */

  {
    cw = FUNCTION(gsl_fft_complex,wavetable_alloc) (n);
    gsl_test (cw == 0, NAME(gsl_fft_complex) "_wavetable_alloc, n = %d", n);
  }

  /* Test initialization */

  {
    status = FUNCTION(gsl_fft_complex,init) (n, cw);
    gsl_test (status, NAME(gsl_fft_complex) "_init, n = %d", n);
  }

  /* Test wavetable generation */
  
  {
    status = FUNCTION(gsl_fft_complex,generate_wavetable) (n, cw);
    gsl_test (status, NAME(gsl_fft_complex) "_generate_wavetable, n = %d", n);
  }

  /* Test mixed radix fft with noise */

  {
    FUNCTION(fft_signal_complex,noise) (n, stride, complex_data, fft_complex_data);
    memcpy (complex_tmp, complex_data, 2 * n * stride * sizeof (BASE));
    
    FUNCTION(gsl_fft_complex,forward) (complex_data, stride, n, cw);
    memcpy (fft_complex_tmp, complex_data, 2 * n * stride * sizeof (BASE));
    
    status = FUNCTION(compare_complex,results) ("dft", fft_complex_data,
						"fft of noise", complex_data,
						n, stride, 1e6);
    gsl_test (status, NAME(gsl_fft_complex) "_forward with signal_noise, n = %d", n);
  }
  
  /* Test the inverse fft */

  {
    status = FUNCTION(gsl_fft_complex,inverse) (complex_data, stride, n, cw);
    status = FUNCTION(compare_complex,results) ("orig", complex_tmp,
						"fft inverse", complex_data,
						n, stride, 1e6);
    gsl_test (status, NAME(gsl_fft_complex) "_inverse with signal_noise, n = %d", n);
  }

  /* Test the backward fft */

  {
    status = FUNCTION(gsl_fft_complex,backward) (fft_complex_tmp, stride, n, cw);

    for (i = 0; i < n; i++)
      {
	REAL(complex_tmp,stride,i) *= n;
	IMAG(complex_tmp,stride,i) *= n;
      }
    status = FUNCTION(compare_complex,results) ("orig", 
						complex_tmp,
						"fft backward", 
						fft_complex_tmp,
						n, stride, 1e6);

    gsl_test (status, NAME(gsl_fft_complex) "_backward with signal_noise, n = %d", n);
  }

  /* Test a pulse signal */
  
  {
    FUNCTION(fft_signal_complex,pulse) (1, n, stride, 1.0, 0.0, complex_data,
					fft_complex_data);
    FUNCTION(gsl_fft_complex,forward) (complex_data, stride, n, cw);
    status = FUNCTION(compare_complex,results) ("analytic", fft_complex_data,
						"fft of pulse", complex_data, 
						n, stride, 1e6);
    gsl_test (status, NAME(gsl_fft_complex) "_forward with signal_pulse, n = %d", n);
  }


  /* Test a constant signal */

  {
    FUNCTION(fft_signal_complex,constant) (n, stride, 1.0, 0.0, complex_data,
					   fft_complex_data);
    FUNCTION(gsl_fft_complex,forward) (complex_data, stride, n, cw);
    status = FUNCTION(compare_complex,results) ("analytic", fft_complex_data,
						"fft of constant", 
						complex_data,
						n, stride, 1e6);
    gsl_test (status, 
	      NAME(gsl_fft_complex) "_forward with signal_constant, n = %d", n);
  }

  /* Test an exponential (cos/sin) signal */
  
  {
    status = 0;
    for (i = 0; i < n; i++)
      {
	FUNCTION(fft_signal_complex,exp) ((int)i, n, stride, 1.0, 0.0, complex_data,
					  fft_complex_data);
	FUNCTION(gsl_fft_complex,forward) (complex_data, stride, n, cw);
	status |= FUNCTION(compare_complex,results) ("analytic", 
						     fft_complex_data,
						     "fft of exp", 
						     complex_data,
						     n, stride, 1e6);
      }
    gsl_test (status, NAME(gsl_fft_complex) "_forward with signal_exp, n = %d", n);
  }

  FUNCTION(gsl_fft_complex,wavetable_free) (cw);
  gsl_test (status, NAME(gsl_fft_complex) "_wavetable_free, n = %d", n);
  
  free (complex_data);
  free (complex_tmp);
  free (fft_complex_data);
  free (fft_complex_tmp);
}

#ifdef JUNK
void test_real (size_t n) 
{
  size_t i ;
  int status ;

  double *real_data, *real_tmp;
  double *fft_real_data, *fft_real_tmp;
  double *complex_data, *complex_tmp;
  double *fft_complex_data, *fft_complex_tmp;

  gsl_fft_real_wavetable * rw ;
  gsl_fft_halfcomplex_wavetable * hcw ;

  real_data = (double *) malloc (n * sizeof (double));
  complex_data = (double *) malloc (n * 2 * sizeof (double));
  complex_tmp = (double *) malloc (n * 2 * sizeof (double));
  fft_complex_data = (double *) malloc (n * 2 * sizeof (double));
  fft_complex_tmp = (double *) malloc (n * 2 * sizeof (double));

  gsl_set_error_handler (NULL);	/* abort on any errors */
  
  /* mixed radix real fft */
  
  rw = gsl_fft_real_wavetable_alloc (n);
  gsl_test (rw == 0, "gsl_fft_real_wavetable_alloc, n = %d", n);
  
  status = gsl_fft_real_init (n, rw);
  gsl_test (status, "gsl_fft_real_wavetable_init, n = %d", n);
  
  status = gsl_fft_real_generate_wavetable (n, rw);
  gsl_test (status, "gsl_fft_real_generate_wavetable, n = %d", n);
  
  real_data = (double *) malloc (n * sizeof (double));
  real_tmp = (double *) malloc (n * sizeof (double));
  fft_real_data = (double *) malloc (n * sizeof (double));
  fft_real_tmp = (double *) malloc (n * sizeof (double));
  
  fft_signal_real_noise (n, complex_data, fft_complex_data);
  memcpy (complex_tmp, complex_data, n * sizeof (gsl_complex));

  for (i = 0; i < n; i++)
    {
      real_data[i] = REAL(complex_data,1,i);
    }
  
  gsl_fft_real (real_data, 1, n, rw);
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
  
  status = gsl_fft_halfcomplex (real_data, 1, n, hcw);
  
  for (i = 0; i < n; i++)
    {
      real_data[i] /= n;
    }
  
  gsl_fft_real_unpack (real_data, complex_data, n);
  
  status = compare_complex_results ("orig", complex_tmp,
				    "fft inverse", complex_data,
				    n, 1e6);
  gsl_test (status, "gsl_fft_halfcomplex with data from signal_noise, n = %d", n);
  
}
#endif
