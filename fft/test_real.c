#include "factorize.h"
#include "signals.h"
#include "compare.h"

void FUNCTION(test_real,func) (size_t stride, size_t n);

void FUNCTION(test_real,func) (size_t stride, size_t n) 
{
  size_t i ;
  int status ;

  TYPE(gsl_fft_wavetable_real) * rw ;
  TYPE(gsl_fft_wavetable_halfcomplex) * hcw ;

  BASE * real_data = (BASE *) malloc (n * stride * sizeof (BASE));
  BASE * complex_data = (BASE *) malloc (2 * n * stride * sizeof (BASE));
  BASE * complex_tmp = (BASE *) malloc (2 * n * stride * sizeof (BASE));
  BASE * fft_complex_data = (BASE *) malloc (2 * n * stride * sizeof (BASE));

  for (i = 0 ; i < n * stride ; i++)
    {
      real_data[i] = i ;
    }

  for (i = 0 ; i < 2 * n * stride ; i++)
    {
      complex_data[i] = i + 1000.0 ;
      complex_tmp[i] = i + 2000.0 ;
      fft_complex_data[i] = i + 3000.0 ;
    }

  gsl_set_error_handler (NULL);	/* abort on any errors */
  
  /* mixed radix real fft */
  
  rw = FUNCTION(gsl_fft_real,wavetable_alloc) (n);
  gsl_test (rw == 0, NAME(gsl_fft_real) 
	    "_wavetable_alloc, n = %d, stride = %d", n, stride);
  
  status = FUNCTION(gsl_fft_real,init) (n, rw);
  gsl_test (status, NAME(gsl_fft_real) 
	    "_wavetable_init, n = %d, stride = %d", n, stride);
  
  status = FUNCTION(gsl_fft_real,generate_wavetable) (n, rw);
  gsl_test (status, NAME(gsl_fft_real) 
	    "_generate_wavetable, n = %d, stride = %d", n, stride);
  
  FUNCTION(fft_signal,real_noise) (n, stride, complex_data, fft_complex_data);
  memcpy (complex_tmp, complex_data, 2 * n * stride * sizeof (BASE));

  for (i = 0; i < n; i++)
    {
      real_data[i*stride] = REAL(complex_data,stride,i);
    }
  
  FUNCTION(gsl_fft_real,transform) (real_data, stride, n, rw);
  FUNCTION(gsl_fft_halfcomplex,unpack) (real_data, complex_data, stride, n);
  
  status = FUNCTION(compare_complex,results) ("dft", fft_complex_data,
					      "fft of noise", complex_data,
					      stride, n, 1e6);
  gsl_test (status, NAME(gsl_fft_real) 
	    " with signal_real_noise, n = %d, stride = %d", n, stride);
  
  /* compute the inverse fft */

  hcw = FUNCTION(gsl_fft_halfcomplex,wavetable_alloc) (n);
  gsl_test (hcw == 0, NAME(gsl_fft_halfcomplex) 
	    "_wavetable_alloc, n = %d, stride = %d", n, stride);
  
  status = FUNCTION(gsl_fft_halfcomplex,init) (n, hcw);
  gsl_test (status, NAME(gsl_fft_halfcomplex) 
	    "_wavetable_init, n = %d, stride = %d", n, stride);
  
  status = FUNCTION(gsl_fft_halfcomplex,generate_wavetable) (n, hcw);
  gsl_test (status, NAME(gsl_fft_halfcomplex) 
	    "_generate_wavetable, n = %d, stride = %d", n, stride);
  
  status = FUNCTION(gsl_fft_halfcomplex,transform) (real_data, stride, n, hcw);
  
  for (i = 0; i < n; i++)
    {
      real_data[i*stride] /= n;
    }
  
  FUNCTION(gsl_fft_real,unpack) (real_data, complex_data, stride, n);
  
  status = FUNCTION(compare_complex,results) ("orig", complex_tmp,
					      "fft inverse", complex_data,
					      stride, n, 1e6);

  gsl_test (status, NAME(gsl_fft_halfcomplex) 
	    " with data from signal_noise, n = %d, stride = %d", n, stride);
  

  free(real_data) ;
  free(complex_data) ;
  free(complex_tmp) ;
  free(fft_complex_data) ;
}



