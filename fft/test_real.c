#include "bitreverse.h"
#include "factorize.h"
#include "signals.h"
#include "compare.h"

void FUNCTION(test_real,func) (size_t stride, size_t n);
void FUNCTION(test_real,bitreverse_order) (size_t stride, size_t n);
void FUNCTION(test_real,radix2) (size_t stride, size_t n);

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
  
  rw = FUNCTION(gsl_fft_real,alloc) (n);
  gsl_test (rw == 0, NAME(gsl_fft_real) 
	    "_alloc, n = %d, stride = %d", n, stride);
    
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

  hcw = FUNCTION(gsl_fft_halfcomplex,alloc) (n);
  gsl_test (hcw == 0, NAME(gsl_fft_halfcomplex) 
	    "_alloc, n = %d, stride = %d", n, stride);
  
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

  FUNCTION(gsl_fft_real,free) (rw);
  FUNCTION(gsl_fft_halfcomplex,free) (hcw);

  free(real_data) ;
  free(complex_data) ;
  free(complex_tmp) ;
  free(fft_complex_data) ;
}


void 
FUNCTION(test_real,bitreverse_order) (size_t stride, size_t n) 
{
  int status ;
  int result ;
  size_t logn, i ;

  BASE * tmp = (BASE *) malloc (n * stride * sizeof (BASE));
  BASE * data = (BASE *) malloc (n * stride * sizeof (BASE));
  BASE * reversed_data = (BASE *) malloc (n * stride * sizeof (BASE));
  
  for (i = 0; i <  stride * n; i++) 
    {
      data[i] = i ;
    }

  memcpy (tmp, data, n * stride * sizeof(BASE)) ;

  result = fft_binary_logn(n) ;
  
  if (result == -1) {
    abort() ;
  } else {
    logn = result ;
  }

  /* do a naive bit reversal as a baseline for testing the other routines */

  for (i = 0; i < n; i++) 
    {
      size_t i_tmp = i ;
      size_t j = 0 ;
      size_t bit ;

      for (bit = 0; bit < logn; bit++)
	{
	  j <<= 1;		/* reverse shift i into j */
	  j |= i_tmp & 1;
	  i_tmp >>= 1;
	}

      reversed_data[j*stride] = data[i*stride] ;
    }

  FUNCTION(fft_real,bitreverse_order) (data, stride, n, logn);

  status = FUNCTION(compare_real,results) ("naive bit reverse", 
					   reversed_data,
					   "gsl_fft_complex_bitreverse_order", 
					   data,
					   stride, n, 1e6);

  gsl_test (status, NAME(gsl_fft_real) "_bitreverse_order, n = %d", n);

  free (reversed_data) ;
  free (data) ;
  free (tmp) ;
}


void FUNCTION(test_real,radix2) (size_t stride, size_t n) 
{
  size_t i ;
  int status ;

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
  
  /* radix-2 real fft */
  
  FUNCTION(fft_signal,real_noise) (n, stride, complex_data, fft_complex_data);
  memcpy (complex_tmp, complex_data, 2 * n * stride * sizeof (BASE));

  for (i = 0; i < n; i++)
    {
      real_data[i*stride] = REAL(complex_data,stride,i);
    }
  
  FUNCTION(gsl_fft_real,radix2_transform) (real_data, stride, n);
  FUNCTION(gsl_fft_halfcomplex,radix2_unpack) (real_data, complex_data, stride, n);
  
  status = FUNCTION(compare_complex,results) ("dft", fft_complex_data,
					      "fft of noise", complex_data,
					      stride, n, 1e6);
  gsl_test (status, NAME(gsl_fft_real) 
	    "_radix2 with signal_real_noise, n = %d, stride = %d", n, stride);
  
  /* compute the inverse fft */
  
  status = FUNCTION(gsl_fft_halfcomplex,radix2_transform) (real_data, stride, n);
  
  for (i = 0; i < n; i++)
    {
      real_data[i*stride] /= n;
    }
  
  FUNCTION(gsl_fft_real,unpack) (real_data, complex_data, stride, n);
  
  status = FUNCTION(compare_complex,results) ("orig", complex_tmp,
					      "fft inverse", complex_data,
					      stride, n, 1e6);

  gsl_test (status, NAME(gsl_fft_halfcomplex) 
	    "_radix2 with data from signal_noise, n = %d, stride = %d", n, stride);


  free(real_data) ;
  free(complex_data) ;
  free(complex_tmp) ;
  free(fft_complex_data) ;
}
