#ifndef GSL_FFT_COMPLEX_FLOAT_H
#define GSL_FFT_COMPLEX_FLOAT_H

#include <stddef.h>

#include <gsl_math.h>
#include <gsl_complex.h>
#include <gsl_fft.h>

typedef struct
  {
    size_t n;
    size_t nf;
    size_t factor[64];
    gsl_complex * twiddle[64];
    gsl_complex * trig;
    float * scratch;
  }
gsl_fft_wavetable_complex_float;

/*  Power of 2 routines  */

int
gsl_fft_complex_float_radix2_forward (float data[], size_t stride, size_t n);

int
gsl_fft_complex_float_radix2_backward (float data[], size_t stride, size_t n);

int
gsl_fft_complex_float_radix2_inverse (float data[], size_t stride, size_t n);

int
gsl_fft_complex_float_radix2 (float data[], 
			      size_t stride, 
			      size_t n,
			      gsl_fft_direction sign);

int
gsl_fft_complex_float_radix2_dif_forward (float data[], 
					  size_t stride, 
					  size_t n);

int
gsl_fft_complex_float_radix2_dif_backward (float data[], 
					   size_t stride, 
					   size_t n);

int
gsl_fft_complex_float_radix2_dif_inverse (float data[], 
					  size_t stride, 
					  size_t n);

int
gsl_fft_complex_float_radix2_dif (float data[], 
				  size_t stride, 
				  size_t n,
				  gsl_fft_direction sign);

int gsl_fft_binary_logn (size_t n);

int gsl_fft_complex_float_bitreverse_order (float data[][], 
					    size_t stride,
					    size_t n,
					    size_t n_bits);

/*  Mixed Radix general-N routines  */

int 
gsl_fft_complex_float_forward (float data[], 
			       size_t stride, 
			       size_t n,
			       const gsl_fft_wavetable_complex_float * wavetable);

int 
gsl_fft_complex_float_backward (float data[], 
				size_t stride, 
				size_t n,
				const gsl_fft_wavetable_complex_float * wavetable);

int 
gsl_fft_complex_float_inverse (float data[], 
			       size_t stride, 
			       size_t n,
			       const gsl_fft_wavetable_complex_float * wavetable);

int
gsl_fft_complex_float_transform (float data[], 
		       size_t stride, 
		       size_t n,
		       const gsl_fft_wavetable_complex_float * wavetable,
		       gsl_fft_direction sign);

int
gsl_fft_complex_float_init (size_t n,
			    gsl_fft_wavetable_complex_float * wavetable);

int
gsl_fft_complex_float_generate_wavetable (size_t n,
					  gsl_fft_wavetable_complex_float * wavetable);

gsl_fft_wavetable_complex_float * gsl_fft_complex_float_wavetable_alloc (size_t n);

void
gsl_fft_complex_float_wavetable_free (gsl_fft_wavetable_complex_float * wavetable);

int
gsl_fft_complex_float_wavetable_cpy (gsl_fft_wavetable_complex_float * dest,
				     gsl_fft_wavetable_complex_float * src) ;

#endif /* GSL_FFT_COMPLEX_FLOAT_H */
