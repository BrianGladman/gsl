#ifndef GSL_FFT_COMPLEX_H
#define GSL_FFT_COMPLEX_H

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
    double * scratch;
  }
gsl_fft_complex_wavetable;

  /*  Power of 2 routines  */


int
gsl_fft_complex_radix2_forward (double data[], size_t stride, size_t n);

int
gsl_fft_complex_radix2_backward (double data[], size_t stride, size_t n);

int
gsl_fft_complex_radix2_inverse (double data[], size_t stride, size_t n);

int
gsl_fft_complex_radix2 (double data[], size_t stride, size_t n,
			gsl_fft_direction sign);

int
gsl_fft_complex_radix2_dif_forward (double data[], size_t stride, size_t n);

int
gsl_fft_complex_radix2_dif_backward (double data[], size_t stride, size_t n);

int
gsl_fft_complex_radix2_dif_inverse (double data[], size_t stride, size_t n);

int
gsl_fft_complex_radix2_dif (double data[], size_t stride, size_t n,
			    gsl_fft_direction sign);

int gsl_fft_binary_logn (size_t n);

int gsl_fft_complex_bitreverse_order (double data[][], 
				      size_t stride,
				      size_t n,
				      size_t n_bits);

  /*  Mixed Radix general-N routines  */

int gsl_fft_complex_forward (double data[], size_t stride, size_t n,
			     const gsl_fft_complex_wavetable * wavetable);

int gsl_fft_complex_backward (double data[], size_t stride, size_t n,
			      const gsl_fft_complex_wavetable * wavetable);

int gsl_fft_complex_inverse (double data[], size_t stride, size_t n,
			     const gsl_fft_complex_wavetable * wavetable);

int
gsl_fft_complex (double data[], size_t stride, size_t n,
		 const gsl_fft_complex_wavetable * wavetable,
		 gsl_fft_direction sign);

int
gsl_fft_complex_init (size_t n,
		      gsl_fft_complex_wavetable * wavetable);

int
gsl_fft_complex_generate_wavetable (size_t n,
				    gsl_fft_complex_wavetable * wavetable);

gsl_fft_complex_wavetable * gsl_fft_complex_wavetable_alloc (size_t n);

void
  gsl_fft_complex_wavetable_free (gsl_fft_complex_wavetable * wavetable);

int
gsl_fft_complex_wavetable_cpy (gsl_fft_complex_wavetable * dest,
			       gsl_fft_complex_wavetable * src) ;

#endif /* GSL_FFT_COMPLEX_H */
