#ifndef GSL_FFT_HALFCOMPLEX_FLOAT_H
#define GSL_FFT_HALFCOMPLEX_FLOAT_H

#include <stddef.h>

#include <gsl_math.h>
#include <gsl_complex.h>
#include <gsl_fft.h>

int gsl_fft_halfcomplex_float_radix2_backward (float data[], size_t stride, size_t n);
int gsl_fft_halfcomplex_float_radix2_inverse (float data[], size_t stride, size_t n);
int gsl_fft_halfcomplex_float_radix2_transform (float data[], size_t stride, size_t n);

typedef struct
  {
    size_t n;
    size_t nf;
    size_t factor[64];
    gsl_complex *twiddle[64];
    gsl_complex *trig;
    float *scratch;
  }
gsl_fft_wavetable_halfcomplex_float;

int gsl_fft_halfcomplex_float_backward (float data[], size_t stride, size_t n,
					const gsl_fft_wavetable_halfcomplex_float * wavetable);

int gsl_fft_halfcomplex_float_inverse (float data[], size_t stride, size_t n,
				       const gsl_fft_wavetable_halfcomplex_float * wavetable);

int gsl_fft_halfcomplex_float_transform (float data[], size_t stride, size_t n,
				   const gsl_fft_wavetable_halfcomplex_float * wavetable);

gsl_fft_wavetable_halfcomplex_float * gsl_fft_halfcomplex_float_alloc (size_t n);

void
gsl_fft_halfcomplex_float_free (gsl_fft_wavetable_halfcomplex_float * wavetable);

int
gsl_fft_halfcomplex_float_unpack (const float halfcomplex_coefficient[],
				  float complex_coefficient[],
				  size_t stride, size_t n);

int
gsl_fft_halfcomplex_float_radix2_unpack (const float halfcomplex_coefficient[],
					 float complex_coefficient[],
					 size_t stride, size_t n);

#endif /* GSL_FFT_HALFCOMPLEX_FLOAT_H */
