#ifndef __GSL_FFT_REAL_FLOAT_H__
#define __GSL_FFT_REAL_FLOAT_H__

#include <stddef.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_fft.h>

int gsl_fft_real_float_radix2_transform (float data[], size_t stride, size_t n) ;

typedef struct
  {
    size_t n;
    size_t nf;
    size_t factor[64];
    gsl_complex *twiddle[64];
    gsl_complex *trig;
    float *scratch;
  }
gsl_fft_wavetable_real_float;

int gsl_fft_real_float_transform (float data[], size_t stride, size_t n,
			    const gsl_fft_wavetable_real_float * wavetable);

gsl_fft_wavetable_real_float * gsl_fft_real_float_alloc (size_t n);

void  gsl_fft_real_float_free (gsl_fft_wavetable_real_float * wavetable);

int gsl_fft_real_float_unpack (const float real_float_coefficient[],
			       float complex_coefficient[],
			       size_t stride, size_t n);

#endif /* __GSL_FFT_REAL_FLOAT_H__ */
