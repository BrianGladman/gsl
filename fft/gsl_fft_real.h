#ifndef GSL_FFT_REAL_H
#define GSL_FFT_REAL_H

#include <stddef.h>

#include <gsl_math.h>
#include <gsl_complex.h>
#include <gsl_fft.h>

int gsl_fft_real_radix2 (double data[], size_t stride, size_t n) ;

typedef struct
  {
    size_t n;
    size_t nf;
    size_t factor[64];
    gsl_complex *twiddle[64];
    gsl_complex *trig;
    double *scratch;
  }
gsl_fft_real_wavetable;

int gsl_fft_real (double data[], size_t stride, size_t n,
		  const gsl_fft_real_wavetable * wavetable);

int gsl_fft_real_init (size_t n, gsl_fft_real_wavetable * wavetable);

int gsl_fft_real_generate_wavetable (size_t n,
				     gsl_fft_real_wavetable * wavetable);

gsl_fft_real_wavetable * gsl_fft_real_wavetable_alloc (size_t n);

void  gsl_fft_real_wavetable_free (gsl_fft_real_wavetable * wavetable);

int gsl_fft_real_unpack (const double real_coefficient[],
			 double complex_coefficient[],
			 size_t n);

#endif /* GSL_FFT_REAL_H */
