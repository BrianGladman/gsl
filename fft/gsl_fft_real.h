#ifndef _GSL_FFT_REAL_H
#define _GSL_FFT_REAL_H

#include <gsl_math.h>
#include <gsl_complex.h>
#include <gsl_fft.h>

int
gsl_fft_real_radix2 (double data[],
		     const unsigned int n) ;

typedef struct
  {
    unsigned int n;
    unsigned int nf;
    unsigned int factor[64];
    complex *twiddle[64];
    complex *trig;
    double *scratch;
  }
gsl_fft_real_wavetable;

int
  gsl_fft_real (double data[],
		const unsigned int n,
		const gsl_fft_real_wavetable * wavetable);

int
  gsl_fft_real_init (unsigned int n,
		     gsl_fft_real_wavetable * wavetable);

int
  gsl_fft_real_generate_wavetable (unsigned int n,
				   gsl_fft_real_wavetable * wavetable);

int
  gsl_fft_real_wavetable_alloc (unsigned int n,
				gsl_fft_real_wavetable * wavetable);

int
  gsl_fft_real_wavetable_free (gsl_fft_real_wavetable * wavetable);

int
  gsl_fft_real_unpack (const double real_coefficient[],
		       complex complex_coefficient[],
		       const unsigned int n);

#endif /* _GSL_FFT_REAL_H */
