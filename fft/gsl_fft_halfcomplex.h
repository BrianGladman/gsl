#ifndef GSL_FFT_HALFCOMPLEX_H
#define GSL_FFT_HALFCOMPLEX_H

#include <stddef.h>

#include <gsl_math.h>
#include <gsl_complex.h>
#include <gsl_vector.h>
#include <gsl_fft.h>

int gsl_fft_halfcomplex_radix2_backward (double data[], size_t stride, size_t n);
int gsl_fft_halfcomplex_radix2_inverse (double data[], size_t stride, size_t n);
int gsl_fft_halfcomplex_radix2_transform (double data[], size_t stride, size_t n);

typedef struct
  {
    size_t n;
    size_t nf;
    size_t factor[64];
    gsl_complex *twiddle[64];
    gsl_complex *trig;
    double *scratch;
  }
gsl_fft_wavetable_halfcomplex;

int gsl_fft_halfcomplex_backward (double data[], size_t stride, size_t n,
				  const gsl_fft_wavetable_halfcomplex * wavetable);

int gsl_fft_halfcomplex_inverse (double data[], size_t stride, size_t n,
				 const gsl_fft_wavetable_halfcomplex * wavetable);

int gsl_fft_halfcomplex_transform (double data[], size_t stride, size_t n,
				   const gsl_fft_wavetable_halfcomplex * wavetable);

int gsl_fft_halfcomplex_init (size_t n,
			      gsl_fft_wavetable_halfcomplex * wavetable);

int gsl_fft_halfcomplex_generate_wavetable (size_t n,
					    gsl_fft_wavetable_halfcomplex * wavetable);

gsl_fft_wavetable_halfcomplex * gsl_fft_halfcomplex_wavetable_alloc (size_t n);

void
gsl_fft_halfcomplex_wavetable_free (gsl_fft_wavetable_halfcomplex * wavetable);

int
gsl_fft_halfcomplex_unpack (const double halfcomplex_coefficient[],
			    double complex_coefficient[],
			    size_t stride, size_t n);

#endif /* GSL_FFT_HALFCOMPLEX_H */
