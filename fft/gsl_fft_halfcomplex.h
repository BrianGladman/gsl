#ifndef GSL_FFT_HALFCOMPLEX_H
#define GSL_FFT_HALFCOMPLEX_H

#include <stddef.h>

#include <gsl_math.h>
#include <gsl_complex.h>
#include <gsl_vector.h>
#include <gsl_fft.h>

int gsl_fft_halfcomplex_radix2_backward (double data[], size_t stride, size_t n);
int gsl_fft_halfcomplex_radix2_inverse (gsl_vector * data);
int gsl_fft_halfcomplex_radix2 (gsl_vector * data);

typedef struct
  {
    size_t n;
    size_t nf;
    size_t factor[64];
    gsl_complex *twiddle[64];
    gsl_complex *trig;
    double *scratch;
  }
gsl_fft_halfcomplex_wavetable;

int
gsl_fft_halfcomplex_backward (gsl_vector * data,
			      const gsl_fft_halfcomplex_wavetable * wavetable);

int
gsl_fft_halfcomplex_inverse (gsl_vector * data,
			     const gsl_fft_halfcomplex_wavetable * wavetable);

int
gsl_fft_halfcomplex (gsl_vector * data,
		     const gsl_fft_halfcomplex_wavetable * wavetable);

int
gsl_fft_halfcomplex_init (size_t n,
			  gsl_fft_halfcomplex_wavetable * wavetable);

int
gsl_fft_halfcomplex_generate_wavetable (size_t n,
					gsl_fft_halfcomplex_wavetable * wavetable);

gsl_fft_halfcomplex_wavetable * gsl_fft_halfcomplex_wavetable_alloc (size_t n);

void
gsl_fft_halfcomplex_wavetable_free (gsl_fft_halfcomplex_wavetable * wavetable);

int
gsl_fft_halfcomplex_unpack (const double halfcomplex_coefficient[],
			    gsl_complex complex_coefficient[],
			    size_t n);

#endif /* GSL_FFT_HALFCOMPLEX_H */
