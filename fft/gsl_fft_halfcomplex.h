#ifndef _GSL_FFT_HALFCOMPLEX_H
#define _GSL_FFT_HALFCOMPLEX_H

#include <gsl_math.h>
#include <gsl_complex.h>
#include <gsl_fft.h>

int
gsl_fft_halfcomplex_radix2_backward (double data[],
				      const unsigned int n) ;

int
gsl_fft_halfcomplex_radix2_inverse (double data[],
				    const unsigned int n) ;

int
gsl_fft_halfcomplex_radix2 (double data[],
			    const unsigned int n) ;

typedef struct
  {
    unsigned int n;
    unsigned int nf;
    unsigned int factor[64];
    gsl_complex *twiddle[64];
    gsl_complex *trig;
    gsl_complex *scratch;
  }
gsl_fft_halfcomplex_wavetable;

int
  gsl_fft_halfcomplex_backward (double *data,
				 const unsigned int n,
			   const gsl_fft_halfcomplex_wavetable * wavetable);

int
  gsl_fft_halfcomplex_inverse (double *data,
			       const unsigned int n,
			   const gsl_fft_halfcomplex_wavetable * wavetable);

int
  gsl_fft_halfcomplex (double *data,
		       const unsigned int n,
		       const gsl_fft_halfcomplex_wavetable * wavetable);

int
  gsl_fft_halfcomplex_init (unsigned int n,
			    gsl_fft_halfcomplex_wavetable * wavetable);

int
  gsl_fft_halfcomplex_generate_wavetable (unsigned int n,
				 gsl_fft_halfcomplex_wavetable * wavetable);

int
  gsl_fft_halfcomplex_wavetable_alloc (unsigned int n,
				 gsl_fft_halfcomplex_wavetable * wavetable);

int
  gsl_fft_halfcomplex_wavetable_free (gsl_fft_halfcomplex_wavetable * wavetable);

int
  gsl_fft_halfcomplex_unpack (const double halfcomplex_coefficient[],
			      gsl_complex complex_coefficient[],
			      const unsigned int n);

#endif /* _GSL_FFT_HALFCOMPLEX_H */
