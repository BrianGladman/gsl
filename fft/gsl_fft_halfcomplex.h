#ifndef _GSL_FFT_HALFCOMPLEX_H
#define _GSL_FFT_HALFCOMPLEX_H

#include <gsl_math.h>
#include <gsl_complex.h>
#include <gsl_fft.h>

typedef struct
  {
    unsigned int n;
    unsigned int nf;
    unsigned int factor[64];
    complex *twiddle[64];
    complex *trig;
    complex *scratch;
  }
gsl_fft_halfcomplex_wavetable;

int
  gsl_fft_halfcomplex_backwards (double *data,
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
  gsl_fft_halfcomplex_pass_2 (const double from[],
			      double to[],
			      const unsigned int product,
			      const unsigned int n,
			      const complex twiddle[]);

int
  gsl_fft_halfcomplex_pass_3 (const double from[], double to[],
			      const unsigned int product,
			      const unsigned int n,
			      const complex twiddle1[],
			      const complex twiddle2[]);

int
  gsl_fft_halfcomplex_pass_4 (const double from[],
			      double to[],
			      const unsigned int product,
			      const unsigned int n,
			      const complex twiddle1[],
			      const complex twiddle2[],
			      const complex twiddle3[]);

int
  gsl_fft_halfcomplex_pass_5 (const double from[],
			      double to[],
			      const unsigned int product,
			      const unsigned int n,
			      const complex twiddle1[],
			      const complex twiddle2[],
			      const complex twiddle3[],
			      const complex twiddle4[]);

int
  gsl_fft_halfcomplex_pass_6 (double *from, double *to,
			      unsigned int product, unsigned int n,
			      complex * twiddle1, complex * twiddle2,
			      complex * twiddle3, complex * twiddle4,
			      complex * twiddle5);

int
  gsl_fft_halfcomplex_pass_n (const double from[],
			      double to[],
			      const unsigned int factor,
			      const unsigned int product,
			      const unsigned int n,
			      const complex twiddle[]);


int
  gsl_fft_halfcomplex_unpack (const double halfcomplex_coefficient[],
			      complex complex_coefficient[],
			      const unsigned int n);


#endif /* _GSL_FFT_HALFCOMPLEX_H */
