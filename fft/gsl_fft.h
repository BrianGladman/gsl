





#ifndef _GSL_FFT_H
#define _GSL_FFT_H

#include <gsl_complex.h>

int
  gsl_fft_complex_factorize (const unsigned int n,
			     unsigned int *nf,
			     unsigned int factors[]);

int
  gsl_fft_halfcomplex_factorize (const unsigned int n,
				 unsigned int *nf,
				 unsigned int factors[]);

int
  gsl_fft_real_factorize (const unsigned int n,
			  unsigned int *nf,
			  unsigned int factors[]);

int
  gsl_fft_factorize (const unsigned int n,
		     const unsigned int implemented_subtransforms[],
		     unsigned int *n_factors,
		     unsigned int factors[]);

typedef enum
  {
    forward = +1, backward = -1
  }
gsl_fft_direction;

#endif /* _GSL_FFT_H */
