#ifndef GSL_DFT_COMPLEX_H
#define GSL_DFT_COMPLEX_H

#include <stddef.h>

#include <gsl_math.h>
#include <gsl_complex.h>
#include <gsl_fft.h>

int
  gsl_dft_complex_forward (const gsl_complex data[],
			   gsl_complex result[],
			   const size_t n);

int
  gsl_dft_complex_backward (const gsl_complex data[],
			    gsl_complex result[],
			    const size_t n);

int
  gsl_dft_complex_inverse (const gsl_complex data[],
			   gsl_complex result[],
			   const size_t n);

int
  gsl_dft_complex (const gsl_complex data[],
		   gsl_complex result[],
		   const size_t n,
		   const gsl_fft_direction sign);


#endif /* GSL_DFT_COMPLEX_H */
