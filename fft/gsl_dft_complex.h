#ifndef _GSL_DFT_COMPLEX_H
#define _GSL_DFT_COMPLEX_H

#include <gsl_math.h>

#include <gsl_complex.h>
#include <gsl_fft.h>

int
  gsl_dft_complex_forward (const gsl_complex data[],
			   gsl_complex result[],
			   const unsigned int n);

int
  gsl_dft_complex_backward (const gsl_complex data[],
			    gsl_complex result[],
			    const unsigned int n);

int
  gsl_dft_complex_inverse (const gsl_complex data[],
			   gsl_complex result[],
			   const unsigned int n);

int
  gsl_dft_complex (const gsl_complex data[],
		   gsl_complex result[],
		   const unsigned int n,
		   const gsl_fft_direction sign);


#endif /* _GSL_DFT_COMPLEX_H */
