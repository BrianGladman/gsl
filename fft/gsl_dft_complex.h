#ifndef GSL_DFT_COMPLEX_H
#define GSL_DFT_COMPLEX_H

#include <stddef.h>

#include <gsl_math.h>
#include <gsl_complex.h>
#include <gsl_fft.h>

int gsl_dft_complex_forward (const double data[], size_t stride, size_t n,
			     double result[]);

int gsl_dft_complex_backward (const double data[], size_t stride, size_t n,
			      double result[]);

int gsl_dft_complex_inverse (const double data[], size_t stride, size_t n,
			     double result[]);

int gsl_dft_complex_transform (const double data[], size_t stride, size_t n,
		     double result[], const gsl_fft_direction sign);

#endif /* GSL_DFT_COMPLEX_H */
