#ifndef __GSL_DFT_COMPLEX_H__
#define __GSL_DFT_COMPLEX_H__

#include <stddef.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_fft.h>

int gsl_dft_complex_forward (const double data[], size_t stride, size_t n,
			     double result[]);

int gsl_dft_complex_backward (const double data[], size_t stride, size_t n,
			      double result[]);

int gsl_dft_complex_inverse (const double data[], size_t stride, size_t n,
			     double result[]);

int gsl_dft_complex_transform (const double data[], size_t stride, size_t n,
		     double result[], const gsl_fft_direction sign);

#endif /* __GSL_DFT_COMPLEX_H__ */
