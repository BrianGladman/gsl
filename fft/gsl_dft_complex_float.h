#ifndef __GSL_DFT_COMPLEX_FLOAT_H__
#define __GSL_DFT_COMPLEX_FLOAT_H__

#include <stddef.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_fft.h>

int gsl_dft_complex_float_forward (const float data[], size_t stride, size_t n,
			     float result[]);

int gsl_dft_complex_float_backward (const float data[], size_t stride, size_t n,
			      float result[]);

int gsl_dft_complex_float_inverse (const float data[], size_t stride, size_t n,
			     float result[]);

int gsl_dft_complex_float_transform (const float data[], size_t stride, size_t n,
		     float result[], const gsl_fft_direction sign);

#endif /* __GSL_DFT_COMPLEX_FLOAT_H__ */
