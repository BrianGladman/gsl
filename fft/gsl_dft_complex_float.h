#ifndef GSL_DFT_COMPLEX_FLOAT_H
#define GSL_DFT_COMPLEX_FLOAT_H

#include <stddef.h>

#include <gsl_math.h>
#include <gsl_complex.h>
#include <gsl_vector.h>
#include <gsl_fft.h>

int gsl_dft_complex_float_forward (const float data[], size_t stride, size_t n,
			     float result[]);

int gsl_dft_complex_float_backward (const float data[], size_t stride, size_t n,
			      float result[]);

int gsl_dft_complex_float_inverse (const float data[], size_t stride, size_t n,
			     float result[]);

int gsl_dft_complex_float_transform (const float data[], size_t stride, size_t n,
		     float result[], const gsl_fft_direction sign);

#endif /* GSL_DFT_COMPLEX_FLOAT_H */
