#ifndef GSL_DFT_COMPLEX_H
#define GSL_DFT_COMPLEX_H

#include <stddef.h>

#include <gsl_math.h>
#include <gsl_complex.h>
#include <gsl_vector.h>
#include <gsl_fft.h>

int
gsl_dft_complex_forward (const gsl_vector_complex data[],
			 gsl_vector_complex result[]);

int
gsl_dft_complex_backward (const gsl_vector_complex data[],
			  gsl_vector_complex result[]);

int
gsl_dft_complex_inverse (const gsl_vector_complex data[],
			 gsl_vector_complex result[]);

int
gsl_dft_complex (const gsl_vector_complex data[],
		 gsl_vector_complex result[],
		 const gsl_fft_direction sign);

#endif /* GSL_DFT_COMPLEX_H */
