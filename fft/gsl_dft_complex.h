#ifndef __GSL_DFT_COMPLEX_H__
#define __GSL_DFT_COMPLEX_H__

#include <stddef.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_fft.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

int gsl_dft_complex_forward (const double data[], size_t stride, size_t n,
			     double result[]);

int gsl_dft_complex_backward (const double data[], size_t stride, size_t n,
			      double result[]);

int gsl_dft_complex_inverse (const double data[], size_t stride, size_t n,
			     double result[]);

int gsl_dft_complex_transform (const double data[], size_t stride, size_t n,
		     double result[], const gsl_fft_direction sign);

__END_DECLS

#endif /* __GSL_DFT_COMPLEX_H__ */
