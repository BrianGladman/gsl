#include "fft.h"

int
gsl_fft_halfcomplex_pass_2 (const double in[],
			    size_t istride,
			    double out[],
			    size_t ostride,
			    size_t product,
			    size_t n,
			    const gsl_complex twiddle[]);

int
gsl_fft_halfcomplex_pass_3 (const double in[], 
			    size_t istride,
			    double out[],
			    size_t ostride,
			    size_t product,
			    size_t n,
			    const gsl_complex twiddle1[],
			    const gsl_complex twiddle2[]);

int
gsl_fft_halfcomplex_pass_4 (const double in[],
			    size_t istride,
			    double out[],
			    size_t ostride,
			    size_t product,
			    size_t n,
			    const gsl_complex twiddle1[],
			    const gsl_complex twiddle2[],
			    const gsl_complex twiddle3[]);

int
gsl_fft_halfcomplex_pass_5 (const double in[],
			    size_t istride,
			    double out[],
			    size_t ostride,
			    size_t product,
			    size_t n,
			    const gsl_complex twiddle1[],
			    const gsl_complex twiddle2[],
			    const gsl_complex twiddle3[],
			    const gsl_complex twiddle4[]);

int
gsl_fft_halfcomplex_pass_6 (const double in[], 
			    size_t istride,
			    const double out[],
			    size_t ostride,
			    size_t product, size_t n,
			    gsl_complex * twiddle1, gsl_complex * twiddle2,
			    gsl_complex * twiddle3, gsl_complex * twiddle4,
			    gsl_complex * twiddle5);

int
gsl_fft_halfcomplex_pass_n (const double in[],
			    size_t istride,
			    double out[],
			    size_t ostride,
			    size_t factor,
			    size_t product,
			    size_t n,
			    const gsl_complex twiddle[]);




