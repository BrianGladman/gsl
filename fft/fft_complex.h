#include "complex_internal.h"

int
fft_complex_pass_2 (const double in[],
			size_t istride,
			double out[],
			size_t ostride,
			gsl_fft_direction sign,
			size_t product,
			size_t n,
			const gsl_complex twiddle[]);

int
fft_complex_pass_3 (const double in[],
			size_t istride,
			double out[],
			size_t ostride,
			gsl_fft_direction sign,
			size_t product,
			size_t n,
			const gsl_complex twiddle1[],
			const gsl_complex twiddle2[]);

int
fft_complex_pass_4 (const double in[],
			size_t istride,
			double out[],
			size_t ostride,
			gsl_fft_direction sign,
			size_t product,
			size_t n,
			const gsl_complex twiddle1[],
			const gsl_complex twiddle2[],
			const gsl_complex twiddle3[]);

int
fft_complex_pass_5 (const double in[],
			size_t istride,
			double out[],
			size_t ostride,
			gsl_fft_direction sign,
			size_t product,
			size_t n,
			const gsl_complex twiddle1[],
			const gsl_complex twiddle2[],
			const gsl_complex twiddle3[],
			const gsl_complex twiddle4[]);

int
fft_complex_pass_6 (const double in[],
			size_t istride,
			double out[],
			size_t ostride,
			gsl_fft_direction sign,
			size_t product,
			size_t n,
			const gsl_complex twiddle1[],
			const gsl_complex twiddle2[],
			const gsl_complex twiddle3[],
			const gsl_complex twiddle4[],
			const gsl_complex twiddle5[]);

int
fft_complex_pass_7 (const double in[],
			size_t istride,
			double out[],
			size_t ostride,
			gsl_fft_direction sign,
			size_t product,
			size_t n,
			const gsl_complex twiddle1[],
			const gsl_complex twiddle2[],
			const gsl_complex twiddle3[],
			const gsl_complex twiddle4[],
			const gsl_complex twiddle5[],
			const gsl_complex twiddle6[]);


int
fft_complex_pass_n (double in[],
			size_t istride,
			double out[],
			size_t ostride,
			gsl_fft_direction sign,
			size_t factor,
			size_t product,
			size_t n,
			const gsl_complex twiddle[]);

