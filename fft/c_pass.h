#include "complex_internal.h"

int
FUNCTION(fft_complex,pass_2) (const BASE in[],
			      size_t istride,
			      BASE out[],
			      size_t ostride,
			      gsl_fft_direction sign,
			      size_t product,
			      size_t n,
			      const gsl_complex twiddle[]);

int
FUNCTION(fft_complex,pass_3) (const BASE in[],
			      size_t istride,
			      BASE out[],
			      size_t ostride,
			      gsl_fft_direction sign,
			      size_t product,
			      size_t n,
			      const gsl_complex twiddle1[],
			      const gsl_complex twiddle2[]);

int
FUNCTION(fft_complex,pass_4) (const BASE in[],
			      size_t istride,
			      BASE out[],
			      size_t ostride,
			      gsl_fft_direction sign,
			      size_t product,
			      size_t n,
			      const gsl_complex twiddle1[],
			      const gsl_complex twiddle2[],
			      const gsl_complex twiddle3[]);

int
FUNCTION(fft_complex,pass_5) (const BASE in[],
			      size_t istride,
			      BASE out[],
			      size_t ostride,
			      gsl_fft_direction sign,
			      size_t product,
			      size_t n,
			      const gsl_complex twiddle1[],
			      const gsl_complex twiddle2[],
			      const gsl_complex twiddle3[],
			      const gsl_complex twiddle4[]);

int
FUNCTION(fft_complex,pass_6) (const BASE in[],
			      size_t istride,
			      BASE out[],
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
FUNCTION(fft_complex,pass_7) (const BASE in[],
			      size_t istride,
			      BASE out[],
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
FUNCTION(fft_complex,pass_n) (BASE in[],
			      size_t istride,
			      BASE out[],
			      size_t ostride,
			      gsl_fft_direction sign,
			      size_t factor,
			      size_t product,
			      size_t n,
			      const gsl_complex twiddle[]);

