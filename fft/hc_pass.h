#include "complex_internal.h"

void
FUNCTION(fft_halfcomplex,pass_2) (const BASE in[],
				  size_t istride,
				  BASE out[],
				  size_t ostride,
				  size_t product,
				  size_t n,
				  const gsl_complex twiddle[]);

void
FUNCTION(fft_halfcomplex,pass_3) (const BASE in[], 
				  size_t istride,
				  BASE out[],
				  size_t ostride,
				  size_t product,
				  size_t n,
				  const gsl_complex twiddle1[],
				  const gsl_complex twiddle2[]);

void
FUNCTION(fft_halfcomplex,pass_4) (const BASE in[],
				  size_t istride,
				  BASE out[],
				  size_t ostride,
				  size_t product,
				  size_t n,
				  const gsl_complex twiddle1[],
				  const gsl_complex twiddle2[],
				  const gsl_complex twiddle3[]);

void
FUNCTION(fft_halfcomplex,pass_5) (const BASE in[],
				  size_t istride,
				  BASE out[],
				  size_t ostride,
				  size_t product,
				  size_t n,
				  const gsl_complex twiddle1[],
				  const gsl_complex twiddle2[],
				  const gsl_complex twiddle3[],
				  const gsl_complex twiddle4[]);

void
FUNCTION(fft_halfcomplex,pass_6) (const BASE in[], 
				  size_t istride,
				  const BASE out[],
				  size_t ostride,
				  size_t product, size_t n,
				  gsl_complex twiddle1[], 
				  gsl_complex twiddle2[],
				  gsl_complex twiddle3[], 
				  gsl_complex twiddle4[],
				  gsl_complex twiddle5[]);

void
FUNCTION(fft_halfcomplex,pass_n) (const BASE in[],
				  size_t istride,
				  BASE out[],
				  size_t ostride,
				  size_t factor,
				  size_t product,
				  size_t n,
				  const gsl_complex twiddle[]);




