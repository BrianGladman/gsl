#include "complex_internal.h"

static void
FUNCTION(fft_halfcomplex,pass_2) (const BASE in[],
				  size_t istride,
				  BASE out[],
				  size_t ostride,
				  size_t product,
				  size_t n,
				  const gsl_complex twiddle[]);

static void
FUNCTION(fft_halfcomplex,pass_3) (const BASE in[], 
				  size_t istride,
				  BASE out[],
				  size_t ostride,
				  size_t product,
				  size_t n,
				  const gsl_complex twiddle1[],
				  const gsl_complex twiddle2[]);

static void
FUNCTION(fft_halfcomplex,pass_4) (const BASE in[],
				  size_t istride,
				  BASE out[],
				  size_t ostride,
				  size_t product,
				  size_t n,
				  const gsl_complex twiddle1[],
				  const gsl_complex twiddle2[],
				  const gsl_complex twiddle3[]);

static void
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

static void
FUNCTION(fft_halfcomplex,pass_n) (const BASE in[],
				  size_t istride,
				  BASE out[],
				  size_t ostride,
				  size_t factor,
				  size_t product,
				  size_t n,
				  const gsl_complex twiddle[]);




