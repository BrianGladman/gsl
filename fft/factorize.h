int
  gsl_fft_complex_factorize (const size_t n,
			     size_t *nf,
			     size_t factors[]);

int
  gsl_fft_halfcomplex_factorize (const size_t n,
				 size_t *nf,
				 size_t factors[]);

int
  gsl_fft_real_factorize (const size_t n,
			  size_t *nf,
			  size_t factors[]);

int
  gsl_fft_factorize (const size_t n,
		     const size_t implemented_subtransforms[],
		     size_t *n_factors,
		     size_t factors[]);

int gsl_fft_binary_logn (const size_t n) ;

