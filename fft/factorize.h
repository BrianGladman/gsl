int
fft_complex_factorize (const size_t n,
			   size_t *nf,
			   size_t factors[]);

int
fft_halfcomplex_factorize (const size_t n,
			       size_t *nf,
			       size_t factors[]);

int
fft_real_factorize (const size_t n,
			size_t *nf,
			size_t factors[]);

int
fft_factorize (const size_t n,
		   const size_t implemented_subtransforms[],
		   size_t *n_factors,
		   size_t factors[]);

int fft_binary_logn (const size_t n) ;

