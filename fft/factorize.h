int
  gsl_fft_complex_factorize (const unsigned int n,
			     unsigned int *nf,
			     unsigned int factors[]);

int
  gsl_fft_halfcomplex_factorize (const unsigned int n,
				 unsigned int *nf,
				 unsigned int factors[]);

int
  gsl_fft_real_factorize (const unsigned int n,
			  unsigned int *nf,
			  unsigned int factors[]);

int
  gsl_fft_factorize (const unsigned int n,
		     const unsigned int implemented_subtransforms[],
		     unsigned int *n_factors,
		     unsigned int factors[]);

int gsl_fft_binary_logn (const unsigned int n) ;

