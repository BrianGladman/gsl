int gsl_fft_complex_bitreverse_order (gsl_complex data[], 
				      const unsigned int n,
				      const unsigned int logn) ;

int gsl_fft_real_bitreverse_order (double data[], 
				   const unsigned int n,
				   const unsigned int logn) ;


int gsl_fft_complex_goldrader_bitreverse_order (gsl_complex data[], 
						const unsigned int n) ;

int gsl_fft_real_goldrader_bitreverse_order (double data[], 
					     const unsigned int n) ;

int gsl_fft_complex_rodriguez_bitreverse_order (gsl_complex data[], 
						const unsigned int n,
						const unsigned int logn) ;

int gsl_fft_real_rodriguez_bitreverse_order (double data[], 
					     const unsigned int n,
					     const unsigned int logn) ;
