#ifndef _GSL_FFT_H
#define _GSL_FFT_H

#include <gsl_complex.h>

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

int gsl_fft_complex_bitreverse_order (complex data[], 
				      const unsigned int n,
				      const unsigned int logn) ;

int gsl_fft_real_bitreverse_order (double data[], 
				   const unsigned int n,
				   const unsigned int logn) ;


int gsl_fft_complex_goldrader_bitreverse_order (complex data[], 
						const unsigned int n) ;

int gsl_fft_real_goldrader_bitreverse_order (double data[], 
					     const unsigned int n) ;

int gsl_fft_complex_rodriguez_bitreverse_order (complex data[], 
						const unsigned int n,
						const unsigned int logn) ;

int gsl_fft_real_rodriguez_bitreverse_order (double data[], 
					     const unsigned int n,
					     const unsigned int logn) ;

typedef enum
  {
    forward = -1, backward = +1   
  }
gsl_fft_direction;

/* this give the sign in the formula

   h(f) = \sum x(t) exp(+/- 2 pi i f t) 
       
   where - is the forward transform direction and + the inverse direction */

#endif /* _GSL_FFT_H */
