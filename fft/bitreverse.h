#include <stddef.h>

int bitreverse_order_complex (double data[], 
			      size_t stride,
			      size_t n,
			      size_t logn) ;

int bitreverse_order_real (double data[], 
			   const size_t n,
			   const size_t logn) ;


int goldrader_bitreverse_order_complex (gsl_complex data[], 
					const size_t n) ;


int rodriguez_bitreverse_order_complex (gsl_complex data[], 
					const size_t n,
					const size_t logn) ;

