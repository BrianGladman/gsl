#include <gsl_complex.h>

int
FUNCTION(compare_complex,results) (const char *name_a, const BASE a[],
				   const char *name_b, const BASE b[],
				   size_t n, size_t stride, 
				   const double allowed_ticks);
