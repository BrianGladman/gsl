#include <gsl_complex.h>

int
compare_complex_results (const char *name_a, const gsl_complex a[],
			 const char *name_b, const gsl_complex b[],
			 size_t n,
			 const double allowed_ticks);
