#include <gsl_complex.h>

int
  compare_complex_results (const char *name_a, const complex a[],
			   const char *name_b, const complex b[],
			   unsigned int n,
			   const double allowed_ticks);
