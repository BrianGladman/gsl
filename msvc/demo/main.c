#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

/* See the file README.txt for information on compiling this program */

int
main (void)
{
	const gsl_rng_type * T = gsl_rng_default;

	gsl_rng * r = gsl_rng_alloc (T);

	int i;

	printf ("Here are ten random numbers in the range 0-99:\n");

	for (i = 0; i < 10; i++)
	{
		int k = gsl_rng_uniform_int (r, 100);
		printf(" %d", k);
	}

	printf("\n");

	{
		double x[5] = { 1.0, 2.0, 3.0, 4.0, 5.0} ;
		double y[5] = { 5.5, 4.4, 3.3, 2.2, 1.1} ; 
		gsl_vector_view v = gsl_vector_view_array(x, 5);
		gsl_vector_view w = gsl_vector_view_array(y, 5);

		printf("blas operation DAXPY\n");

		printf("x:"); for (i = 0; i < 5; i++) { printf(" %g", x[i]); } ; 
		printf("  y:"); for (i = 0; i < 5; i++) { printf(" %g", y[i]); } ;
		
		gsl_blas_daxpy (1.0, &v.vector, &w.vector);
		
		printf("  a x + y:"); for (i = 0; i < 5; i++) { printf(" %g", y[i]); } ;
	}

	printf("\n");

	return 0;
}
