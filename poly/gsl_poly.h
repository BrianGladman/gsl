#ifndef __GSL_POLY_H__
#define __GSL_POLY_H__

#include <stdlib.h>
#include <gsl/gsl_complex.h>

/* Solve for real or complex roots of the standard quadratic equation,
 * returning the number of real roots.
 *
 * Roots are returned ordered.
 */
int gsl_poly_solve_quadratic (double a, double b, double c, 
			      double * x0, double * x1);

int 
gsl_poly_complex_solve_quadratic (double a, double b, double c, 
				  gsl_complex * z0, gsl_complex * z1);


/* Solve for real roots of the cubic equation
 * x^3 + a x^2 + b x + c = 0, returning the
 * number of real roots.
 *
 * Roots are returned ordered.
 */
int gsl_poly_solve_cubic (double a, double b, double c, 
			  double * x0, double * x1, double * x2);

int 
gsl_poly_complex_solve_cubic (double a, double b, double c, 
			      gsl_complex * z0, gsl_complex * z1, 
			      gsl_complex * z2);


/* Solve for the complex roots of a general real polynomial */

typedef struct 
{ 
  size_t nc ;
  double * matrix ; 
} 
gsl_poly_complex_workspace ;

gsl_poly_complex_workspace * gsl_poly_complex_workspace_alloc (size_t n);
void gsl_poly_complex_workspace_free (gsl_poly_complex_workspace * w);

int
gsl_poly_complex_solve (const double * a, size_t n, 
                        gsl_poly_complex_workspace * w,
                        gsl_complex_packed_ptr z);

#endif /* __GSL_POLY_H__ */
