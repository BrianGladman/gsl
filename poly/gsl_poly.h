#ifndef GSL_POLY_H
#define GSL_POLY_H

#include <gsl_complex.h>

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

#endif /* GSL_POLY_H */
