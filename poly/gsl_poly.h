#ifndef GSL_POLY_H
#define GSL_POLY_H

#include <gsl_complex.h>

/* Solve for real or complex roots of the standard quadratic equation,
 * returning the number of real roots.
 * x[] is assumed big enough.
 * Roots are returned ordered.
 */
int gsl_root_solve_quadratic (double a, double b, double c, double x[]);

int 
gsl_root_complex_solve_quadratic (double a, double b, double c, gsl_complex z[]);


/* Solve for real roots of the cubic equation
 * x^3 + a x^2 + b x + c = 0, returning the
 * number of real roots.
 * x[] is assumed big enough.
 * Roots are returned ordered.
 */
int gsl_root_solve_cubic (double a, double b, double c, double x[]);

int 
gsl_root_complex_solve_cubic (double a, double b, double c, gsl_complex z[]);

#endif /* GSL_POLY_H */
