/* poly/gsl_poly.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef __GSL_POLY_H__
#define __GSL_POLY_H__

#include <stdlib.h>
#include <gsl/gsl_complex.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

/* Struct for a polynomial */
struct gsl_poly_struct {

  double * c;    /* coefficients   */
  size_t size;   /* allocated size */
};
typedef struct gsl_poly_struct gsl_poly;
 
/* Allocate polynomial */
gsl_poly * gsl_poly_alloc (const size_t size);

/* Allocate polynomial and init to zero polynomial*/
gsl_poly * gsl_poly_calloc (const size_t size);

/* Free polynomial */
void gsl_poly_free (gsl_poly * p);

/* Operations */
double gsl_poly_get (const gsl_poly * p, const size_t i);
void gsl_poly_set (gsl_poly * p, const size_t i, double a);

int gsl_poly_set_zero (gsl_poly * p);
int gsl_poly_set_all (gsl_poly * p, size_t d, double x);
int gsl_poly_set_basis (gsl_poly * p, size_t i);
int gsl_poly_scale (gsl_poly * p, double a);
int gsl_poly_consistent (const gsl_poly * p, double tol);
 
int gsl_poly_memcpy (gsl_poly * dest, const gsl_poly * src);

int gsl_poly_add (gsl_poly * p1, const gsl_poly * p2);
int gsl_poly_sub (gsl_poly * p1, const gsl_poly * p2);
int gsl_poly_mul (gsl_poly * q, const gsl_poly * p1, const gsl_poly * p2);
int gsl_poly_div (gsl_poly * q, gsl_poly *r, const gsl_poly * u, const gsl_poly * v);

int gsl_poly_diff (gsl_poly * dp, const gsl_poly * p);

double gsl_poly_eval2 (gsl_poly * p, double x);

void gsl_poly_dump (gsl_poly * p);

/* Evaluate polynomial
 *
 * c[0] + c[1] x + c[2] x^2 + ... + c[len-1] x^(len-1)
 *
 * exceptions: none
 */
double gsl_poly_eval(const double c[], const int len, const double x);


#ifdef HAVE_INLINE
extern inline
double gsl_poly_eval(const double c[], const int len, const double x)
{
  int i;
  double ans = c[len-1];
  for(i=len-1; i>0; i--) ans = c[i-1] + x * ans;
  return ans;
}
#endif /* HAVE_INLINE */

/* Work with divided-difference polynomials, Abramowitz & Stegun 25.2.26 */

int
gsl_poly_dd_init (double dd[], const double x[], const double y[],
                  size_t size);

double
gsl_poly_dd_eval (const double dd[], const double xa[], const size_t size, const double x);

#ifdef HAVE_INLINE
extern inline
double gsl_poly_dd_eval(const double dd[], const double xa[], const size_t size, const double x)
{
  size_t i;
  double y = dd[size - 1];
  for (i = size - 1; i--;) y = dd[i] + (x - xa[i]) * y;
  return y;
}
#endif /* HAVE_INLINE */


int
gsl_poly_dd_taylor (double c[], double xp,
                    const double dd[], const double x[], size_t size,
                    double w[]);

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

__END_DECLS

#endif /* __GSL_POLY_H__ */
