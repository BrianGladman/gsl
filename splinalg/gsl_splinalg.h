/* gsl_splinalg.h
 * 
 * Copyright (C) 2012-2014 Patrick Alken
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_SPLINALG_H__
#define __GSL_SPLINALG_H__

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spmatrix.h>

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

typedef struct
{
  size_t n;        /* size of linear system */
  size_t m;        /* dimension of Krylov subspace K_m */
  gsl_vector *r;   /* residual vector r = b - A*x */
  gsl_matrix *H;   /* Hessenberg matrix n-by-(m+1) */
  gsl_vector *tau; /* householder scalars */
  gsl_vector *y;   /* least squares rhs and solution vector */

  double *c;       /* Givens rotations */
  double *s;

  double normr;    /* residual norm ||r|| */
} gsl_splinalg_gmres_workspace;

/*
 * Prototypes
 */

gsl_splinalg_gmres_workspace *gsl_splinalg_gmres_alloc(const size_t n,
                                                       const size_t krylov_m);
void gsl_splinalg_gmres_free(gsl_splinalg_gmres_workspace *w);
int gsl_splinalg_gmres_solve(const gsl_spmatrix *A, const gsl_vector *b,
                             gsl_vector *x,
                             gsl_splinalg_gmres_workspace *w);
int gsl_splinalg_gmres_solve_x(const gsl_spmatrix *A,
                               const gsl_vector *b, const double tol,
                               gsl_vector *x,
                               gsl_splinalg_gmres_workspace *w);

__END_DECLS

#endif /* __GSL_SPLINALG_H__ */
