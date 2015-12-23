/* multifit_nlinear/gsl_multifit_nlinear.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * Copyright (C) 2015 Patrick Alken
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

#ifndef __GSL_MULTIFIT_NLINEAR_H__
#define __GSL_MULTIFIT_NLINEAR_H__

#include <stdlib.h>
#include <gsl/gsl_types.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>

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

/* Definition of vector-valued functions and gradient with parameters
   based on gsl_vector */

typedef struct
{
  int (* f) (const gsl_vector * x, void * params, gsl_vector * f);
  int (* df) (const gsl_vector * x, void * params, gsl_matrix * df);
  size_t n;       /* number of functions */
  size_t p;       /* number of independent variables */
  void * params;  /* user parameters */
  size_t nevalf;  /* number of function evaluations */
  size_t nevaldf; /* number of Jacobian evaluations */
} gsl_multifit_function_fdf;

typedef struct
{
  const char *name;
  size_t size;
  int (*alloc) (void *state, size_t n, size_t p);
  int (*set) (void *state, const gsl_vector * wts,
              gsl_multifit_function_fdf * fdf, gsl_vector * x,
              gsl_vector * f, gsl_vector * dx);
  int (*iterate) (void *state, const gsl_vector * wts,
                  gsl_multifit_function_fdf * fdf, gsl_vector * x,
                  gsl_vector * f, gsl_vector * dx);
  int (*gradient) (void *state, gsl_vector * g);
  int (*jac) (void *state, gsl_matrix * J);
  void (*free) (void *state);
} gsl_multifit_nlinear_type;

typedef struct
{
  const gsl_multifit_nlinear_type * type;
  gsl_multifit_function_fdf * fdf ;
  gsl_vector * x;        /* parameter values x */
  gsl_vector * f;        /* residual vector f(x) */
  gsl_vector * dx;       /* step dx */
  gsl_vector * g;        /* gradient J^T f */
  gsl_vector * sqrt_wts; /* sqrt(wts) */
  size_t niter;          /* number of iterations performed */
  void *state;
} gsl_multifit_nlinear_workspace;


gsl_multifit_nlinear_workspace *
gsl_multifit_nlinear_alloc (const gsl_multifit_nlinear_type * T, 
                            size_t n, size_t p);

int
gsl_multifit_nlinear_set (gsl_multifit_nlinear_workspace * s, 
                          gsl_multifit_function_fdf * fdf,
                          const gsl_vector * x);
int gsl_multifit_nlinear_wset (gsl_multifit_nlinear_workspace * s, 
                               gsl_multifit_function_fdf * f, 
                               const gsl_vector * x,
                               const gsl_vector * wts);

int
gsl_multifit_nlinear_iterate (gsl_multifit_nlinear_workspace * s);

int gsl_multifit_nlinear_driver (gsl_multifit_nlinear_workspace * s,
                                 const size_t maxiter,
                                 const double xtol,
                                 const double gtol,
                                 const double ftol,
                                 int *info);

int gsl_multifit_nlinear_jac (gsl_multifit_nlinear_workspace * s,
                              gsl_matrix * J);

void
gsl_multifit_nlinear_free (gsl_multifit_nlinear_workspace * s);

const char * gsl_multifit_nlinear_name (const gsl_multifit_nlinear_workspace * s);
gsl_vector * gsl_multifit_nlinear_position (const gsl_multifit_nlinear_workspace * s);
gsl_vector * gsl_multifit_nlinear_residual (const gsl_multifit_nlinear_workspace * s);
size_t gsl_multifit_nlinear_niter (const gsl_multifit_nlinear_workspace * s);
int gsl_multifit_eval_wf(gsl_multifit_function_fdf *fdf,
                         const gsl_vector *x, const gsl_vector *wts,
                         gsl_vector *y);
int gsl_multifit_eval_wdf(gsl_multifit_function_fdf *fdf,
                          const gsl_vector *x, const gsl_vector *wts,
                          gsl_matrix *dy);

int gsl_multifit_nlinear_test (const gsl_multifit_nlinear_workspace * s,
                               const double xtol,
                               const double gtol,
                               const double ftol, int *info);
int gsl_multifit_test_delta (const gsl_vector * dx, const gsl_vector * x, 
                             double epsabs, double epsrel);

int gsl_multifit_test_gradient (const gsl_vector * g, double epsabs);

int gsl_multifit_nlinear_dif_df(const gsl_vector *x,
                                const gsl_vector *wts,
                                gsl_multifit_function_fdf *fdf,
                                const gsl_vector *f, gsl_matrix *J);

GSL_VAR const gsl_multifit_nlinear_type * gsl_multifit_nlinear_lmsder;
GSL_VAR const gsl_multifit_nlinear_type * gsl_multifit_nlinear_lmder;
GSL_VAR const gsl_multifit_nlinear_type * gsl_multifit_nlinear_lmniel;

__END_DECLS

#endif /* __GSL_MULTIFIT_NLINEAR_H__ */
