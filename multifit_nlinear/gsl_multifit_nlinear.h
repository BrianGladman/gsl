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

/* scaling matrix specification */
typedef enum
{
  GSL_MULTIFIT_NLINEAR_SCALE_LEVENBERG, /* D^T D = I */
  GSL_MULTIFIT_NLINEAR_SCALE_MARQUARDT, /* D^T D = diag(J^T J) */
  GSL_MULTIFIT_NLINEAR_SCALE_MORE       /* D^T D = max(diag(J^T J)) */
} gsl_multifit_nlinear_scale_t;

typedef enum
{
  GSL_MULTIFIT_NLINEAR_SOLVER_NORMAL, /* normal equations approach */
  GSL_MULTIFIT_NLINEAR_SOLVER_QR      /* QR decomposition of J */
} gsl_multifit_nlinear_solver_t;

/* tunable parameters for Levenberg-Marquardt method */
typedef struct
{
  gsl_multifit_nlinear_scale_t scale;   /* scaling method */
  gsl_multifit_nlinear_solver_t solver; /* solver method */
} gsl_multifit_nlinear_parameters;

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
} gsl_multifit_nlinear_fdf;

typedef struct
{
  const char *name;
  void * (*alloc) (const gsl_multifit_nlinear_parameters * params,
                   const size_t n, const size_t p);
  int (*init) (void *state, const gsl_vector * wts,
               gsl_multifit_nlinear_fdf * fdf, gsl_vector * x,
               gsl_vector * f, gsl_matrix * J);
  int (*iterate) (void *state, const gsl_vector * wts,
                  gsl_multifit_nlinear_fdf * fdf, gsl_vector * x,
                  gsl_vector * f, gsl_matrix * J, gsl_vector * dx);
  int (*gradient) (void *state, gsl_vector * g);
  void (*free) (void *state);
} gsl_multifit_nlinear_type;

typedef struct
{
  const gsl_multifit_nlinear_type * type;
  gsl_multifit_nlinear_fdf * fdf ;
  gsl_vector * x;        /* parameter values x */
  gsl_vector * f;        /* residual vector f(x) */
  gsl_vector * dx;       /* step dx */
  gsl_vector * g;        /* gradient J^T f */
  gsl_matrix * J;        /* Jacobian J(x) */
  gsl_vector * sqrt_wts; /* sqrt(W) */
  size_t niter;          /* number of iterations performed */
  void *state;
} gsl_multifit_nlinear_workspace;


gsl_multifit_nlinear_workspace *
gsl_multifit_nlinear_alloc (const gsl_multifit_nlinear_type * T,
                            const gsl_multifit_nlinear_parameters * params,
                            size_t n, size_t p);

void gsl_multifit_nlinear_free (gsl_multifit_nlinear_workspace * w);

int
gsl_multifit_nlinear_init (gsl_multifit_nlinear_fdf * fdf,
                           const gsl_vector * x,
                           gsl_multifit_nlinear_workspace * w);

int gsl_multifit_nlinear_winit (gsl_multifit_nlinear_fdf * f, 
                                const gsl_vector * x,
                                const gsl_vector * wts,
                                gsl_multifit_nlinear_workspace * w);

int
gsl_multifit_nlinear_iterate (gsl_multifit_nlinear_workspace * w);

int gsl_multifit_nlinear_driver (const size_t maxiter,
                                 const double xtol,
                                 const double gtol,
                                 const double ftol,
                                 int *info,
                                 gsl_multifit_nlinear_workspace * w);

gsl_matrix *
gsl_multifit_nlinear_jac (const gsl_multifit_nlinear_workspace * w);

const char *
gsl_multifit_nlinear_name (const gsl_multifit_nlinear_workspace * w);

gsl_vector *
gsl_multifit_nlinear_position (const gsl_multifit_nlinear_workspace * w);

gsl_vector *
gsl_multifit_nlinear_residual (const gsl_multifit_nlinear_workspace * w);

size_t
gsl_multifit_nlinear_niter (const gsl_multifit_nlinear_workspace * w);

int gsl_multifit_nlinear_eval_f(gsl_multifit_nlinear_fdf *fdf,
                                const gsl_vector *x,
                                const gsl_vector *swts,
                                gsl_vector *y);

int gsl_multifit_nlinear_eval_df(gsl_multifit_nlinear_fdf *fdf,
                                 const gsl_vector *x,
                                 const gsl_vector *f,
                                 const gsl_vector *swts,
                                 gsl_matrix *df);

/* covar.c */
int
gsl_multifit_nlinear_covar (const gsl_matrix * J, const double epsrel,
                            gsl_matrix * covar);

/* convergence.c */
int
gsl_multifit_nlinear_test (const double xtol, const double gtol,
                           const double ftol, int *info,
                           const gsl_multifit_nlinear_workspace * w);

/* fdjac.c */
int
gsl_multifit_nlinear_df(const gsl_vector *x, const gsl_vector *wts,
                        gsl_multifit_nlinear_fdf *fdf,
                        const gsl_vector *f, gsl_matrix *J);

GSL_VAR const gsl_multifit_nlinear_type * gsl_multifit_nlinear_lm;

__END_DECLS

#endif /* __GSL_MULTIFIT_NLINEAR_H__ */
