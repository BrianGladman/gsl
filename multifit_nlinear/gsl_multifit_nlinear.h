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

typedef enum
{
  GSL_MULTIFIT_NLINEAR_FWDIFF,
  GSL_MULTIFIT_NLINEAR_CTRDIFF
} gsl_multifit_nlinear_fdtype;

/* Definition of vector-valued functions and gradient with parameters
   based on gsl_vector */

typedef struct
{
  int (* f) (const gsl_vector * x, void * params, gsl_vector * f);
  int (* df) (const gsl_vector * x, void * params, gsl_matrix * df);
  int (* fvv) (const gsl_vector * x, const gsl_vector * v, void * params,
               gsl_vector * fvv);
  size_t n;        /* number of functions */
  size_t p;        /* number of independent variables */
  void * params;   /* user parameters */
  size_t nevalf;   /* number of function evaluations */
  size_t nevaldf;  /* number of Jacobian evaluations */
  size_t nevalfvv; /* number of fvv evaluations */
} gsl_multifit_nlinear_fdf;

/* LM parameter updating method */
typedef struct
{
  const char *name;
  void * (*alloc) (void);
  int (*init) (const gsl_matrix * J, const gsl_vector * diag,
               const gsl_vector * x, double * mu, void * vstate);
  int (*accept) (const double rho, double * mu, void * vstate);
  int (*reject) (double * mu, void * vstate);
  void (*free) (void * vstate);
} gsl_multifit_nlinear_update;

/* scaling matrix specification */
typedef struct
{
  const char *name;
  int (*init) (const gsl_matrix * J, gsl_vector * diag);
  int (*update) (const gsl_matrix * J, gsl_vector * diag);
} gsl_multifit_nlinear_scale;

/* linear least squares solvers */
typedef struct
{
  const char *name;
  void * (*alloc) (const size_t n, const size_t p);
  int (*init) (const gsl_vector * f, const gsl_matrix * J,
               const gsl_vector * g, void * vstate);
  int (*init_mu) (const double mu, const gsl_vector * diag,
                  void * vstate);
  int (*solve_vel) (gsl_vector * v, void * vstate);
  int (*solve_acc) (const gsl_matrix * J, const gsl_vector * fvv,
                    gsl_vector * a, void * vstate);
  void (*free) (void * vstate);
} gsl_multifit_nlinear_solver;

/* tunable parameters for Levenberg-Marquardt method */
typedef struct
{
  const gsl_multifit_nlinear_update *update;  /* LM parameter update method */
  const gsl_multifit_nlinear_scale *scale;    /* scaling method */
  const gsl_multifit_nlinear_solver *solver;  /* solver method */
  gsl_multifit_nlinear_fdtype fdtype;         /* finite difference method */
  int accel;                                  /* use geodesic acceleration */
  double avmax;                               /* max allowed |a|/|v| */
  double h_df;                                /* step size for finite difference Jacobian */
  double h_fvv;                               /* step size for finite difference fvv */
} gsl_multifit_nlinear_parameters;

typedef struct
{
  const char *name;
  void * (*alloc) (const gsl_multifit_nlinear_parameters * params,
                   const size_t n, const size_t p);
  int (*init) (void * state, const gsl_vector * wts,
               gsl_multifit_nlinear_fdf * fdf, const gsl_vector * x,
               gsl_vector * f, gsl_matrix * J, gsl_vector * g);
  int (*iterate) (void * state, const gsl_vector * wts,
                  gsl_multifit_nlinear_fdf * fdf, gsl_vector * x,
                  gsl_vector * f, gsl_matrix * J, gsl_vector * g,
                  gsl_vector * dx);
  int (*rcond) (const gsl_matrix * J, double * rcond, void * state);
  double (*avratio) (void * state);
  void (*free) (void * state);
} gsl_multifit_nlinear_type;

typedef struct
{
  const gsl_multifit_nlinear_type * type;
  gsl_multifit_nlinear_fdf * fdf ;
  gsl_vector * x;             /* parameter values x */
  gsl_vector * f;             /* residual vector f(x) */
  gsl_vector * dx;            /* step dx */
  gsl_vector * g;             /* gradient J^T f */
  gsl_matrix * J;             /* Jacobian J(x) */
  gsl_vector * sqrt_wts_work; /* sqrt(W) */
  gsl_vector * sqrt_wts;      /* ptr to sqrt_wts_work, or NULL if not using weights */
  size_t niter;               /* number of iterations performed */
  void *state;
} gsl_multifit_nlinear_workspace;

gsl_multifit_nlinear_workspace *
gsl_multifit_nlinear_alloc (const gsl_multifit_nlinear_type * T,
                            const gsl_multifit_nlinear_parameters * params,
                            size_t n, size_t p);

void gsl_multifit_nlinear_free (gsl_multifit_nlinear_workspace * w);

gsl_multifit_nlinear_parameters gsl_multifit_nlinear_default_parameters(void);

int
gsl_multifit_nlinear_init (const gsl_vector * x,
                           gsl_multifit_nlinear_fdf * fdf,
                           gsl_multifit_nlinear_workspace * w);

int gsl_multifit_nlinear_winit (const gsl_vector * x,
                                const gsl_vector * wts,
                                gsl_multifit_nlinear_fdf * fdf,
                                gsl_multifit_nlinear_workspace * w);

int
gsl_multifit_nlinear_iterate (gsl_multifit_nlinear_workspace * w);

double
gsl_multifit_nlinear_avratio (const gsl_multifit_nlinear_workspace * w);

int
gsl_multifit_nlinear_driver (const size_t maxiter,
                             const double xtol,
                             const double gtol,
                             const double ftol,
                             void (*callback)(const size_t iter, void *params,
                                              const gsl_multifit_nlinear_workspace *w),
                             void *callback_params,
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

int
gsl_multifit_nlinear_rcond (double *rcond, const gsl_multifit_nlinear_workspace * w);

int gsl_multifit_nlinear_eval_f(gsl_multifit_nlinear_fdf *fdf,
                                const gsl_vector *x,
                                const gsl_vector *swts,
                                gsl_vector *y);

int gsl_multifit_nlinear_eval_df(const gsl_vector *x,
                                 const gsl_vector *f,
                                 const gsl_vector *swts,
                                 const double h,
                                 const gsl_multifit_nlinear_fdtype fdtype,
                                 gsl_multifit_nlinear_fdf *fdf,
                                 gsl_matrix *df, gsl_vector *work);

int
gsl_multifit_nlinear_eval_fvv(const double h,
                              const gsl_vector *x,
                              const gsl_vector *v,
                              const gsl_vector *f,
                              const gsl_matrix *J,
                              const gsl_vector *swts,
                              gsl_multifit_nlinear_fdf *fdf,
                              gsl_vector *yvv, gsl_vector *work);

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
gsl_multifit_nlinear_df(const double h, const gsl_multifit_nlinear_fdtype fdtype,
                        const gsl_vector *x, const gsl_vector *wts,
                        gsl_multifit_nlinear_fdf *fdf,
                        const gsl_vector *f, gsl_matrix *J, gsl_vector *work);

/* fdfvv.c */
int
gsl_multifit_nlinear_fdfvv(const double h, const gsl_vector *x, const gsl_vector *v,
                           const gsl_vector *f, const gsl_matrix *J,
                           const gsl_vector *swts,
                           gsl_multifit_nlinear_fdf *fdf,
                           gsl_vector *fvv, gsl_vector *work);

/* top-level methods */
GSL_VAR const gsl_multifit_nlinear_type * gsl_multifit_nlinear_lm;

/* parameter update methods */
GSL_VAR const gsl_multifit_nlinear_update * gsl_multifit_nlinear_update_trust;
GSL_VAR const gsl_multifit_nlinear_update * gsl_multifit_nlinear_update_nielsen;

/* scaling matrix strategies */
GSL_VAR const gsl_multifit_nlinear_scale * gsl_multifit_nlinear_scale_levenberg;
GSL_VAR const gsl_multifit_nlinear_scale * gsl_multifit_nlinear_scale_marquardt;
GSL_VAR const gsl_multifit_nlinear_scale * gsl_multifit_nlinear_scale_more;

/* linear solvers */
GSL_VAR const gsl_multifit_nlinear_solver * gsl_multifit_nlinear_solver_normal;
GSL_VAR const gsl_multifit_nlinear_solver * gsl_multifit_nlinear_solver_qr;

__END_DECLS

#endif /* __GSL_MULTIFIT_NLINEAR_H__ */
