/* gsl_multilarge_nlin.h
 * 
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

#ifndef __GSL_MULTILARGE_NLIN_H__
#define __GSL_MULTILARGE_NLIN_H__

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_types.h>
#include <gsl/gsl_multilarge.h>

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
  int (* fdf) (const gsl_vector * x, gsl_matrix * JTJ,
               gsl_vector * JTf, double * normf,
               void * params);
  size_t p;              /* number of independent variables */
  void * params;         /* user parameters */
  size_t nevalf;         /* number of evaluations of f */
  size_t nevaldf;        /* number of evaluations of Jacobian */
} gsl_multilarge_function_fdf;

#define GSL_MULTILARGE_EVAL_FDF(F, x, A, b, normf, status) \
       do { \
       status = (*((F)->fdf)) (x, A, b, normf, (F)->params); \
       ++(F)->nevalf; \
       if (evaldf) { ++(F)->nevaldf; } \
       } while (0)

typedef struct
{
  const char *name;
  void * (*alloc) (const size_t p);
  int (*init) (const gsl_vector * x,
               gsl_multilarge_function_fdf * fdf,
               void * vstate);
  int (*iterate) (gsl_vector * x, gsl_vector * dx,
                  gsl_multilarge_function_fdf * fdf,
                  void * vstate);
  gsl_vector * (*gradient) (void *vstate);
  double (*normf) (void *vstate);
  int (*rcond) (double * rcond, void *vstate);
  void (*free) (void * vstate);
} gsl_multilarge_nlinear_type;

typedef struct
{
  const gsl_multilarge_nlinear_type * type;
  gsl_multilarge_function_fdf * callback;
  gsl_vector * x;        /* parameter values x */
  gsl_vector * dx;       /* step dx */
  size_t p;              /* number of model parameters */
  size_t niter;          /* number of iterations performed */
  void *state;           /* solver workspace */
} gsl_multilarge_nlinear_workspace;

/* available solvers */
GSL_VAR const gsl_multilarge_nlinear_type * gsl_multilarge_nlinear_lm;
GSL_VAR const gsl_multilarge_nlinear_type * gsl_multilarge_nlinear_lms;

gsl_multilarge_nlinear_workspace *
gsl_multilarge_nlinear_alloc (const gsl_multilarge_nlinear_type * T, 
                              const size_t p);

void gsl_multilarge_nlinear_free (gsl_multilarge_nlinear_workspace * w);

const char *gsl_multilarge_nlinear_name (const gsl_multilarge_nlinear_workspace * w);

size_t gsl_multilarge_nlinear_niter (const gsl_multilarge_nlinear_workspace * w);

int
gsl_multilarge_nlinear_init (const gsl_vector * x, gsl_multilarge_function_fdf * fdf,
                             gsl_multilarge_nlinear_workspace * w);

int gsl_multilarge_nlinear_iterate (gsl_multilarge_nlinear_workspace * w);

double gsl_multilarge_nlinear_normf (const gsl_multilarge_nlinear_workspace * w);

int gsl_multilarge_nlinear_rcond (double * rcond, const gsl_multilarge_nlinear_workspace * w);

gsl_vector *gsl_multilarge_nlinear_position (const gsl_multilarge_nlinear_workspace * w);

int gsl_multilarge_nlinear_driver (const size_t maxiter,
                                   const double xtol,
                                   const double gtol,
                                   const double ftol,
                                   int *info,
                                   gsl_multilarge_nlinear_workspace *w);

int gsl_multilarge_nlinear_applyW(const gsl_vector * w, gsl_matrix * J, gsl_vector * f);

int gsl_multilarge_nlinear_test (const double xtol, const double gtol,
                                 const double ftol, int *info,
                                 const gsl_multilarge_nlinear_workspace * w);

__END_DECLS

#endif /* __GSL_MULTILARGE_NLIN_H__ */
