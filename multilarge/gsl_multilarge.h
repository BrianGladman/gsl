/* gsl_multilarge.h
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

#ifndef __GSL_MULTILARGE_H__
#define __GSL_MULTILARGE_H__

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_types.h>

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

/* iteration solver type */
typedef struct
{
  const char *name;
  void * (*alloc) (const size_t n, const size_t m);
  int (*reset) (void *);
  int (*accumulate) (const gsl_matrix * X, const gsl_vector * y,
                     void *);
  int (*solve) (const double lambda, gsl_vector * c,
                double * rnorm, double * snorm, void *);
  void (*free) (void *);
} gsl_multilarge_linear_type;

typedef struct
{
  const gsl_multilarge_linear_type * type;
  void * state;
} gsl_multilarge_linear_workspace;

/* available types */
GSL_VAR const gsl_multilarge_linear_type * gsl_multilarge_linear_normal;
GSL_VAR const gsl_multilarge_linear_type * gsl_multilarge_linear_tsqr;

/*
 * Prototypes
 */
gsl_multilarge_linear_workspace *
gsl_multilarge_linear_alloc(const gsl_multilarge_linear_type * T,
                            const size_t nmax, const size_t p);
void gsl_multilarge_linear_free(gsl_multilarge_linear_workspace * w);
const char *gsl_multilarge_name(const gsl_multilarge_linear_workspace * w);
int gsl_multilarge_linear_reset(gsl_multilarge_linear_workspace * w);
int gsl_multilarge_linear_accumulate(const gsl_matrix * X,
                                     const gsl_vector * y,
                                     gsl_multilarge_linear_workspace * w);
int gsl_multilarge_linear_solve(const double lambda, gsl_vector * c,
                                double * rnorm, double * snorm,
                                gsl_multilarge_linear_workspace * w);

__END_DECLS

#endif /* __GSL_MULTILARGE_H__ */
