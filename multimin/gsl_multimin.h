/* multimin/gsl_multimin.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Fabrice Rossi
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

#ifndef __GSL_MULTIMIN_H__
#define __GSL_MULTIMIN_H__

#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_min.h>

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

/* Definition of an arbitrary real-valued function with gsl_vector input and */
/* parameters */
struct gsl_multimin_function_struct 
{
  double (* f) (const gsl_vector * x, void * params);
  size_t n;
  void * params;
};

typedef struct gsl_multimin_function_struct gsl_multimin_function;

#define GSL_MULTIMIN_FN_EVAL(F,x) (*((F)->f))(x,(F)->params)

/* Definition of an arbitrary differentiable real-valued function */
/* with gsl_vector input and parameters */
struct gsl_multimin_function_fdf_struct 
{
  double (* f) (const gsl_vector  * x, void * params);
  void (* df) (const gsl_vector * x, void * params,gsl_vector * df);
  void (* fdf) (const gsl_vector * x, void * params,double *f,gsl_vector * df);
  size_t n;
  void * params;
};

typedef struct gsl_multimin_function_fdf_struct gsl_multimin_function_fdf;

#define GSL_MULTIMIN_FN_EVAL_F(F,x) (*((F)->f))(x,(F)->params)
#define GSL_MULTIMIN_FN_EVAL_DF(F,x,g) (*((F)->df))(x,(F)->params,(g))
#define GSL_MULTIMIN_FN_EVAL_F_DF(F,x,y,g) (*((F)->fdf))(x,(F)->params,(y),(g))

int gsl_diff_gradient(const gsl_multimin_function *f,
		      gsl_vector *x,
		      gsl_vector *g);

/* wrapper for using a vector input function as a real input function */
typedef struct 
{
  gsl_multimin_function f;
  gsl_vector * starting_point;
  gsl_vector * direction;
  gsl_vector * evaluation_point;
}
gsl_multimin_to_single;

gsl_multimin_to_single *
gsl_multimin_to_single_alloc(const gsl_multimin_function *f,
			     gsl_vector * starting_point,
			     gsl_vector * direction);
int
gsl_multimin_to_single_set(gsl_multimin_to_single *w,
			   const gsl_multimin_function *f,
			   gsl_vector * starting_point,
			   gsl_vector * direction);

void
gsl_multimin_to_single_free(gsl_multimin_to_single *w);

gsl_function *
gsl_multimin_to_single_function_alloc(gsl_multimin_to_single *w);

void
gsl_multimin_to_single_function_free(gsl_function *f);

gsl_multimin_to_single *
gsl_multimin_to_single_alloc_fdf(const gsl_multimin_function_fdf *f,
				 gsl_vector * starting_point,
				 gsl_vector * direction);
int
gsl_multimin_to_single_set_fdf(gsl_multimin_to_single *w,
			       const gsl_multimin_function_fdf *f,
			       gsl_vector * starting_point,
			       gsl_vector * direction);

void
gsl_multimin_compute_ep(gsl_multimin_to_single *params,double x);

void
gsl_multimin_compute_evaluation_point(gsl_vector *evaluation_point,
				      const gsl_vector *starting_point,
				      double x,
				      const gsl_vector *direction);

/* minimisation of differentiable functions */

typedef struct 
{
  double f; /* function value */
  double f1; /* previous function value */
  gsl_vector *x; /* minimum estimate */
  gsl_vector *x1; /* previous minimum estimate */
  gsl_vector *g; /* gradient */
  gsl_vector *g1; /* previous gradient */
}
gsl_multimin_fdf_history;

gsl_multimin_fdf_history *
gsl_multimin_fdf_history_alloc(gsl_multimin_function_fdf *fdf,
			       const gsl_vector * x);

int
gsl_multimin_fdf_history_set(gsl_multimin_fdf_history *h,
			     gsl_multimin_function_fdf *fdf,
			     const gsl_vector * x);

void
gsl_multimin_fdf_history_free(gsl_multimin_fdf_history *h);

int
gsl_multimin_fdf_history_step(gsl_multimin_fdf_history *h,
			      gsl_multimin_function_fdf *fdf,
			      const gsl_vector * direction,
			      double step);

int
gsl_multimin_fdf_history_step_with_value(gsl_multimin_fdf_history *h,
					 gsl_multimin_function_fdf *fdf,
					 const gsl_vector * direction,
					 double step,double fx);

typedef struct 
{
  const char *name;
  size_t size;
  int (*alloc) (void *state, size_t n);
  int (*restart) (void *state);
  int (*direction) (void *state,gsl_multimin_fdf_history *h,gsl_vector * dir);
  void (*free) (void *state);
}
gsl_multimin_fdf_minimizer_type;

typedef struct 
{
  /* multi dimensional part */
  const gsl_multimin_fdf_minimizer_type *type;
  gsl_multimin_function_fdf *fdf;
  gsl_multimin_fdf_history *history;
  /* one dimensional part */
  gsl_min_bracketing_function bracketing;
  gsl_function *f_directional;
  const gsl_min_fminimizer_type * line_search_type;
  gsl_min_fminimizer *line_search;
  void *state;
}
gsl_multimin_fdf_minimizer;

gsl_multimin_fdf_minimizer *
gsl_multimin_fdf_minimizer_alloc(const gsl_multimin_fdf_minimizer_type *T,
				 gsl_multimin_function_fdf *fdf,
				 const gsl_vector * x,
				 gsl_min_bracketing_function bracket,
				 const gsl_min_fminimizer_type * T_line);
void
gsl_multimin_fdf_minimizer_free(gsl_multimin_fdf_minimizer *s);

const gsl_vector *
gsl_multimin_fdf_minimizer_direction(gsl_multimin_fdf_minimizer *s);

int
gsl_multimin_fdf_minimizer_next_direction(gsl_multimin_fdf_minimizer *s);

int
gsl_multimin_fdf_minimizer_bracket(gsl_multimin_fdf_minimizer *s,
				   double first_step,size_t eval_max);

int
gsl_multimin_fdf_minimizer_iterate(gsl_multimin_fdf_minimizer *s);

int
gsl_multimin_fdf_minimizer_step(gsl_multimin_fdf_minimizer *s,
				double step);
int
gsl_multimin_fdf_minimizer_best_step(gsl_multimin_fdf_minimizer *s);

int
gsl_multimin_fdf_minimizer_restart(gsl_multimin_fdf_minimizer *s);

int
gsl_multimin_fdf_minimizer_step_with_value(gsl_multimin_fdf_minimizer *s,
					   double step,double f_at_end);

extern const gsl_multimin_fdf_minimizer_type *gsl_multimin_fdf_minimizer_steepest_descent;

extern const gsl_multimin_fdf_minimizer_type *gsl_multimin_fdf_minimizer_conjugate_pr;

extern const gsl_multimin_fdf_minimizer_type *gsl_multimin_fdf_minimizer_conjugate_fr;

extern const gsl_multimin_fdf_minimizer_type *gsl_multimin_fdf_minimizer_vector_bfgs;

int
gsl_multimin_test_gradient_sqr_norm(gsl_multimin_fdf_history *h,double epsabs);

__END_DECLS

#endif /* __GSL_MULTIMIN_H__ */
