#ifndef GSL_ROOTS_H
#define GSL_ROOTS_H

#include <stdlib.h>
#include <gsl_math.h>

typedef struct
  {
    const char *name;
    size_t size;
    int (*set) (void *state, gsl_function * f, double * root, gsl_interval * x);
    int (*iterate) (void *state, gsl_function * f, double * root, gsl_interval * x);
  }
gsl_root_f_solver_type;

typedef struct
  {
    const char *name;
    size_t size;
    int (*set) (void *state, gsl_function * f, double * root, gsl_interval * x);
    int (*iterate) (void *state, gsl_function * f, double * root, gsl_interval * x);
    gsl_function * function ;
    double root ;
    gsl_interval interval ;
    void *state;
  }
gsl_root_f_solver;

typedef struct
  {
    const char *name;
    size_t size;
    int (*set) (void *state, gsl_fdf * f, double * root);
    int (*iterate) (void *state, gsl_fdf * f, double * root);
  }
gsl_root_fdf_solver_type;

typedef struct
  {
    const char *name;
    size_t size;
    int (*set) (void *state, gsl_fdf * f, double * root);
    int (*iterate) (void *state, gsl_fdf * f, double * root);
    gsl_fdf * fdf ;
    double root ;
    void *state;
  }
gsl_root_fdf_solver;

gsl_root_f_solver *
gsl_root_f_solver_alloc (const gsl_root_f_solver_type * T, 
			 gsl_function * f, gsl_interval x);

int
gsl_root_f_solver_set (gsl_root_f_solver * s, 
		       gsl_function * f, gsl_interval x);

int
gsl_root_f_solver_iterate (gsl_root_f_solver * s);

void
gsl_root_f_solver_free (gsl_root_f_solver * s);


gsl_root_fdf_solver *
gsl_root_fdf_solver_alloc (const gsl_root_fdf_solver_type * T, 
			   gsl_fdf * fdf, double root);

int
gsl_root_fdf_solver_set (gsl_root_fdf_solver * s, 
			 gsl_fdf * fdf, double root);

int
gsl_root_fdf_solver_iterate (gsl_root_fdf_solver * s);

void
gsl_root_fdf_solver_free (gsl_root_fdf_solver * s);


int
gsl_root_test_interval (gsl_interval x, double rel_epsilon, double abs_epsilon);

int
gsl_root_test_residual (double f, double abs_epsilon);

extern const gsl_root_f_solver_type  * gsl_root_f_solver_bisection;
extern const gsl_root_f_solver_type  * gsl_root_f_solver_brent;
extern const gsl_root_f_solver_type  * gsl_root_f_solver_falsepos;
extern const gsl_root_fdf_solver_type  * gsl_root_fdf_solver_newton;
extern const gsl_root_fdf_solver_type  * gsl_root_fdf_solver_secant;
extern const gsl_root_fdf_solver_type  * gsl_root_fdf_solver_steffenson;

/* Requested epsilon must be greater than GSL_DBL_EPSILON by this
   factor to protect against roundoff problems. */

#define GSL_ROOT_EPSILON_BUFFER 10.0

#endif /* GSL_ROOTS_H */
