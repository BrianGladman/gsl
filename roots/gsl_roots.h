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
gsl_root_fsolver_type;

typedef struct
  {
    const gsl_root_fsolver_type * type;
    gsl_function * function ;
    double root ;
    gsl_interval interval ;
    void *state;
  }
gsl_root_fsolver;

typedef struct
  {
    const char *name;
    size_t size;
    int (*set) (void *state, gsl_function_fdf * f, double * root);
    int (*iterate) (void *state, gsl_function_fdf * f, double * root);
  }
gsl_root_fdfsolver_type;

typedef struct
  {
    const gsl_root_fdfsolver_type * type;
    gsl_function_fdf * fdf ;
    double root ;
    void *state;
  }
gsl_root_fdfsolver;

gsl_root_fsolver *
gsl_root_fsolver_alloc (const gsl_root_fsolver_type * T, 
			 gsl_function * f, gsl_interval x);
void gsl_root_fsolver_free (gsl_root_fsolver * s);

int gsl_root_fsolver_set (gsl_root_fsolver * s, 
			   gsl_function * f, gsl_interval x);

int gsl_root_fsolver_iterate (gsl_root_fsolver * s);

const char * gsl_root_fsolver_name (const gsl_root_fsolver * s);
double gsl_root_fsolver_root (const gsl_root_fsolver * s);
gsl_interval gsl_root_fsolver_interval (const gsl_root_fsolver * s);


gsl_root_fdfsolver *
gsl_root_fdfsolver_alloc (const gsl_root_fdfsolver_type * T, 
			   gsl_function_fdf * fdf, double root);

int
gsl_root_fdfsolver_set (gsl_root_fdfsolver * s, 
			 gsl_function_fdf * fdf, double root);

int
gsl_root_fdfsolver_iterate (gsl_root_fdfsolver * s);

void
gsl_root_fdfsolver_free (gsl_root_fdfsolver * s);

const char * gsl_root_fdfsolver_name (const gsl_root_fdfsolver * s);
double gsl_root_fdfsolver_root (const gsl_root_fdfsolver * s);

int
gsl_root_test_interval (gsl_interval x, double epsabs, double epsrel);

int
gsl_root_test_residual (double f, double epsabs);

int
gsl_root_test_delta (double x1, double x0, double epsabs, double epsrel);

extern const gsl_root_fsolver_type  * gsl_root_fsolver_bisection;
extern const gsl_root_fsolver_type  * gsl_root_fsolver_brent;
extern const gsl_root_fsolver_type  * gsl_root_fsolver_falsepos;
extern const gsl_root_fdfsolver_type  * gsl_root_fdfsolver_newton;
extern const gsl_root_fdfsolver_type  * gsl_root_fdfsolver_secant;
extern const gsl_root_fdfsolver_type  * gsl_root_fdfsolver_steffenson;

/* Requested epsilon must be greater than GSL_DBL_EPSILON by this
   factor to protect against roundoff problems. */

#define GSL_ROOT_EPSILON_BUFFER 10.0

#endif /* GSL_ROOTS_H */
