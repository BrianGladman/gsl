#ifndef GSL_MULTIROOTS_H
#define GSL_MULTIROOTS_H

#include <stdlib.h>
#include <gsl_math.h>

typedef struct
  {
    const char *name;
    size_t size;
    int (*set) (void *state, gsl_function * f, double * root, gsl_interval * x);
    int (*iterate) (void *state, gsl_function * f, double * root, gsl_interval * x);
  }
gsl_multiroot_fsolver_type;

typedef struct
  {
    const gsl_multiroot_fsolver_type * type;
    gsl_function * function ;
    double root ;
    gsl_interval interval ;
    void *state;
  }
gsl_multiroot_fsolver;

typedef struct
  {
    const char *name;
    size_t size;
    int (*set) (void *state, gsl_function_fdf * f, double * root);
    int (*iterate) (void *state, gsl_function_fdf * f, double * root);
  }
gsl_multiroot_fdfsolver_type;

typedef struct
  {
    const gsl_multiroot_fdfsolver_type * type;
    gsl_function_fdf * fdf ;
    double root ;
    void *state;
  }
gsl_multiroot_fdfsolver;

gsl_multiroot_fsolver *
gsl_multiroot_fsolver_alloc (const gsl_multiroot_fsolver_type * T, 
			 gsl_function * f, gsl_interval x);
void gsl_multiroot_fsolver_free (gsl_multiroot_fsolver * s);

int gsl_multiroot_fsolver_set (gsl_multiroot_fsolver * s, 
			   gsl_function * f, gsl_interval x);

int gsl_multiroot_fsolver_iterate (gsl_multiroot_fsolver * s);

const char * gsl_multiroot_fsolver_name (const gsl_multiroot_fsolver * s);
double gsl_multiroot_fsolver_root (const gsl_multiroot_fsolver * s);
gsl_interval gsl_multiroot_fsolver_interval (const gsl_multiroot_fsolver * s);


gsl_multiroot_fdfsolver *
gsl_multiroot_fdfsolver_alloc (const gsl_multiroot_fdfsolver_type * T, 
			   gsl_function_fdf * fdf, double root);

int
gsl_multiroot_fdfsolver_set (gsl_multiroot_fdfsolver * s, 
			 gsl_function_fdf * fdf, double root);

int
gsl_multiroot_fdfsolver_iterate (gsl_multiroot_fdfsolver * s);

void
gsl_multiroot_fdfsolver_free (gsl_multiroot_fdfsolver * s);

const char * gsl_multiroot_fdfsolver_name (const gsl_multiroot_fdfsolver * s);
double gsl_multiroot_fdfsolver_root (const gsl_multiroot_fdfsolver * s);

int
gsl_multiroot_test_interval (gsl_interval x, double epsabs, double epsrel);

int
gsl_multiroot_test_residual (double f, double epsabs);

int
gsl_multiroot_test_delta (double x1, double x0, double epsabs, double epsrel);

extern const gsl_multiroot_fsolver_type  * gsl_multiroot_fsolver_bisection;
extern const gsl_multiroot_fdfsolver_type  * gsl_multiroot_fdfsolver_newton;


#endif /* GSL_MULTIROOTS_H */
