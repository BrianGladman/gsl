#ifndef GSL_MIN_H
#define GSL_MIN_H

#include <stdlib.h>
#include <gsl_math.h>

typedef struct
  {
    const char *name;
    size_t size;
    int (*set) (void *state, gsl_function * f, double * minimum, gsl_interval * x);
    int (*iterate) (void *state, gsl_function * f, double * minimum, gsl_interval * x);
  }
gsl_min_fsolver_type;

typedef struct
  {
    const gsl_min_fsolver_type * type;
    gsl_function * function ;
    double mininum ;
    gsl_interval interval ;
    void *state;
  }
gsl_min_fsolver;

gsl_min_fsolver *
gsl_min_fsolver_alloc (const gsl_min_fsolver_type * T, 
			 gsl_function * f, gsl_interval x);
void gsl_min_fsolver_free (gsl_min_fsolver * s);

int gsl_min_fsolver_set (gsl_min_fsolver * s, 
			   gsl_function * f, gsl_interval x);

int gsl_min_fsolver_iterate (gsl_min_fsolver * s);

const char * gsl_min_fsolver_name (const gsl_min_fsolver * s);
double gsl_min_fsolver_min (const gsl_min_fsolver * s);
gsl_interval gsl_min_fsolver_interval (const gsl_min_fsolver * s);


int
gsl_min_test_interval (gsl_interval x, double epsabs, double epsrel);

int
gsl_min_test_residual (double f, double epsabs);

int
gsl_min_test_delta (double x1, double x0, double epsabs, double epsrel);

extern const gsl_min_fsolver_type  * gsl_min_fsolver_bisection;
extern const gsl_min_fsolver_type  * gsl_min_fsolver_brent;
extern const gsl_min_fsolver_type  * gsl_min_fsolver_falsepos;

#endif /* GSL_MIN_H */
