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
gsl_min_fminimizer_type;

typedef struct
  {
    const gsl_min_fminimizer_type * type;
    gsl_function * function ;
    double minimum ;
    gsl_interval interval ;
    void *state;
  }
gsl_min_fminimizer;

gsl_min_fminimizer *
gsl_min_fminimizer_alloc (const gsl_min_fminimizer_type * T, 
			 gsl_function * f, double minimum, gsl_interval x);
void gsl_min_fminimizer_free (gsl_min_fminimizer * s);

int gsl_min_fminimizer_set (gsl_min_fminimizer * s, 
			   gsl_function * f, double minimum, gsl_interval x);

int gsl_min_fminimizer_iterate (gsl_min_fminimizer * s);

const char * gsl_min_fminimizer_name (const gsl_min_fminimizer * s);
double gsl_min_fminimizer_minimum (const gsl_min_fminimizer * s);
gsl_interval gsl_min_fminimizer_interval (const gsl_min_fminimizer * s);

int
gsl_min_test_interval (gsl_interval x, double epsabs, double epsrel);

extern const gsl_min_fminimizer_type  * gsl_min_fminimizer_goldensection;
extern const gsl_min_fminimizer_type  * gsl_min_fminimizer_brent;

#endif /* GSL_MIN_H */
