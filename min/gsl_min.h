#ifndef __GSL_MIN_H__
#define __GSL_MIN_H__

#include <stdlib.h>
#include <gsl/gsl_math.h>

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
    const char *name;
    size_t size;
    int (*set) (void *state, gsl_function * f, double minimum, double f_minimum, gsl_interval x, double f_lower, double f_upper);
    int (*iterate) (void *state, gsl_function * f, double * minimum, double * f_minimum, gsl_interval * x, double * f_lower, double * f_upper);
  }
gsl_min_fminimizer_type;

typedef struct
  {
    const gsl_min_fminimizer_type * type;
    gsl_function * function ;
    double minimum ;
    gsl_interval interval ;
    double f_minimum, f_lower, f_upper;
    void *state;
  }
gsl_min_fminimizer;

gsl_min_fminimizer *
gsl_min_fminimizer_alloc (const gsl_min_fminimizer_type * T, 
			 gsl_function * f, double minimum, gsl_interval x);

gsl_min_fminimizer *
gsl_min_fminimizer_alloc_with_values (const gsl_min_fminimizer_type * T, 
                                      gsl_function * f, 
                                      double minimum, double f_minimum,
                                      gsl_interval x, 
                                      double f_lower, double f_upper);

void gsl_min_fminimizer_free (gsl_min_fminimizer * s);

int gsl_min_fminimizer_set (gsl_min_fminimizer * s, 
			   gsl_function * f, double minimum, gsl_interval x);

int gsl_min_fminimizer_set_with_values (gsl_min_fminimizer * s, 
                                        gsl_function * f, 
                                        double minimum, double f_minimum,
                                        gsl_interval x,
                                        double f_lower, double f_upper);

int gsl_min_fminimizer_iterate (gsl_min_fminimizer * s);

const char * gsl_min_fminimizer_name (const gsl_min_fminimizer * s);
double gsl_min_fminimizer_minimum (const gsl_min_fminimizer * s);
gsl_interval gsl_min_fminimizer_interval (const gsl_min_fminimizer * s);

int
gsl_min_test_interval (const gsl_interval x, double epsabs, double epsrel);

extern const gsl_min_fminimizer_type  * gsl_min_fminimizer_goldensection;
extern const gsl_min_fminimizer_type  * gsl_min_fminimizer_brent;

typedef
int (*gsl_min_bracketing_function)(gsl_function *f,
				   double *minimum,double * f_minimum,
				   gsl_interval *x, double * f_lower,
				   double * f_upper,
				   size_t eval_max);

int 
gsl_min_find_bracket(gsl_function *f,double *minimum,double * f_minimum,
		     gsl_interval *x, double * f_lower,double * f_upper,
		     size_t eval_max);

__END_DECLS

#endif /* __GSL_MIN_H__ */
