/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_INTERP_H_
#define GSL_INTERP_H_


/* evaluation accelerator */
typedef struct {
  size_t  cache;        /* cache of index   */
  size_t  miss_count;   /* keep statistics  */
  size_t  hit_count;
}
gsl_interp_accel;


/* general interpolation object */
struct _gsl_interp_obj_struct {
  int     (*eval_impl)   (const struct _gsl_interp_obj_struct *, const double xa[], const double ya[], double x, gsl_interp_accel *, double * y);
  int     (*eval_d_impl) (const struct _gsl_interp_obj_struct *, const double xa[], const double ya[], double x, gsl_interp_accel *, double * dydx);
  void    (*free)        (struct _gsl_interp_obj_struct *);
  double  xmin;
  double  xmax;
  size_t     size;
};
typedef  struct _gsl_interp_obj_struct  gsl_interp_obj;


/* interpolation object factory */
typedef struct {
  const char * name;
  gsl_interp_obj *  (*create) (const double x_array[], const double y_array[], size_t size);
}
gsl_interp_factory;


/* available factories */
extern const gsl_interp_factory   gsl_interp_factory_linear;
extern const gsl_interp_factory   gsl_interp_factory_cspline_natural;
extern const gsl_interp_factory   gsl_interp_factory_cspline_periodic;
extern const gsl_interp_factory   gsl_interp_factory_akima_natural;
extern const gsl_interp_factory   gsl_interp_factory_akima_periodic;


gsl_interp_accel *
gsl_interp_accel_new(void);

size_t
gsl_interp_accel_find(gsl_interp_accel * a, const double x_array[], size_t size, double x);

void
gsl_interp_accel_free(gsl_interp_accel * a);


int
gsl_interp_eval_impl(const gsl_interp_obj * obj,
                     const double xa[], const double ya[], double x,
                     gsl_interp_accel * a, double * y
                     );

int
gsl_interp_eval_e(const gsl_interp_obj * obj,
                  const double xa[], const double ya[], double x,
                  gsl_interp_accel * a, double * y
                  );

double
gsl_interp_eval(const gsl_interp_obj * obj,
                const double xa[], const double ya[], double x,
                gsl_interp_accel * a
                );

int
gsl_interp_eval_deriv_impl(const gsl_interp_obj * obj,
                           const double xa[], const double ya[], double x,
			   gsl_interp_accel * a,
                           double * y
                           );

int
gsl_interp_eval_deriv_e(const gsl_interp_obj * obj,
                        const double xa[], const double ya[], double x,
			gsl_interp_accel * a,
                        double * y
                        );

double
gsl_interp_eval_deriv(const gsl_interp_obj * obj,
                      const double xa[], const double ya[], double x,
		      gsl_interp_accel * a
                      );

void
gsl_interp_obj_free(gsl_interp_obj * interp_obj);


#ifdef HAVE_INLINE
#include "bsearch.h"
extern
inline
size_t
gsl_interp_accel_find(gsl_interp_accel * a, const double xa[], size_t len, double x)
{
  size_t x_index = a->cache;
 
  if(x < xa[x_index]) {
    a->miss_count++;
    a->cache = interp_bsearch(xa, x, 0, x_index);
  }
  else if(x > xa[x_index + 1]) {
    a->miss_count++;
    a->cache = interp_bsearch(xa, x, x_index, len-1);
  }
  else {
    a->hit_count++;
  }
  
  return a->cache;
}
#endif /* HAVE_INLINE */


#endif  /* !GSL_INTERP_H_ */
