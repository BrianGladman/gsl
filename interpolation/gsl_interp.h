/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef __GSL_INTERP_H__
#define __GSL_INTERP_H__
#include <stdlib.h>

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
  int     (*eval_i_impl) (const struct _gsl_interp_obj_struct *, const double xa[], const double ya[], gsl_interp_accel *, double a, double b, double * result);
  void    (*free)        (struct _gsl_interp_obj_struct *);
  double  xmin;
  double  xmax;
  size_t  size;
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

int
gsl_interp_eval_integ_impl(const gsl_interp_obj * obj,
                           const double xa[], const double ya[],
                           double a, double b,
			   gsl_interp_accel * acc,
                           double * y
                           );

int
gsl_interp_eval_integ_e(const gsl_interp_obj * obj,
                        const double xa[], const double ya[],
                        double a, double b,
			gsl_interp_accel * acc,
                        double * y
                        );

double
gsl_interp_eval_integ(const gsl_interp_obj * obj,
                      const double xa[], const double ya[],
                      double a, double b,
		      gsl_interp_accel * acc
                      );

void
gsl_interp_obj_free(gsl_interp_obj * interp_obj);

size_t gsl_interp_bsearch(const double x_array[], double x,
                          size_t index_lo, size_t index_hi );

#ifdef HAVE_INLINE
extern inline size_t
gsl_interp_bsearch(const double x_array[], double x,
               size_t index_lo, size_t index_hi );

extern inline size_t
gsl_interp_bsearch(const double x_array[], double x,
               size_t index_lo, size_t index_hi )
{
  size_t ilo = index_lo;
  size_t ihi = index_hi;
  while(ihi > ilo + 1) {
    size_t i = (ihi + ilo)/2;
    if(x_array[i] > x)
      ihi = i;
    else
      ilo = i;
  }
  
  return ilo;
}
#endif

size_t
gsl_interp_accel_find(gsl_interp_accel * a, const double xa[], size_t len, double x);

#ifdef HAVE_INLINE
extern inline size_t
gsl_interp_accel_find(gsl_interp_accel * a, const double xa[], size_t len, double x)
{
  size_t x_index = a->cache;
 
  if(x < xa[x_index]) {
    a->miss_count++;
    a->cache = gsl_interp_bsearch(xa, x, 0, x_index);
  }
  else if(x > xa[x_index + 1]) {
    a->miss_count++;
    a->cache = gsl_interp_bsearch(xa, x, x_index, len-1);
  }
  else {
    a->hit_count++;
  }
  
  return a->cache;
}
#endif /* HAVE_INLINE */


#endif /* __GSL_INTERP_H__ */
