/* interpolation/gsl_interp.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
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

/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef __GSL_INTERP_H__
#define __GSL_INTERP_H__
#include <stdlib.h>

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

/* evaluation accelerator */
typedef struct {
  size_t  cache;        /* cache of index   */
  size_t  miss_count;   /* keep statistics  */
  size_t  hit_count;
}
gsl_interp_accel;


/* general interpolation object */
struct _gsl_interp_obj_struct {
  int     (*eval)    (const struct _gsl_interp_obj_struct *, const double xa[], const double ya[], double x, gsl_interp_accel *, double * y);
  int     (*eval_d)  (const struct _gsl_interp_obj_struct *, const double xa[], const double ya[], double x, gsl_interp_accel *, double * y_p);
  int     (*eval_d2) (const struct _gsl_interp_obj_struct *, const double xa[], const double ya[], double x, gsl_interp_accel *, double * y_pp);
  int     (*eval_i)  (const struct _gsl_interp_obj_struct *, const double xa[], const double ya[], gsl_interp_accel *, double a, double b, double * result);
  void    (*free)         (struct _gsl_interp_obj_struct *);
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
gsl_interp_eval_deriv2_e(const gsl_interp_obj * obj,
                         const double xa[], const double ya[], double x,
			 gsl_interp_accel * a,
                         double * y
                         );

double
gsl_interp_eval_deriv2(const gsl_interp_obj * obj,
                       const double xa[], const double ya[], double x,
		       gsl_interp_accel * a
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


__END_DECLS

#endif /* __GSL_INTERP_H__ */
