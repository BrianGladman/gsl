/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <stdlib.h>
#include <gsl_errno.h>
#include "bsearch.h"
#include "gsl_interp.h"


/* linear interpolation object */
typedef struct {
  int     (*eval_impl)   (const gsl_interp_obj *, const double xa[], const double ya[], double x, gsl_interp_accel *, double * y);
  int     (*eval_d_impl) (const gsl_interp_obj *, const double xa[], const double ya[], double x, gsl_interp_accel *, double * dydx);
  void    (*free)        (gsl_interp_obj *);
  double  xmin;
  double  xmax;
  int     size;
}
gsl_interp_linear;


static
gsl_interp_obj *
linear_create(const double xa[], const double ya[], int size);

static
void
linear_free(gsl_interp_obj * interp);

static
int
linear_eval_impl(const gsl_interp_obj *, const double xa[], const double ya[], double x, gsl_interp_accel *, double * y);

static
int
linear_eval_d_impl(const gsl_interp_obj *, const double xa[], const double ya[], double x, gsl_interp_accel *, double * dydx);



const gsl_interp_factory gsl_interp_factory_linear = {
  "linear",
  linear_create
};



static
gsl_interp_obj *
linear_create(const double x_array[], const double y_array[], int size)
{
  if(size <= 1)
    return 0;
  else {
    gsl_interp_linear * interp = (gsl_interp_linear *) malloc(sizeof(gsl_interp_linear));
    if(interp != 0) {
      interp->eval_impl   = linear_eval_impl;
      interp->eval_d_impl = linear_eval_d_impl;
      interp->free        = linear_free;
      interp->xmin = x_array[0];
      interp->xmax = x_array[size-1];
      interp->size = size;
    }
    return (gsl_interp_obj *) interp;
  }
}


static
void
linear_free(gsl_interp_obj * linear_interp)
{
  if(linear_interp != 0) free((gsl_interp_linear * ) linear_interp);
}


static
int
linear_eval_impl(const gsl_interp_obj * linear_interp,
                 const double x_array[], const double y_array[],
		 double x,
		 gsl_interp_accel * a,
		 double * y
		 )
{
  const gsl_interp_linear * interp = (const gsl_interp_linear * ) linear_interp;

  if(x < interp->xmin) {
    *y = y_array[0];
    return GSL_EDOM;
  }
  else if(x > interp->xmax) {
    *y = y_array[interp->size - 1];
    return GSL_EDOM;
  }
  else {
    double x_lo, x_hi;
    double y_lo, y_hi;
    double dx;
    unsigned long index;

    if(a != 0) {
      index = gsl_interp_accel_find(a, x_array, interp->size, x);
    }
    else {
      index = interp_bsearch(x_array, x, 0, interp->size - 1);
    }

    /* evaluate */
    x_lo = x_array[index];
    x_hi = x_array[index + 1];
    y_lo = y_array[index];
    y_hi = y_array[index + 1];
    dx = x_hi - x_lo;
    if(dx > 0.0) {
      *y = y_lo + (x - x_lo)/dx * (y_hi - y_lo);
      return GSL_SUCCESS;
    }
    else {
      *y = 0.0;
      return GSL_EINVAL;
    }
  }
}


static
int
linear_eval_d_impl(const gsl_interp_obj * linear_interp,
                   const double x_array[], const double y_array[],
                   double x,
		   gsl_interp_accel * a,
                   double * dydx
                   )
{
  const gsl_interp_linear * interp = (const gsl_interp_linear * ) linear_interp;

  if(x < interp->xmin) {
    *dydx = 0.0;
    return GSL_EDOM;
  }
  else if(x > interp->xmax) {
    *dydx = 0.0;
    return GSL_EDOM;
  }
  else {
    double x_lo, x_hi;
    double y_lo, y_hi;
    double dx;
    double dy;
    unsigned long index;

    if(a != 0) {
      index = gsl_interp_accel_find(a, x_array, interp->size, x);
    }
    else {
      index = interp_bsearch(x_array, x, 0, interp->size - 1);
    }

    /* evaluate */
    x_lo = x_array[index];
    x_hi = x_array[index + 1];
    y_lo = y_array[index];
    y_hi = y_array[index + 1];
    dx = x_hi - x_lo;
    dy = y_hi - y_lo;
    if(dx > 0.0) {
      *dydx = dy/dx;;
      return GSL_SUCCESS;
    }
    else {
      *dydx = 0.0;
      return GSL_EINVAL;
    }
  } 
}
