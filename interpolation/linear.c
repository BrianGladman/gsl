/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include "gsl_interp.h"


/* linear interpolation object */
typedef struct
  {
    int (*eval_impl) (const gsl_interp_obj *, const double xa[], const double ya[], double x, gsl_interp_accel *, double *y);
    int (*eval_d_impl) (const gsl_interp_obj *, const double xa[], const double ya[], double x, gsl_interp_accel *, double *dydx);
    int (*eval_i_impl) (const struct _gsl_interp_obj_struct *, const double xa[], const double ya[], gsl_interp_accel *, double a, double b, double * result);    
    void (*free) (gsl_interp_obj *);
    double xmin;
    double xmax;
    size_t size;
  }
gsl_interp_linear;


static
gsl_interp_obj *
linear_create (const double xa[], const double ya[], size_t size);

static
void
linear_free (gsl_interp_obj * interp);

static
int
linear_eval_impl (const gsl_interp_obj *, const double xa[], const double ya[], double x, gsl_interp_accel *, double *y);

static
int
linear_eval_d_impl (const gsl_interp_obj *, const double xa[], const double ya[], double x, gsl_interp_accel *, double *dydx);

static
int
linear_eval_i_impl (const gsl_interp_obj *, const double xa[], const double ya[], gsl_interp_accel *, double, double, double *);


const gsl_interp_factory gsl_interp_factory_linear =
{
  "linear",
  linear_create
};



static
gsl_interp_obj *
linear_create (const double x_array[], const double y_array[], size_t size)
{
  y_array = 0;			/* prevent warning about unused parameter */

  if (size <= 1)
    return 0;
  else
    {
      gsl_interp_linear *interp = (gsl_interp_linear *) malloc (sizeof (gsl_interp_linear));
      if (interp != 0)
	{
	  interp->eval_impl = linear_eval_impl;
	  interp->eval_d_impl = linear_eval_d_impl;
	  interp->eval_i_impl = linear_eval_i_impl;
	  interp->free = linear_free;
	  interp->xmin = x_array[0];
	  interp->xmax = x_array[size - 1];
	  interp->size = size;
	}
      return (gsl_interp_obj *) interp;
    }
}


static
void
linear_free (gsl_interp_obj * linear_interp)
{
  if (linear_interp != 0)
    free ((gsl_interp_linear *) linear_interp);
}


static
int
linear_eval_impl (const gsl_interp_obj * linear_interp,
		  const double x_array[], const double y_array[],
		  double x,
		  gsl_interp_accel * a,
		  double *y)
{
  const gsl_interp_linear *interp = (const gsl_interp_linear *) linear_interp;

  if (x < interp->xmin)
    {
      *y = y_array[0];
      return GSL_EDOM;
    }
  else if (x > interp->xmax)
    {
      *y = y_array[interp->size - 1];
      return GSL_EDOM;
    }
  else
    {
      double x_lo, x_hi;
      double y_lo, y_hi;
      double dx;
      size_t index;

      if (a != 0)
	{
	  index = gsl_interp_accel_find (a, x_array, interp->size, x);
	}
      else
	{
	  index = gsl_interp_bsearch (x_array, x, 0, interp->size - 1);
	}

      /* evaluate */
      x_lo = x_array[index];
      x_hi = x_array[index + 1];
      y_lo = y_array[index];
      y_hi = y_array[index + 1];
      dx = x_hi - x_lo;
      if (dx > 0.0)
	{
	  *y = y_lo + (x - x_lo) / dx * (y_hi - y_lo);
	  return GSL_SUCCESS;
	}
      else
	{
	  *y = 0.0;
	  return GSL_EINVAL;
	}
    }
}


static
int
linear_eval_d_impl (const gsl_interp_obj * linear_interp,
		    const double x_array[], const double y_array[],
		    double x,
		    gsl_interp_accel * a,
		    double *dydx)
{
  const gsl_interp_linear *interp = (const gsl_interp_linear *) linear_interp;

  if (x < interp->xmin)
    {
      *dydx = 0.0;
      return GSL_EDOM;
    }
  else if (x > interp->xmax)
    {
      *dydx = 0.0;
      return GSL_EDOM;
    }
  else
    {
      double x_lo, x_hi;
      double y_lo, y_hi;
      double dx;
      double dy;
      size_t index;

      if (a != 0)
	{
	  index = gsl_interp_accel_find (a, x_array, interp->size, x);
	}
      else
	{
	  index = gsl_interp_bsearch (x_array, x, 0, interp->size - 1);
	}

      /* evaluate */
      x_lo = x_array[index];
      x_hi = x_array[index + 1];
      y_lo = y_array[index];
      y_hi = y_array[index + 1];
      dx = x_hi - x_lo;
      dy = y_hi - y_lo;
      if (dx > 0.0)
	{
	  *dydx = dy / dx;;
	  return GSL_SUCCESS;
	}
      else
	{
	  *dydx = 0.0;
	  return GSL_EINVAL;
	}
    }
}


static
int
linear_eval_i_impl (const gsl_interp_obj * linear_interp,
		    const double x_array[], const double y_array[],
		    gsl_interp_accel * acc,
                    double a, double b,
		    double * result)
{
  const gsl_interp_linear *interp = (const gsl_interp_linear *) linear_interp;

  if (a > b || a < interp->xmin || b > interp->xmax)
    {
      *result = 0.0;
      return GSL_EDOM;
    }
  else if(a == b)
    {
      *result = 0.0;
      return GSL_SUCCESS;
    }
  else
    {
      size_t index_a, index_b;

      if (acc != 0)
	{
	  index_a = gsl_interp_accel_find (acc, x_array, interp->size, a);
	  index_b = gsl_interp_accel_find (acc, x_array, interp->size, b);
	}
      else
	{
	  index_a = gsl_interp_bsearch (x_array, a, 0, interp->size - 1);
	  index_b = gsl_interp_bsearch (x_array, b, 0, interp->size - 1);
	}

      if(index_a == index_b) {
        /* endpoints inside same interval */
        const double x_hi = x_array[index_a + 1];
        const double x_lo = x_array[index_a];
	const double y_lo = y_array[index_a];
	const double y_hi = y_array[index_a + 1];
	const double dx = x_hi - x_lo;
	const double dy = y_hi - y_lo;
	if(dx != 0.0) {
	  const double D = dy/dx;
	  *result = (b-a) * (y_lo + 0.5*D*(b+a - 2.0*x_lo));
	  return GSL_SUCCESS;
	}
	else {
	  *result = 0.0;
	  return GSL_FAILURE;
	}
      }
      else {
        /* endpoints span more than one interval */
	size_t i;
	*result = 0.0;

	/* interior intervals */
	for(i=index_a+1; i<index_b; i++) {
	  const double x_hi = x_array[i + 1];
          const double x_lo = x_array[i];
	  const double y_lo = y_array[i];
	  const double y_hi = y_array[i + 1];
	  const double dx = x_hi - x_lo;
	  *result += 0.5 * dx * (y_lo + y_hi);
	}

        /* lower end interval */
	{
          const double x_hi = x_array[index_a + 1];
          const double x_lo = x_array[index_a];
	  const double y_lo = y_array[index_a];
	  const double y_hi = y_array[index_a + 1];
	  const double dx = x_hi - x_lo;
	  const double dy = y_hi - y_lo;
	  if(dx != 0.0) {
	    const double D = dy/dx;
	    *result += (x_hi-a) * (y_lo + 0.5*D*(x_hi-a));
          }
	  else {
            *result = 0.0;
            return GSL_FAILURE;
	  }
	}

        /* upper end interval */
	{
          const double x_hi = x_array[index_b + 1];
          const double x_lo = x_array[index_b];
	  const double y_lo = y_array[index_b];
	  const double y_hi = y_array[index_b + 1];
	  const double dx = x_hi - x_lo;
	  const double dy = y_hi - y_lo;
	  if(dx != 0.0) {
	    const double D = dy/dx;
	    *result += (b-x_lo) * (y_lo + 0.5*D*(b-x_lo));
          }
	  else {
            *result = 0.0;
            return GSL_FAILURE;
	  }
	}

        return GSL_SUCCESS;
      }
    }
}
