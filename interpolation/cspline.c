/* interpolation/cspline.c
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
#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include "integ_eval_macro.h"
#include "gsl_interp.h"


/* cubic spline interpolation object */
typedef struct
  {
    int (*eval_impl) (const gsl_interp_obj *, const double xa[], const double ya[], double x, gsl_interp_accel *, double *y);
    int (*eval_d_impl) (const gsl_interp_obj *, const double xa[], const double ya[], double x, gsl_interp_accel *, double *dydx);
    int (*eval_d2_impl) (const gsl_interp_obj *, const double xa[], const double ya[], double x, gsl_interp_accel *, double *dydx);
    int (*eval_i_impl) (const struct _gsl_interp_obj_struct *, const double xa[], const double ya[], gsl_interp_accel *, double a, double b, double * result);
    void (*free) (gsl_interp_obj *);
    double xmin;
    double xmax;
    size_t size;
    double *c;
  }
gsl_interp_cspline;

static
gsl_interp_obj *
cspline_natural_create (const double xa[], const double ya[], size_t size);

static
gsl_interp_obj *
cspline_periodic_create (const double xa[], const double ya[], size_t size);

static
void
cspline_free (gsl_interp_obj * interp);

static
int
cspline_eval_impl (const gsl_interp_obj *, const double xa[], const double ya[], double x, gsl_interp_accel *, double *y);

static
int
cspline_eval_d_impl (const gsl_interp_obj *, const double xa[], const double ya[], double x, gsl_interp_accel *, double *dydx);

static
int
cspline_eval_d2_impl (const gsl_interp_obj *, const double xa[], const double ya[], double x, gsl_interp_accel *, double *dydx);

static
int
cspline_eval_i_impl (const gsl_interp_obj *, const double xa[], const double ya[], gsl_interp_accel *, double a, double b, double * result);


/* global cubic spline factory objects */
const gsl_interp_factory gsl_interp_factory_cspline_natural =
{
  "cspline_natural",
  cspline_natural_create
};
const gsl_interp_factory gsl_interp_factory_cspline_periodic =
{
  "cspline_periodic",
  cspline_periodic_create
};



/* common initialization */
static
gsl_interp_cspline *
cspline_new (const double xa[], size_t size)
{
  gsl_interp_cspline *interp = (gsl_interp_cspline *) malloc (sizeof (gsl_interp_cspline));
  if (interp != 0)
    {
      interp->eval_impl = cspline_eval_impl;
      interp->eval_d_impl = cspline_eval_d_impl;
      interp->eval_i_impl = cspline_eval_i_impl;
      interp->eval_d2_impl = cspline_eval_d2_impl;
      interp->free = cspline_free;
      interp->xmin = xa[0];
      interp->xmax = xa[size - 1];
      interp->size = size;
      interp->c = (double *) malloc (size * sizeof (double));
      if (interp->c == 0)
	{
	  free (interp);
	  return 0;
	}
    }
  return interp;
}


/* natural spline calculation
 * see [Engeln-Mullges + Uhlig, p. 254]
 */
static
int
cspline_calc_natural (
  gsl_interp_cspline * interp,
  const double xa[],
  const double ya[])
{
  size_t i;
  size_t num_points = interp->size;
  size_t max_index = num_points - 1;	/* Engeln-Mullges + Uhlig "n" */
  size_t sys_size = max_index - 1;	/* linear system is sys_size x sys_size */
  int status;

  gsl_vector * g = gsl_vector_alloc(sys_size);
  gsl_vector * diag = gsl_vector_alloc(sys_size);
  gsl_vector * offdiag = gsl_vector_alloc(sys_size);

  gsl_vector solution_vec = { sys_size, 1, interp->c + 1, 0 };

  if (g == 0 || diag == 0 || offdiag == 0)
    {
      status = GSL_ENOMEM;
    }
  else
    {
      interp->c[0] = 0.0;
      interp->c[max_index] = 0.0;

      for (i = 0; i < sys_size; i++)
	{
	  const double h_i   = xa[i + 1] - xa[i];
	  const double h_ip1 = xa[i + 2] - xa[i + 1];
	  const double ydiff_i   = ya[i + 1] - ya[i];
	  const double ydiff_ip1 = ya[i + 2] - ya[i + 1];
          gsl_vector_set(offdiag, i, h_ip1);
          gsl_vector_set(diag, i, 2.0 * (h_ip1 + h_i));
          gsl_vector_set(g, i, 3.0 * (ydiff_ip1 / h_ip1  -  ydiff_i / h_i));
	}

      status = gsl_linalg_solve_symm_tridiag(diag, offdiag, g, &solution_vec);
    }

  if (g != 0) gsl_vector_free (g);
  if (diag != 0) gsl_vector_free (diag);
  if (offdiag != 0) gsl_vector_free (offdiag);

  return status;
}


/* periodic spline calculation
 * see [Engeln-Mullges + Uhlig, p. 256]
 */
static
int
cspline_calc_periodic (
  gsl_interp_cspline * interp,
  const double xa[],
  const double ya[])
{
  size_t i;
  size_t num_points = interp->size;
  size_t max_index = num_points - 1;  /* Engeln-Mullges + Uhlig "n" */
  size_t sys_size = max_index;        /* linear system is sys_size x sys_size */
  int status;

  if (sys_size == 2) {
    /* solve 2x2 system */

    const double h0 = xa[1] - xa[0];
    const double h1 = xa[2] - xa[1];
    const double h2 = xa[3] - xa[2];
    const double A = 2.0*(h0 + h1);
    const double B = h0 + h1;
    double g[2];
    double det;

    g[0] = 3.0 * ((ya[2] - ya[1]) / h1 - (ya[1] - ya[0]) / h0);
    g[1] = 3.0 * ((ya[1] - ya[2]) / h2 - (ya[2] - ya[1]) / h1);

    det = 3.0 * (h0 + h1) * (h0 + h1);
    interp->c[1] = ( A * g[0] - B * g[1])/det;
    interp->c[2] = (-B * g[0] + A * g[1])/det;
    interp->c[0] = interp->c[2];

    status = GSL_SUCCESS;
  }
  else {
    gsl_vector * g = gsl_vector_alloc(sys_size);
    gsl_vector * diag = gsl_vector_alloc(sys_size);
    gsl_vector * offdiag = gsl_vector_alloc(sys_size);
    gsl_vector solution_vec = { sys_size, 1, interp->c + 1, 0 };

    if (g == 0 || diag == 0 || offdiag == 0) {
      status = GSL_ENOMEM;
    }
    else {
      gsl_vector_set(offdiag, max_index, xa[1]-xa[0]);

      for (i = 0; i < sys_size; i++) {
    	const double h_i   = xa[i + 1] - xa[i];
        const double h_ip1 = xa[i + 2] - xa[i + 1];
        const double ydiff_i   = ya[i + 1] - ya[i];
        const double ydiff_ip1 = ya[(i + 2) % num_points] - ya[i + 1];
        gsl_vector_set(offdiag, i, h_ip1);
        gsl_vector_set(diag, i, 2.0 * (h_ip1 + h_i));
        gsl_vector_set(g, i, 3.0 * (ydiff_ip1 / h_ip1 - ydiff_i / h_i));
      }

      status = gsl_linalg_solve_symm_cyc_tridiag(diag, offdiag, g, &solution_vec);
      interp->c[0] = interp->c[max_index];
    }

    if (g != 0) gsl_vector_free (g);
    if (diag != 0) gsl_vector_free (diag);
    if (offdiag != 0) gsl_vector_free (offdiag);
  }

  return status;
}


/* factory method */
static
gsl_interp_obj *
cspline_natural_create (const double x_array[], const double y_array[], size_t size)
{
  if (size <= 1)
    return 0;
  else
    {
      gsl_interp_cspline *interp = cspline_new (x_array, size);
      if (interp != 0)
	{
	  int status = cspline_calc_natural (interp, x_array, y_array);
	  if(status != GSL_SUCCESS) interp->free((gsl_interp_obj *)interp);
	}
      return (gsl_interp_obj *) interp;
    }
}


/* factory method */
static
gsl_interp_obj *
cspline_periodic_create (const double x_array[], const double y_array[], size_t size)
{
  if (size <= 1)
    return 0;
  else
    {
      gsl_interp_cspline *interp = cspline_new (x_array, size);
      if (interp != 0)
	{
	  int status = cspline_calc_periodic (interp, x_array, y_array);
	  if(status != GSL_SUCCESS) interp->free((gsl_interp_obj *)interp);
	}
      return (gsl_interp_obj *) interp;
    }
}


static
void
cspline_free (gsl_interp_obj * interp_cspline)
{
  gsl_interp_cspline *interp = (gsl_interp_cspline *) interp_cspline;

  if (interp != 0)
    {
      if (interp->c != 0)
	free (interp->c);
      free (interp);
    }
}

/* disgusting macro for common coefficient determination
 */
#define COEFF_CALC(c_array, dy, dx, index,  b_i, c_i, d_i)   \
do {                                                         \
  const double _c_ip1 = c_array[(index) + 1];                \
  c_i = c_array[index];                                      \
  b_i = (dy) / (dx) - (dx) * (_c_ip1 + 2.0 * c_i) / 3.0;     \
  d_i = (_c_ip1 - c_i) / (3.0 * (dx));                       \
} while(0)


static
int
cspline_eval_impl (const gsl_interp_obj * cspline_interp,
		   const double x_array[], const double y_array[],
		   double x,
		   gsl_interp_accel * a,
		   double *y
)
{
  const gsl_interp_cspline *interp = (const gsl_interp_cspline *) cspline_interp;

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
      x_hi = x_array[index + 1];
      x_lo = x_array[index];
      dx = x_hi - x_lo;
      if (dx > 0.0)
	{
	  const double y_lo = y_array[index];
	  const double y_hi = y_array[index + 1];
	  const double dy = y_hi - y_lo;
	  double delx = x - x_lo;
	  double b_i, c_i, d_i; 
	  COEFF_CALC(interp->c, dy, dx, index,  b_i, c_i, d_i);
	  *y = y_lo + delx * (b_i + delx * (c_i + delx * d_i));
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
cspline_eval_d_impl (const gsl_interp_obj * cspline_interp,
		     const double x_array[], const double y_array[],
		     double x,
		     gsl_interp_accel * a,
		     double *dydx)
{
  const gsl_interp_cspline *interp = (const gsl_interp_cspline *) cspline_interp;

  if (x > interp->xmax)
    {
      *dydx = 0.0;
      return GSL_EDOM;
    }
  else if (x < interp->xmin)
    {
      *dydx = 0.0;
      return GSL_EDOM;
    }
  else
    {
      double x_lo, x_hi;
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
      x_hi = x_array[index + 1];
      x_lo = x_array[index];
      dx = x_hi - x_lo;
      if (dx > 0.0)
	{
	  const double y_lo = y_array[index];
	  const double y_hi = y_array[index + 1];
	  const double dy = y_hi - y_lo;
	  double delx = x - x_lo;
	  double b_i, c_i, d_i; 
	  COEFF_CALC(interp->c, dy, dx, index,  b_i, c_i, d_i);
	  *dydx = b_i + delx * (2.0 * c_i + 3.0 * d_i * delx);
	  return GSL_SUCCESS;
	}
      else
	{
	  *dydx = 0.0;
	  return GSL_FAILURE;
	}
    }
}


static
int
cspline_eval_d2_impl (const gsl_interp_obj * cspline_interp,
		      const double x_array[], const double y_array[],
		      double x,
		      gsl_interp_accel * a,
		      double * y_pp)
{
  const gsl_interp_cspline *interp = (const gsl_interp_cspline *) cspline_interp;

  if (x > interp->xmax)
    {
      *y_pp = 0.0;
      return GSL_EDOM;
    }
  else if (x < interp->xmin)
    {
      *y_pp = 0.0;
      return GSL_EDOM;
    }
  else
    {
      double x_lo, x_hi;
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
      x_hi = x_array[index + 1];
      x_lo = x_array[index];
      dx = x_hi - x_lo;
      if (dx > 0.0)
	{
	  const double y_lo = y_array[index];
	  const double y_hi = y_array[index + 1];
	  const double dy = y_hi - y_lo;
	  double delx = x - x_lo;
	  double b_i, c_i, d_i;
	  COEFF_CALC(interp->c, dy, dx, index,  b_i, c_i, d_i);
	  *y_pp = 2.0 * c_i + 6.0 * d_i * delx;
	  return GSL_SUCCESS;
	}
      else
	{
	  *y_pp = 0.0;
	  return GSL_FAILURE;
	}
    }
}


static
int
cspline_eval_i_impl (const gsl_interp_obj * cspline_interp,
		     const double x_array[], const double y_array[],
		     gsl_interp_accel * acc,
                     double a, double b,
		     double * result)
{
  const gsl_interp_cspline *interp = (const gsl_interp_cspline *) cspline_interp;

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
	  double b_i, c_i, d_i; 
	  COEFF_CALC(interp->c, dy, dx, index_a,  b_i, c_i, d_i);
	  INTEG_EVAL(y_lo, b_i, c_i, d_i, x_lo, a, b, *result);
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
	  const double dy = y_hi - y_lo;
	  if(dx != 0.0) {
	    double b_i, c_i, d_i; 
	    COEFF_CALC(interp->c, dy, dx, i,  b_i, c_i, d_i);
	    *result += dx * (y_lo + dx*(0.5*b_i + dx*(c_i/3.0 + 0.25*d_i*dx)));
	  }
	  else {
	    *result = 0.0;
	    return GSL_FAILURE;
	  }
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
	    double b_i, c_i, d_i; 
	    double tmp;
	    COEFF_CALC(interp->c, dy, dx, index_a,  b_i, c_i, d_i);
	    INTEG_EVAL(y_lo, b_i, c_i, d_i, x_lo, a, x_hi, tmp);
	    *result += tmp;
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
	    double b_i, c_i, d_i; 
	    double tmp;
	    COEFF_CALC(interp->c, dy, dx, index_a,  b_i, c_i, d_i);
	    INTEG_EVAL(y_lo, b_i, c_i, d_i, x_lo, x_lo, b, tmp);
	    *result += tmp;
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
