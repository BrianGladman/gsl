/* interpolation/akima.c
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
#include <math.h>
#include <gsl/gsl_errno.h>
#include "integ_eval_macro.h"
#include "gsl_interp.h"


/* akima interpolation object */
typedef struct
  {
    int (*eval_impl) (const gsl_interp_obj *,
		      const double xa[], const double ya[],
		      double x,
		      gsl_interp_accel *, double *y);
    int (*eval_d_impl) (const gsl_interp_obj *,
			const double xa[], const double ya[],
			double x,
			gsl_interp_accel *, double *dydx);
    int (*eval_d2_impl) (const gsl_interp_obj *,
			 const double xa[], const double ya[],
			 double x,
			 gsl_interp_accel *, double *y_pp);
    int (*eval_i_impl) (const gsl_interp_obj *,
			const double xa[], const double ya[],
			gsl_interp_accel *, double a, double b, double * result);
    void (*free) (gsl_interp_obj *);
    double xmin;
    double xmax;
    size_t size;
    double *b;
    double *c;
    double *d;
  }
gsl_interp_akima;


static
gsl_interp_obj *
akima_natural_create (const double xa[], const double ya[], size_t size);

static
gsl_interp_obj *
akima_periodic_create (const double xa[], const double ya[], size_t size);


static
void
akima_free (gsl_interp_obj * interp);

static
int
akima_eval_impl (const gsl_interp_obj *, const double xa[], const double ya[], double x, gsl_interp_accel *, double *y);

static
int
akima_eval_d_impl (const gsl_interp_obj *, const double xa[], const double ya[], double x, gsl_interp_accel *, double *y_p);

static
int
akima_eval_d2_impl (const gsl_interp_obj *, const double xa[], const double ya[], double x, gsl_interp_accel *, double *y_pp);

static
int
akima_eval_i_impl (const gsl_interp_obj *, const double xa[], const double ya[], gsl_interp_accel *, double, double, double *);


const gsl_interp_factory gsl_interp_factory_akima_natural =
{
  "akima_natural",
  akima_natural_create
};
const gsl_interp_factory gsl_interp_factory_akima_periodic =
{
  "akima_periodic",
  akima_periodic_create
};


/* common creation */
static
gsl_interp_akima *
interp_akima_new (const double x_array[], const double y_array[], size_t size)
{
  y_array = 0;			/* prevent warning about unused parameter */

  if (size <= 4)
    return 0;
  else
    {
      gsl_interp_akima *interp = (gsl_interp_akima *) malloc (sizeof (gsl_interp_akima));
      if (interp != 0)
	{
	  interp->eval_impl = akima_eval_impl;
	  interp->eval_d_impl = akima_eval_d_impl;
	  interp->eval_d2_impl = akima_eval_d2_impl;
	  interp->eval_i_impl = akima_eval_i_impl;
	  interp->free = akima_free;
	  interp->xmin = x_array[0];
	  interp->xmax = x_array[size - 1];
	  interp->size = size;
	  interp->b = (double *) malloc (size * sizeof (double));
	  interp->c = (double *) malloc (size * sizeof (double));
	  interp->d = (double *) malloc (size * sizeof (double));
	  if (interp->b == 0 || interp->c == 0 || interp->d == 0)
	    {
	      if (interp->b != 0)
		free (interp->b);
	      if (interp->c != 0)
		free (interp->c);
	      if (interp->d != 0)
		free (interp->d);
	      free (interp);
	      return 0;
	    }
	}
      return interp;
    }
}


/* common calculation */
static
void
interp_akima_calc (gsl_interp_akima * interp, const double x_array[], double *m)
{
  size_t i;
  for (i = 0; i < interp->size; i++)
    {
      double NE = fabs (m[i + 1] - m[i]) + fabs (m[i - 1] - m[i - 2]);
      if (NE == 0.0)
	{
	  interp->b[i] = m[i];
	  interp->c[i] = 0.0;
	  interp->d[i] = 0.0;
	}
      else
	{
	  double h_i = x_array[i + 1] - x_array[i];
	  double NE_next = fabs (m[i + 2] - m[i + 1]) + fabs (m[i] - m[i - 1]);
	  double alpha_i = fabs (m[i - 1] - m[i - 2]) / NE;
	  double alpha_ip1;
	  double tL_ip1;
	  if (NE_next == 0.0)
	    {
	      tL_ip1 = m[i];
	    }
	  else
	    {
	      alpha_ip1 = fabs (m[i] - m[i - 1]) / NE_next;
	      tL_ip1 = (1 - alpha_ip1) * m[i] + alpha_ip1 * m[i + 1];
	    }
	  interp->b[i] = (1 - alpha_i) * m[i - 1] + alpha_i * m[i];
	  interp->c[i] = (3 * m[i] - 2 * interp->b[i] - tL_ip1) / h_i;
	  interp->d[i] = (interp->b[i] + tL_ip1 - 2 * m[i]) / (h_i * h_i);
	}
    }
}

static
gsl_interp_obj *
akima_natural_create (const double x_array[],
		      const double y_array[],
		      size_t size)
{
  gsl_interp_akima *interp = interp_akima_new (x_array, y_array, size);

  if (interp != 0)
    {
      double *_m = (double *) malloc ((size + 4) * sizeof (double));
      if (_m != 0)
	{
	  double *m = _m + 2;
	  size_t i;
	  for (i = 0; i < size - 2; i++)
	    {
	      m[i] = (y_array[i + 1] - y_array[i]) / (x_array[i + 1] - x_array[i]);
	    }

	  /* non-periodic boundary conditions */
	  m[-2] = 3 * m[0] - 2 * m[1];
	  m[-1] = 2 * m[0] - m[1];
	  m[size - 1] = 2 * m[size - 2] - m[size - 3];
	  m[size] = 3 * m[size - 2] - 2 * m[size - 3];

	  interp_akima_calc (interp, x_array, m);
	  free (_m);
	}
      else
	{
	  free (interp);
	  return 0;
	}
    }

  return (gsl_interp_obj *) interp;
}

static
gsl_interp_obj *
akima_periodic_create (const double x_array[],
		       const double y_array[],
		       size_t size)
{
  gsl_interp_akima *interp = interp_akima_new (x_array, y_array, size);

  if (interp != 0)
    {
      double *_m = (double *) malloc ((size + 4) * sizeof (double));
      if (_m != 0)
	{
	  double *m = _m + 2;
	  size_t i;
	  for (i = 0; i < size - 2; i++)
	    {
	      m[i] = (y_array[i + 1] - y_array[i]) / (x_array[i + 1] - x_array[i]);
	    }

	  /* periodic boundary conditions */
	  m[-2] = m[size - 1 - 2];
	  m[-1] = m[size - 1 - 1];
	  m[size - 1] = m[0];
	  m[size] = m[1];

	  interp_akima_calc (interp, x_array, m);
	  free (_m);
	}
      else
	{
	  free (interp);
	  return 0;
	}
    }

  return (gsl_interp_obj *) interp;
}

static
void
akima_free (gsl_interp_obj * akima_interp)
{
  gsl_interp_akima *interp = (gsl_interp_akima *) akima_interp;
  if (interp != 0)
    {
      if (interp->b != 0)
	free (interp->b);
      if (interp->c != 0)
	free (interp->c);
      if (interp->d != 0)
	free (interp->d);
      free (interp);
    }
}


static
int
akima_eval_impl (const gsl_interp_obj * akima_interp,
		 const double x_array[], const double y_array[],
		 double x,
		 gsl_interp_accel * a,
		 double *y)
{
  const gsl_interp_akima *interp = (const gsl_interp_akima *) akima_interp;

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
      {
	double x_lo = x_array[index];
	double delx = x - x_lo;
	double b = interp->b[index];
	double c = interp->c[index];
	double d = interp->d[index];
	*y = y_array[index] + delx * (b + delx * (c + d * delx));
	return GSL_SUCCESS;
      }
    }
}


static
int
akima_eval_d_impl (const gsl_interp_obj * akima_interp,
		   const double x_array[], const double y_array[],
		   double x,
		   gsl_interp_accel * a,
		   double *dydx)
{
  const gsl_interp_akima *interp = (const gsl_interp_akima *) akima_interp;

  y_array = 0;			/* prevent warning about unused parameter */

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
      {
	double x_lo = x_array[index];
	double delx = x - x_lo;
	double b = interp->b[index];
	double c = interp->c[index];
	double d = interp->d[index];
	*dydx = b + delx * (2.0 * c + 3.0 * d * delx);
	return GSL_SUCCESS;
      }
    }
}


static
int
akima_eval_d2_impl (const gsl_interp_obj * akima_interp,
		    const double x_array[], const double y_array[],
		    double x,
		    gsl_interp_accel * a,
		    double *y_pp)
{
  const gsl_interp_akima *interp = (const gsl_interp_akima *) akima_interp;

  y_array = 0;			/* prevent warning about unused parameter */

  if (x < interp->xmin)
    {
      *y_pp = 0.0;
      return GSL_EDOM;
    }
  else if (x > interp->xmax)
    {
      *y_pp = 0.0;
      return GSL_EDOM;
    }
  else
    {
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
      {
	const double x_lo = x_array[index];
	const double delx = x - x_lo;
	const double c = interp->c[index];
	const double d = interp->d[index];
	*y_pp = 2.0 * c + 6.0 * d * delx;
	return GSL_SUCCESS;
      }
    }
}


static
int
akima_eval_i_impl (const gsl_interp_obj * akima_interp,
		   const double x_array[], const double y_array[],
		   gsl_interp_accel * acc,
                   double a, double b,
		   double * result)
{
  const gsl_interp_akima *interp = (const gsl_interp_akima *) akima_interp;

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
	const double dx = x_hi - x_lo;
	if(dx != 0.0) {
	  const double b_i = interp->b[index_a];
	  const double c_i = interp->c[index_a];
	  const double d_i = interp->d[index_a];
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
	  const double dx = x_hi - x_lo;
	  if(dx != 0.0) {
	    const double b_i = interp->b[index_a];
	    const double c_i = interp->c[index_a];
	    const double d_i = interp->d[index_a];
	    INTEG_EVAL(y_lo, b_i, c_i, d_i, x_lo, a, b, *result);
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
	  const double dx = x_hi - x_lo;
	  if(dx != 0.0) {
	    const double b_i = interp->b[index_a];
	    const double c_i = interp->c[index_a];
	    const double d_i = interp->d[index_a];
	    double tmp;
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
	  const double dx = x_hi - x_lo;
	  if(dx != 0.0) {
	    const double b_i = interp->b[index_a];
	    const double c_i = interp->c[index_a];
	    const double d_i = interp->d[index_a];
	    double tmp;
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
