/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <gsl_errno.h>
#include "bsearch.h"
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
    void (*free) (gsl_interp_obj *);
    double xmin;
    double xmax;
    size_t size;
    double *b;
    double *c;
    double *d;
  }
gsl_interp_akima;


gsl_interp_obj *
  gsl_interp_akima_natural_create (const double xa[], const double ya[], size_t size);

gsl_interp_obj *
  gsl_interp_akima_periodic_create (const double xa[], const double ya[], size_t size);


void
  gsl_interp_akima_free (gsl_interp_obj * interp);

int
  gsl_interp_akima_eval_impl (const gsl_interp_obj *, const double xa[], const double ya[], double x, gsl_interp_accel *, double *y);

int
  gsl_interp_akima_eval_d_impl (const gsl_interp_obj *, const double xa[], const double ya[], double x, gsl_interp_accel *, double *dydx);



const gsl_interp_factory gsl_interp_factory_akima_natural =
{
  "akima_natural",
  gsl_interp_akima_natural_create
};
const gsl_interp_factory gsl_interp_factory_akima_periodic =
{
  "akima_periodic",
  gsl_interp_akima_periodic_create
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
	  interp->eval_impl = gsl_interp_akima_eval_impl;
	  interp->eval_d_impl = gsl_interp_akima_eval_d_impl;
	  interp->free = gsl_interp_akima_free;
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


gsl_interp_obj *
gsl_interp_akima_natural_create (const double x_array[],
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

gsl_interp_obj *
gsl_interp_akima_periodic_create (const double x_array[],
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

void
gsl_interp_akima_free (gsl_interp_obj * akima_interp)
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


int
gsl_interp_akima_eval_impl (const gsl_interp_obj * akima_interp,
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
	  index = interp_bsearch (x_array, x, 0, interp->size - 1);
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


int
gsl_interp_akima_eval_d_impl (const gsl_interp_obj * akima_interp,
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
	  index = interp_bsearch (x_array, x, 0, interp->size - 1);
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
