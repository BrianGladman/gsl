/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdlib.h>
#include <gsl_errno.h>
#include "tridiag.h"
#include "gsl_interp.h"


/* cubic spline interpolation object */
typedef struct
  {
    int (*eval_impl) (const gsl_interp_obj *, const double xa[], const double ya[], double x, gsl_interp_accel *, double *y);
    int (*eval_d_impl) (const gsl_interp_obj *, const double xa[], const double ya[], double x, gsl_interp_accel *, double *dydx);
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
cspline_calc_natural (gsl_interp_cspline * interp,
		      const double xa[], const double ya[]
)
{
  size_t i;
  size_t num_points = interp->size;
  size_t max_index = num_points - 1;	/* Engeln-Mullges + Uhlig "n" */
  size_t sys_size = max_index - 1;	/* linear system is sys_size x sys_size */
  int status;
  double *g = (double *) malloc (sys_size * sizeof (double));
  double *diag = (double *) malloc (sys_size * sizeof (double));
  double *offdiag = (double *) malloc (sys_size * sizeof (double));

  if (g == 0 || diag == 0 || offdiag == 0)
    {
      status = GSL_ENOMEM;
    }
  else
    {
      interp->c[0] = 0.0;
      interp->c[max_index] = 0.0;

      for (i = 1; i < sys_size; i++)
	{
	  offdiag[i] = xa[i + 1] - xa[i];
	}

      for (i = 0; i < sys_size; i++)
	{
	  double h_i = xa[i + 1] - xa[i];
	  double h_ip1 = xa[i + 2] - xa[i + 1];
	  diag[i] = 2.0 * (h_ip1 + h_i);
	  g[i] = 3.0 * ((ya[i + 2] - ya[i + 1]) / h_ip1 - (ya[i + 1] - ya[i]) / h_i);
	}

      status = solve_tridiag (diag, offdiag, g, interp->c + 1, sys_size);
    }

  if (g != 0)
    free (g);
  if (diag != 0)
    free (diag);
  if (offdiag != 0)
    free (offdiag);
  return status;
}


/* periodic spline calculation
 * see [Engeln-Mullges + Uhlig, p. 256]
 */
static
int
cspline_calc_periodic (gsl_interp_cspline * interp,
		       const double xa[], const double ya[]
)
{
  size_t i;
  size_t num_points = interp->size;
  size_t max_index = num_points - 1;	/* Engeln-Mullges + Uhlig "n" */
  size_t sys_size = max_index;	/* linear system is sys_size x sys_size */
  int status;
  double *g = (double *) malloc (sys_size * sizeof (double));
  double *diag = (double *) malloc (sys_size * sizeof (double));
  double *offdiag = (double *) malloc (sys_size * sizeof (double));

  if (g == 0 || diag == 0 || offdiag == 0)
    {
      status = GSL_ENOMEM;
    }
  else
    {
      if (sys_size == 2)
	{
	  double h0 = xa[1] - xa[0];
	  double h1 = xa[2] - xa[1];
	  double h2 = xa[3] - xa[2];
	  offdiag[0] = offdiag[1] = h0 + h1;
	  diag[0] = diag[1] = 2.0 * (h0 + h1);
	  g[0] = 3.0 * ((ya[2] - ya[1]) / h1 - (ya[1] - ya[0]) / h0);
	  g[1] = 3.0 * ((ya[1] - ya[2]) / h2 - (ya[2] - ya[1]) / h1);

	  /* FIXME: solve 2x2 system */
	  status = GSL_SUCCESS;
	}
      else
	{
	  for (i = 0; i < sys_size; i++)
	    {
	      offdiag[i] = xa[i + 1] - xa[i];
	    }

	  for (i = 0; i < sys_size; i++)
	    {
	      double h_i = xa[i + 1] - xa[i];
	      double h_ip1 = xa[i + 2] - xa[i + 1];
	      diag[i] = 2.0 * (h_ip1 + h_i);
	      g[i] = 3.0 * ((ya[(i + 2) % num_points] - ya[i + 1]) / h_ip1 - (ya[i + 1] - ya[i]) / h_i);
	    }

	  status = solve_cyctridiag (diag, offdiag, g, interp->c + 1, sys_size);
	  interp->c[0] = interp->c[max_index];
	}
    }

  if (g != 0)
    free (g);
  if (diag != 0)
    free (diag);
  if (offdiag != 0)
    free (offdiag);
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
	  cspline_calc_natural (interp, x_array, y_array);
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
	  cspline_calc_periodic (interp, x_array, y_array);
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
	  double y_lo = y_array[index];
	  double y_hi = y_array[index + 1];
	  double dy = y_hi - y_lo;
	  double c_ip1 = interp->c[index + 1];
	  double c_i = interp->c[index];
	  double b_i = dy / dx - dx * (c_ip1 + 2.0 * c_i) / 3.0;
	  double d_i = (c_ip1 - c_i) / (3.0 * dx);
	  double delx = x - x_lo;
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
		     double *dydx
)
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
	  double y_lo = y_array[index];
	  double y_hi = y_array[index + 1];
	  double dy = y_hi - y_lo;
	  double c_ip1 = interp->c[index + 1];
	  double c_i = interp->c[index];
	  double b_i = dy / dx - dx * (c_ip1 + 2.0 * c_i) / 3.0;
	  double d_i = (c_ip1 - c_i) / (3.0 * dx);
	  double delx = x - x_lo;
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
