/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <gsl_errno.h>
#include "bsearch.h"
#include "tridiag.h"
#include "gsl_interp.h"

#define locMax(a, b)  ((a) > (b) ? (a) : (b))
#define locMin(a, b)  ((a) < (b) ? (a) : (b))


/* cubic spline interpolation object */
typedef struct {
  int       (*eval_impl)   (const gsl_interp_obj *, const double xa[], const double ya[], double x, gsl_interp_accel *, double * y);
  int       (*eval_na_impl) (const gsl_interp_obj *, const double xa[], const double ya[], double x, double * y);
  void      (*free)        (gsl_interp_obj *);
  double    xmin;
  double    xmax;
  int       size;
  double *  c;
}
gsl_interp_cspline;


gsl_interp_obj *
gsl_interp_cspline_natural_create(const double xa[], const double ya[], int size);

gsl_interp_obj *
gsl_interp_cspline_periodic_create(const double xa[], const double ya[], int size);

gsl_interp_obj *
gsl_interp_cspline_notnode_create(const double xa[], const double ya[], int size);

gsl_interp_obj *
gsl_interp_cspline_fixed_create(const double xa[], const double ya[], int size);

void
gsl_interp_cspline_free(gsl_interp_obj * interp);

int
gsl_interp_cspline_eval_impl(const gsl_interp_obj *, const double xa[], const double ya[], double x, gsl_interp_accel *, double * y);

int
gsl_interp_cspline_eval_na_impl(const gsl_interp_obj *, const double xa[], const double ya[], double x, double * y);


/* global cubic spline factory objects */
const gsl_interp_obj_factory gsl_interp_cspline_natural_factory = {
  "cspline_natural",
  gsl_interp_cspline_natural_create
};
const gsl_interp_obj_factory gsl_interp_cspline_periodic_factory = {
  "cspline_periodic",
  gsl_interp_cspline_periodic_create
};
const gsl_interp_obj_factory gsl_interp_cspline_notnode_factory = {
  "cspline_notnode",
  gsl_interp_cspline_notnode_create
};
const gsl_interp_obj_factory gsl_interp_cspline_fixed_factory = {
  "cspline_fixed",
  gsl_interp_cspline_fixed_create
};


/* common initialization */
static
gsl_interp_cspline * cspline_new(const double xa[], int size)
{
  gsl_interp_cspline * interp = (gsl_interp_cspline *) malloc(sizeof(gsl_interp_cspline));
  if(interp != 0) {
    interp->eval_impl	= gsl_interp_cspline_eval_impl;
    interp->eval_na_impl = gsl_interp_cspline_eval_na_impl;
    interp->xmin = xa[0];
    interp->xmax = xa[size-1];
    interp->size = size;
    interp->c = (double *) malloc(size * sizeof(double));
    if(interp->c == 0) {
      free(interp);
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
cspline_calc_natural(gsl_interp_cspline * interp,
                     const double xa[], const double ya[]
                     )
{
  int i;
  int num_points   = interp->size;
  int max_index    = num_points - 1; /* Engeln-Mullges + Uhlig "n" */
  int sys_size     = max_index  - 1; /* linear system is sys_size x sys_size */
  int status;
  double * g       = (double *) malloc(sys_size * sizeof(double));
  double * diag    = (double *) malloc(sys_size * sizeof(double));
  double * offdiag = (double *) malloc(sys_size * sizeof(double));
  
  if(g == 0 || diag == 0 || offdiag == 0) {
    status = GSL_ENOMEM;
  }
  else {
    interp->c[0]         = 0.0;
    interp->c[max_index] = 0.0;

    for(i=1; i < sys_size; i++) {
      offdiag[i] = xa[i+1] - xa[i];
    }

    for(i=0; i < sys_size; i++) {
      double h_i   = xa[i+1] - xa[i];
      double h_ip1 = xa[i+2] - xa[i+1];
      diag[i] = 2.0 * (h_ip1 + h_i);
      g[i]    = 3.0*((ya[i+2] - ya[i+1])/h_ip1 - (ya[i+1] - ya[i])/h_i);
    }

    status = solve_tridiag(diag, offdiag, g, interp->c + 1, sys_size);
  }

  if(g != 0)       free(g);
  if(diag != 0)    free(diag);
  if(offdiag != 0) free(offdiag);
  return status;
}


/* periodic spline calculation
 * see [Engeln-Mullges + Uhlig, p. 256]
 */
static
int
cspline_calc_periodic(gsl_interp_cspline * interp,
                      const double xa[], const double ya[]
                      )
{
  int i;
  int num_points   = interp->size;
  int max_index    = num_points - 1; /* Engeln-Mullges + Uhlig "n" */
  int sys_size     = max_index;      /* linear system is sys_size x sys_size */
  int status;
  double * g       = (double *) malloc(sys_size * sizeof(double));
  double * diag    = (double *) malloc(sys_size * sizeof(double));
  double * offdiag = (double *) malloc(sys_size * sizeof(double));
  
  if(g == 0 || diag == 0 || offdiag == 0) {
    status = GSL_ENOMEM;
  }
  else {
    if(sys_size == 2) {
      double h0 = xa[1] - xa[0];
      double h1 = xa[2] - xa[1];
      double h2 = xa[3] - xa[2];
      offdiag[0] = offdiag[1] = h0 + h1;
      diag[0]    =    diag[1] = 2.0*(h0 + h1);
      g[0] = 3.0*((ya[2] - ya[1])/h1 - (ya[1] - ya[0])/h0);
      g[1] = 3.0*((ya[1] - ya[2])/h2 - (ya[2] - ya[1])/h1);
      /* FIXME: solve system */
      status == GSL_SUCCESS;  
    }
    else {
      for(i=0; i < sys_size; i++) {
        offdiag[i] = xa[i+1] - xa[i];
      }

      for(i=0; i < sys_size; i++) {
        double h_i   = xa[i+1] - xa[i];
        double h_ip1 = xa[i+2] - xa[i+1];
        diag[i] = 2.0 * (h_ip1 + h_i);
        g[i]    = 3.0*((ya[(i+2) % num_points] - ya[i+1])/h_ip1 - (ya[i+1] - ya[i])/h_i);
      }

      status = solve_cyctridiag(diag, offdiag, g, interp->c + 1, sys_size);
      interp->c[0] = interp->c[max_index];
    }
  }

  if(g != 0)       free(g);
  if(diag != 0)    free(diag);
  if(offdiag != 0) free(offdiag);
  return status;
}


/* factory method */
gsl_interp_obj *
gsl_interp_cspline_natural_create(const double x_array[], const double y_array[], int size)
{
  if(size <= 1)
    return 0;
  else {
    gsl_interp_cspline * interp = cspline_new(x_array, size);
    if(interp != 0) {
      cspline_calc_natural(interp, x_array, y_array);
    }
    return (gsl_interp_obj *) interp;
  }
}


/* factory method */
gsl_interp_obj *
gsl_interp_cspline_periodic_create(const double x_array[], const double y_array[], int size)
{
  if(size <= 1)
    return 0;
  else {
    gsl_interp_cspline * interp = cspline_new(x_array, size);
    if(interp != 0) {
      cspline_calc_periodic(interp, x_array, y_array);
    }
    return (gsl_interp_obj *) interp;
  }
}

/* factory method */
gsl_interp_obj *
gsl_interp_cspline_notnode_create(const double x_array[], const double y_array[], int size)
{
  if(size <= 1)
    return 0;
  else {
    gsl_interp_cspline * interp = cspline_new(x_array, size);
    /* FIXME: calculate for notnode case */
    return (gsl_interp_obj *) interp;
  }
}

void
gsl_interp_cspline_free(gsl_interp_obj * interp_cspline)
{
  gsl_interp_cspline * interp = (gsl_interp_cspline *) interp_cspline;
  
  if(interp != 0) {
    if(interp->c != 0) free(interp->c);
    free(interp);
  }
}


int
gsl_interp_cspline_eval_na_impl(const gsl_interp_obj * cspline_interp,
                                const double x_array[], const double y_array[],
                                double x,
                                double * y
                                )
{
  gsl_interp_cspline * interp = (gsl_interp_cspline * )cspline_interp;

  if(x > interp->xmax) {
    *y = y_array[interp->size - 1];
    return GSL_EDOM;
  }
  else if(x < interp->xmin) {
    *y = y_array[0];
    return GSL_EDOM;
  }
  else {
    int    index = bsearch(x_array, x, 0, interp->size - 1);
    double x_lo = x_array[index];
    double y_lo = y_array[index];
    double x_hi = x_array[index + 1];
    double y_hi = y_array[index + 1];
    double dx = x_hi - x_lo;
    double dy = y_hi - y_lo;
    if(dx > 0.0) {
      double c_ip1 = interp->c[index+1];
      double c_i   = interp->c[index];
      double b_i = dy/dx - dx*(c_ip1 + 2.0*c_i)/3.0;
      double d_i = (c_ip1 - c_i)/(3.0*dx);
      *y = y_array[index] + dx*(b_i + dx*(c_i + dx*d_i));
      return GSL_SUCCESS;
    }
    else {
      *y = 0.0;
      return GSL_FAILURE;
    }
  }
}


int
gsl_interp_cspline_eval_impl(const gsl_interp_obj * cspline_interp,
                             const double x_array[], const double y_array[],
		             double x,
			     gsl_interp_accel * a,
			     double * y
		             )
{
  gsl_interp_cspline * interp = (gsl_interp_cspline * ) cspline_interp;

  if(x < interp->xmin) {
    *y = y_array[0];
    return GSL_EDOM;
  }
  else if(x > interp->xmax) {
    *y = y_array[interp->size - 1];
    return GSL_EDOM;
  }
  else {
    int     index_lo = a->cache_lo;
    int     index_hi = a->cache_hi;
    int     index;
    double x_hi, x_lo, dx;

    /* find correct bin; check for accelerator cache hits */

    if(x < x_array[index_lo]) {
      /* cache miss: lo side */

      if(a->heuristic == GSL_INTERP_LOCALSTEP) {
        /* reverse step heuristic */
	index_lo = locMax(0, index_lo - a->cache_size);
	index_hi = locMin(interp->size - 1, index_lo + a->cache_size);
	if(x >= x_array[index_lo] && x <= x_array[index_hi]) {
	  ++a->hit_count;
	  index = CHECK_BSEARCH(x_array, x, index_lo, index_hi);
	}
	else {
	  ++a->miss_count;
	  index = bsearch(x_array, x, 0, index_lo);
	}
      }
      else {
        ++a->miss_count;
	index = bsearch(x_array, x, 0, index_lo);
      }
    }
    else if(x > x_array[index_hi]) {
      /* cache miss: hi side */

      if(a->heuristic == GSL_INTERP_LOCALSTEP) {
        /* forward step heuristic */
	index_hi = locMin(interp->size - 1, index_hi + a->cache_size);
	index_lo = locMax(0, index_hi - a->cache_size);
	if(x >= x_array[index_lo] && x <= x_array[index_hi]) {
	  ++a->hit_count;
	  index = CHECK_BSEARCH(x_array, x, index_lo, index_hi);
	}
	else {
	  ++a->miss_count;
	  index = bsearch(x_array, x, index_hi, interp->size - 1);
	}
      }
      else {
        ++a->miss_count;
	index = bsearch(x_array, x, index_hi, interp->size - 1);
      }
    }
    else {
      /* cache hit */
      ++a->hit_count;
      index = CHECK_BSEARCH(x_array, x, index_lo, index_hi);
    }

      
    /* adjust accelerator cache */
    index_lo = locMax(0, index - a->cache_size/2);
    index_hi = locMin(interp->size - 1, index_lo + a->cache_size);

    /* evaluate */
    x_hi = x_array[index + 1];
    x_lo = x_array[index];
    dx = x_hi - x_lo;
    if(dx > 0.0) {
      double y_lo = y_array[index];
      double y_hi = y_array[index + 1];
      double dy = y_hi - y_lo;
      double c_ip1 = interp->c[index+1];
      double c_i   = interp->c[index];
      double b_i = dy/dx - dx*(c_ip1 + 2.0*c_i)/3.0;
      double d_i = (c_ip1 - c_i)/(3.0*dx);
      *y = y_array[index] + dx*(b_i + dx*(c_i + dx*d_i));
      return GSL_SUCCESS;
    }
    else {
      *y = 0.0;
      return GSL_FAILURE;
    }
  }
}
