/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_errno.h>
#include "gsl_interp.h"


int
gsl_interp_eval_impl (const gsl_interp_obj * obj,
		      const double xa[], const double ya[], double x,
		      gsl_interp_accel * a, double *y)
{
  return obj->eval_impl (obj, xa, ya, x, a, y);
}

int
gsl_interp_eval_e (const gsl_interp_obj * obj,
		   const double xa[], const double ya[], double x,
		   gsl_interp_accel * a, double *y)
{
  int status = obj->eval_impl (obj, xa, ya, x, a, y);
  if (status != GSL_SUCCESS)
    {
      GSL_ERROR ("gsl_interp_eval_e", status);
    }
  return status;
}

double
gsl_interp_eval (const gsl_interp_obj * obj,
		 const double xa[], const double ya[], double x,
		 gsl_interp_accel * a)
{
  double y;
  int status = obj->eval_impl (obj, xa, ya, x, a, &y);
  if (status != GSL_SUCCESS)
    {
      GSL_WARNING ("gsl_interp_eval", status);
    }
  return y;
}


int
gsl_interp_eval_deriv_impl (const gsl_interp_obj * obj,
			    const double xa[], const double ya[], double x,
			    gsl_interp_accel * a,
			    double *dydx)
{
  return obj->eval_d_impl (obj, xa, ya, x, a, dydx);
}

int
gsl_interp_eval_deriv_e (const gsl_interp_obj * obj,
			 const double xa[], const double ya[], double x,
			 gsl_interp_accel * a,
			 double *dydx)
{
  int status = obj->eval_d_impl (obj, xa, ya, x, a, dydx);
  if (status != GSL_SUCCESS)
    {
      GSL_ERROR ("gsl_interp_eval_deriv_e", status);
    }
  return status;
}

double
gsl_interp_eval_deriv (const gsl_interp_obj * obj,
		       const double xa[], const double ya[], double x,
		       gsl_interp_accel * a)
{
  double dydx;
  int status = obj->eval_d_impl (obj, xa, ya, x, a, &dydx);
  if (status != GSL_SUCCESS)
    {
      GSL_WARNING ("gsl_interp_eval_deriv", status);
    }
  return dydx;
}


int
gsl_interp_eval_integ_impl (const gsl_interp_obj * obj,
			    const double xa[], const double ya[],
                            double a, double b,
			    gsl_interp_accel * acc,
			    double * result)
{
  return obj->eval_i_impl (obj, xa, ya, acc, a, b, result);
}

int
gsl_interp_eval_integ_e (const gsl_interp_obj * obj,
			 const double xa[], const double ya[],
                         double a, double b,
			 gsl_interp_accel * acc,
			 double * result)
{
  int status = obj->eval_i_impl (obj, xa, ya, acc, a, b, result);
  if (status != GSL_SUCCESS)
    {
      GSL_ERROR ("gsl_interp_eval_integ_e", status);
    }
  return status;
}

double
gsl_interp_eval_integ (const gsl_interp_obj * obj,
		       const double xa[], const double ya[],
                       double a, double b,
		       gsl_interp_accel * acc)
{
  double result;
  int status = obj->eval_i_impl (obj, xa, ya, acc, a, b, &result);
  if (status != GSL_SUCCESS)
    {
      GSL_WARNING ("gsl_interp_eval_integ", status);
    }
  return result;
}


void
gsl_interp_obj_free (gsl_interp_obj * obj)
{
  if (obj != 0)
    obj->free (obj);
}
