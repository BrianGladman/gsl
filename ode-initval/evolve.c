/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdlib.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_odeiv.h"


gsl_odeiv_evolve *
gsl_odeiv_evolve_new(void)
{
  gsl_odeiv_evolve * e = (gsl_odeiv_evolve *) malloc(sizeof(gsl_odeiv_evolve));
  if(e != 0) {
    e->y0 = 0;
    e->yerr = 0;
    e->count = 0;
    e->dimension = 0;
  }
  return e;
}


/* Evolution framework method.
 *
 * Uses an adaptive step control object
 * and/or a monitor object if given.
 */
int
gsl_odeiv_evolve_impl(
  gsl_odeiv_evolve * e, 
  gsl_odeiv_evolve_mon * mon,
  gsl_odeiv_evolve_control * con,
  gsl_odeiv_step * step,
  const gsl_odeiv_system * dydt,
  double t0, double t1, double hstart,
  double y[])
{
  double h = fabs(hstart);
  double t = t0;

  if(e == 0 || step == 0 || dydt == 0) return GSL_EFAULT;

  if(e->yerr == 0 || e->y0 == 0 || e->dimension != dydt->dimension) {
    if(e->yerr != 0) free(e->yerr);
    if(e->y0 != 0) free(e->y0);
    e->yerr = (double * ) malloc(dydt->dimension * sizeof(double));
    e->y0   = (double * ) malloc(dydt->dimension * sizeof(double));
    e->dydt_in  = (double * ) malloc(dydt->dimension * sizeof(double));
    e->dydt_out = (double * ) malloc(dydt->dimension * sizeof(double));
    e->dimension = dydt->dimension;

    if(e->yerr == 0 || e->y0 == 0 || e->dydt_in == 0 || e->dydt_out == 0) {
      if(e->yerr != 0) free(e->yerr);
      if(e->y0 != 0) free(e->y0);
      if(e->dydt_in != 0) free(e->dydt_in);
      if(e->dydt_out != 0) free(e->dydt_out);
      e->dimension = 0;
      e->yerr = 0;
      e->y0 = 0;
      e->dydt_in = 0;
      e->dydt_out = 0;
      return GSL_ENOMEM;
    }
  }

  e->count = 0;

  while(t < t1) {

    /* No need to copy if we cannot control the step size. */
    if(con!= 0) memcpy(e->y0, y, e->dimension * sizeof(double));

    /* Calculate initial dydt once if the method can benefit. */
    if(step->can_use_dydt) {
      GSL_ODEIV_FN_EVAL(dydt, t, y, e->dydt_in);
    }

    if(mon != 0 && mon->pre_step != 0) mon->pre_step(t, y);

    while(1) {
      const double t_save = t;
      int step_stat;

      h = GSL_MIN_DBL(t1-t, h);
      if(h < 8.0 * GSL_DBL_EPSILON) return GSL_EUNDRFLW; /* FIXME */

      if(step->can_use_dydt) {
        step_stat = gsl_odeiv_step_impl(step, t, h, y, e->yerr, e->dydt_in, e->dydt_out, dydt);
      }
      else {
        step_stat = gsl_odeiv_step_impl(step, t, h, y, e->yerr, 0, e->dydt_out, dydt);
      }
      t += h;

      if(con == 0) {
        /* We do not have the ability to adjust the step.
         * Continue with prescribed evolution.
         */
        break;
      }
      else {
        /* Check error and attempt to adjust the step. */
        const int hadj_stat = GSL_ODEIV_CONTROL_HADJ(con, e->dimension, step->order, y, e->yerr, e->dydt_out, &h);
	if(hadj_stat == GSL_ODEIV_HADJ_INC || hadj_stat == GSL_ODEIV_HADJ_NIL) {
	  /* Step was increased or unchanged. Continue with evolution. */
	  break;
	}
        else if(hadj_stat == GSL_ODEIV_HADJ_DEC) {
	  /* Step was decreased. Undo and go back to try again. */
          memcpy(y, y0, dydt->dimension * sizeof(double));
	  t = t_save;
	}
      }
    }

    ++e->count;
    if(mon != 0 && mon->post_step != 0) mon->post_step(t, y, e->yerr);
  }

  return GSL_SUCCESS;
}


void
gsl_odeiv_evolve_free(gsl_odeiv_evolve * e)
{
  if(e != 0) {
    if(e->y0 != 0) free(e->y0);
    if(e->yerr != 0) free(e->yerr);
    free(e);
  }
}
