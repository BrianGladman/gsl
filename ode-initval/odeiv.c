/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdlib.h>
#include <gsl_errno.h>
#include "gsl_odeiv.h"


int
gsl_odeiv_step_impl(
  gsl_odeiv_step * s,
  double t,
  double h,
  double y[],
  double yerr[],
  const gsl_odeiv_system * dydt)
{
  return s->_step(s->_state, s->_work, s->dimension, t, h, y, yerr, dydt);
}


int
gsl_odeiv_step_reset(gsl_odeiv_step * s)
{
  if(s->_reset != 0) {
    return s->_reset(s->_state);
  }
  else {
    return GSL_SUCCESS;
  }
}


void
gsl_odeiv_step_free(gsl_odeiv_step * s)
{
  if(s != 0) {
    if(s->_free != 0)  s->_free(s->_state, s->_work);
    if(s->_state != 0) free(s->_state);
    if(s->_work != 0)  free(s->_work);
    free(s);
  }
}
