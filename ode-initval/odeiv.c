/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdlib.h>
#include <gsl_errno.h>
#include "gsl_odeiv.h"


int
gsl_odeiv_step_impl(gsl_odeiv_step * s, double t, double h, double y[], double yerr[], gsl_odeiv_function * dydt)
{
  return s->_step(s->_state, s->dimension, t, h, y, yerr, dydt);
}


int
gsl_odeiv_step_reset(gsl_odeiv_step * s)
{
  return s->_reset(s->_state);
}


void
gsl_odeiv_step_free(gsl_odeiv_step * s)
{
  if(s != 0) {
    s->_free(s->_state);
    free(s);
  }
}
