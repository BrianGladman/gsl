/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef ODEIV_UTIL_H
#define ODEIV_UTIL_H

#include <stdlib.h>
#include "gsl_odeiv.h"


/* Low-level allocator for gsl_odeiv_step objects.
 */
 /*
gsl_odeiv_step *
gsl_odeiv_step_new(
  const char * name,
  unsigned int dim,
  unsigned int ord,
  size_t state_size,
  size_t work_size);
*/


void gsl_odeiv_step_construct(
  gsl_odeiv_step *,
  const char * _name,
  int  (*_step)  (void *, double, double, double *, double *, const double *, double *, const gsl_odeiv_system *),
  int  (*_reset) (void *),
  void (*_free)  (void *),
  int can_use_dydt,
  size_t dimension,
  unsigned int order
  );


#endif /* !ODEIV_UTIL_H */
