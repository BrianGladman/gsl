/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef ODEIV_UTIL_H
#define ODEIV_UTIL_H

#include <stdlib.h>
#include "gsl_odeiv.h"


/* Low-level allocator for gsl_odeiv_step objects.
 */
gsl_odeiv_step *
gsl_odeiv_step_new(
  const char * name,
  unsigned int dim,
  unsigned int ord,
  size_t state_size,
  size_t work_size);

#endif  /* !ODEIV_UTIL_H */
