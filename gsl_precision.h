/* Author:  B. Gough and G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_PRECISION_H
#define GSL_PRECISION_H
#include <gsl_machine.h>

/* Here are the constants which control the
 * precision/error system. These should be
 * expressed in terms of the underlying
 * machine constants, so this file is
 * platform-independent.
 */

/* Users can specify the desired precision
 * by passing one of the following. This
 * provides a simple and easily maintained
 * way to overload the functions for various
 * desired precisions.
 */
typedef
enum { GSL_APPROX_PREC=0, GSL_SINGLE_PREC=1, GSL_DOUBLE_PREC=2 }
gsl_prec_t;
#define _GSL_PREC_T_NUM 3


/* Arrays containing derived
 * precision constants for the
 * different precision levels.
 */
extern const double gsl_prec_eps[];
extern const double gsl_sqrt_prec_eps[];
extern const double gsl_root3_prec_eps[];
extern const double gsl_root4_prec_eps[];
extern const double gsl_root5_prec_eps[];
extern const double gsl_root6_prec_eps[];


#endif /* !GSL_PRECISION_H */
