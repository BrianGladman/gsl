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
enum { GSL_APPROX_PREC, GSL_SINGLE_PREC, GSL_DOUBLE_PREC }
gsl_prec_t;


#define GSL_APPROX_EPS     5.0e-04

#define GSL_SINGLE_EPS     1.0e-07

#define GSL_DOUBLE_EPS       (3.0*GSL_MACH_EPS)


#endif /* !GSL_PRECISION_H */
