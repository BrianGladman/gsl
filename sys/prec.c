/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_precision.h>

const double gsl_prec_eps[_GSL_PREC_T_NUM] = {
  GSL_DBL_EPSILON,
  GSL_FLT_EPSILON,
  GSL_SFLT_EPSILON
};

const double gsl_prec_sqrt_eps[_GSL_PREC_T_NUM] = {
  GSL_SQRT_DBL_EPSILON,
  GSL_SQRT_FLT_EPSILON,
  GSL_SQRT_SFLT_EPSILON
};

const double gsl_prec_root3_eps[_GSL_PREC_T_NUM] = {
  GSL_ROOT3_DBL_EPSILON,
  GSL_ROOT3_FLT_EPSILON,
  GSL_ROOT3_SFLT_EPSILON
};

const double gsl_prec_root4_eps[_GSL_PREC_T_NUM] = {
  GSL_ROOT4_DBL_EPSILON,
  GSL_ROOT4_FLT_EPSILON,
  GSL_ROOT4_SFLT_EPSILON
};

const double gsl_prec_root5_eps[_GSL_PREC_T_NUM] = {
  GSL_ROOT5_DBL_EPSILON,
  GSL_ROOT5_FLT_EPSILON,
  GSL_ROOT5_SFLT_EPSILON
};

const double gsl_prec_root6_eps[_GSL_PREC_T_NUM] = {
  GSL_ROOT6_DBL_EPSILON,
  GSL_ROOT6_FLT_EPSILON,
  GSL_ROOT6_SFLT_EPSILON
};


/* We need this somewhere, in case the inline is ignored.
 */

typedef unsigned int gsl_mode_t;

unsigned int GSL_MODE_PREC(gsl_mode_t mt);

unsigned int
GSL_MODE_PREC(gsl_mode_t mt)
{ 
  return  (mt & (unsigned int)7); 
}

