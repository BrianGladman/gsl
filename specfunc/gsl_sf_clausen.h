/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_CLAUSEN_H_
#define GSL_SF_CLAUSEN_H_

#include <gsl/gsl_sf_result.h>


/* Calculate the Clausen integral:
 *   Cl_2(x) := Integrate[-Log[2 Sin[t/2]], {t,0,x}]
 *
 * Relation to dilogarithm:
 *   Cl_2(theta) = Im[ Li_2(e^(i theta)) ]
 */
int gsl_sf_clausen_impl(double x, gsl_sf_result * result);
int gsl_sf_clausen_e(double x, gsl_sf_result * result);


#endif /* !GSL_SF_CLAUSEN_H_ */
