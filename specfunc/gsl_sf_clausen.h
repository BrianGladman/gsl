/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef __GSL_SF_CLAUSEN_H__
#define __GSL_SF_CLAUSEN_H__

#include <gsl/gsl_sf_result.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS


/* Calculate the Clausen integral:
 *   Cl_2(x) := Integrate[-Log[2 Sin[t/2]], {t,0,x}]
 *
 * Relation to dilogarithm:
 *   Cl_2(theta) = Im[ Li_2(e^(i theta)) ]
 */
int gsl_sf_clausen_impl(double x, gsl_sf_result * result);
int gsl_sf_clausen_e(const double x, gsl_sf_result * result);


__END_DECLS

#endif /* __GSL_SF_CLAUSEN_H__ */
