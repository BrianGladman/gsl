/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef __GSL_SF_ELLJAC_H__
#define __GSL_SF_ELLJAC_H__

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


/* Jacobian elliptic functions sn, dn, cn,
 * by descending Landen transformations
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_elljac_impl(double u, double m, double * sn, double * cn, double * dn);
int gsl_sf_elljac_e(double u, double m, double * sn, double * cn, double * dn);


__END_DECLS

#endif /* __GSL_SF_ELLJAC_H__ */
