/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef __GSL_SF_ELLJAC_H__
#define __GSL_SF_ELLJAC_H__


/* Jacobian elliptic functions sn, dn, cn,
 * by descending Landen transformations
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_elljac_impl(double u, double m, double * sn, double * cn, double * dn);
int gsl_sf_elljac_e(double u, double m, double * sn, double * cn, double * dn);


#endif /* __GSL_SF_ELLJAC_H__ */
