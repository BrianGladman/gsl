/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_ELLJAC_H_
#define GSL_SF_ELLJAC_H_


/* Jacobian elliptic functions sn, dn, cn,
 * by descending Landen transformations
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_elljac_impl(double u, double m, double * sn, double *cn, double * dn);
int gsl_sf_elljac_e(double u, double m, double * sn, double *cn, double * dn);


#endif  /* !GSL_SF_ELLJAC_H_ */
