/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_CLAUSEN_H_
#define GSL_SF_CLAUSEN_H_


/* Calculate the Clausen integral:
 *   Cl_2(x) := Integrate[-Log[2 Sin[t/2]], {t,0,x}]
 *
 * Relation to dilogarithm:
 *  Cl_2(theta) = Im[ Li_2(e^(i theta)) ]
 */
int     gsl_sf_clausen_impl(double x, double * result);
int     gsl_sf_clausen_e(double x, double * result);
double  gsl_sf_clausen(double x);



#endif /* !GSL_SF_CLAUSEN_H_ */
