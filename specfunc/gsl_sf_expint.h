/* Author: G. Jungman
 * RCS: $Id$
 */
#ifndef __GSL_EXPINT_H__
#define __GSL_EXPINT_H__

#include <gsl/gsl_sf_result.h>


/* E_1(x) := Re[ Integrate[ Exp[-xt]/t, {t,1,Infinity}] ]
 *
 * x != 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_expint_E1_impl(double x, gsl_sf_result * result);
int     gsl_sf_expint_E1_e(double x, gsl_sf_result * result);


/* E_2(x) := Re[ Integrate[ Exp[-xt]/t^2, {t,1,Infinity}] ]
 *
 * x != 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_expint_E2_impl(double x, gsl_sf_result * result);
int     gsl_sf_expint_E2_e(double x, gsl_sf_result * result);


/* Ei(x) := PV Integrate[ Exp[-t]/t, {t,-x,Infinity}]
 *
 * x != 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_expint_Ei_impl(double x, gsl_sf_result * result);
int     gsl_sf_expint_Ei_e(double x, gsl_sf_result * result);


/* Shi(x) := Integrate[ Sinh[t]/t, {t,0,x}]
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_Shi_impl(double x, gsl_sf_result * result);
int     gsl_sf_Shi_e(double x, gsl_sf_result * result);


/* Chi(x) := Re[ M_EULER + log(x) + Integrate[(Cosh[t]-1)/t, {t,0,x}] ]
 *
 * x != 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_Chi_impl(double x, gsl_sf_result * result);
int     gsl_sf_Chi_e(double x, gsl_sf_result * result);


/* Ei_3(x) := Integral[ Exp[-t^3], {t,0,x}]
 *
 * x >= 0.0
 * exceptions: GSL_EDOM
 */
int     gsl_sf_expint_3_impl(double x, gsl_sf_result * result);
int     gsl_sf_expint_3_e(double x, gsl_sf_result * result);


/* Si(x) := Integrate[ Sin[t]/t, {t,0,x}]
 *
 * exceptions: none
 */
int     gsl_sf_Si_impl(const double x, gsl_sf_result * result);
int     gsl_sf_Si_e(double x, gsl_sf_result * result);


/* Ci(x) := -Integrate[ Cos[t]/t, {t,x,Infinity}]
 *
 * x > 0.0
 * exceptions: GSL_EDOM 
 */
int     gsl_sf_Ci_impl(const double x, gsl_sf_result * result);
int     gsl_sf_Ci_e(double x, gsl_sf_result * result);


/* AtanInt(x) := Integral[ Arctan[t]/t, {t,0,x}]
 *
 *
 * exceptions:
 */
int gsl_sf_atanint_impl(double x, gsl_sf_result * result);
int gsl_sf_atanint_e(double x, gsl_sf_result * result);


#endif /* __GSL_EXPINT_H__ */
