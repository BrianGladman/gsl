/* Author: G. Jungman
 * RCS: $Id$
 */
#ifndef GSL_ELLINT_H_
#define GSL_ELLINT_H_
#include <gsl_precision.h>


/* Legendre form of complete elliptic integrals
 *
 * K(k) = Integral[1/Sqrt[1 - k^2 Sin[t]^2], {t, 0, Pi/2}]
 * E(k) = Integral[  Sqrt[1 - k^2 Sin[t]^2], {t, 0, Pi/2}]
 *
 * exceptions: GSL_EDOM
 */
int     gsl_sf_ellint_Kcomp_impl(double k, double * result, gsl_prec_t precision_goal);
int     gsl_sf_ellint_Kcomp_e(double k, double * result, gsl_prec_t precision_goal);
double  gsl_sf_ellint_Kcomp(double k, gsl_prec_t precision_goal);

int     gsl_sf_ellint_Ecomp_impl(double k, double * result, gsl_prec_t precision_goal);
int     gsl_sf_ellint_Ecomp_e(double k, double * result, gsl_prec_t precision_goal);
double  gsl_sf_ellint_Ecomp(double k, gsl_prec_t precision_goal);


/* Legendre form of incomplete elliptic integrals
 *
 * F(phi,k)   = Integral[1/Sqrt[1 - k^2 Sin[t]^2], {t, 0, phi}]
 * E(phi,k)   = Integral[  Sqrt[1 - k^2 Sin[t]^2], {t, 0, phi}]
 * P(phi,k,n) = Integral[(1 + n Sin[t]^2)^(-1)/Sqrt[1 - k^2 Sin[t]^2], {t, 0, phi}]
 *
 * exceptions: GSL_EDOM
 */
int     gsl_sf_ellint_F_impl(double phi, double k, double * result, gsl_prec_t precision_goal);
int     gsl_sf_ellint_F_e(double phi, double k, double * result, gsl_prec_t precision_goal);
double  gsl_sf_ellint_F(double phi, double k, gsl_prec_t precision_goal);

int     gsl_sf_ellint_E_impl(double phi, double k, double * result, gsl_prec_t precision_goal);
int     gsl_sf_ellint_E_e(double phi, double k, double * result, gsl_prec_t precision_goal);
double  gsl_sf_ellint_E(double phi, double k, gsl_prec_t precision_goal);

int     gsl_sf_ellint_P_impl(double phi, double k, double n, double * result, gsl_prec_t precision_goal);
int     gsl_sf_ellint_P_e(double phi, double k, double n, double * result, gsl_prec_t precision_goal);
double  gsl_sf_ellint_P(double phi, double k, double n, gsl_prec_t precision_goal);

int     gsl_sf_ellint_D_impl(double phi, double k, double n, double * result, gsl_prec_t precision_goal);
int     gsl_sf_ellint_D_e(double phi, double k, double n, double * result, gsl_prec_t precision_goal);
double  gsl_sf_ellint_D(double phi, double k, double n, gsl_prec_t precision_goal);


/* Carlson's symmetric basis of functions
 *
 * RC(x,y)   = 1/2 Integral[(t+x)^(-1/2) (t+y)^(-1)], {t,0,Inf}]
 * RD(x,y,z) = 3/2 Integral[(t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-3/2), {t,0,Inf}]
 * RF(x,y,z) = 1/2 Integral[(t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-1/2), {t,0,Inf}]
 * RJ(x,y,z,p) = 3/2 Integral[(t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-1/2) (t+p)^(-1), {t,0,Inf}]
 *
 * exceptions: GSL_EDOM
 */
int     gsl_sf_ellint_RC_impl(double x, double y, double * result, gsl_prec_t precision_goal);
int     gsl_sf_ellint_RC_e(double x, double y, double * result, gsl_prec_t precision_goal);
double  gsl_sf_ellint_RC(double x, double y, gsl_prec_t precision_goal);

int     gsl_sf_ellint_RD_impl(double x, double y, double z, double * result, gsl_prec_t precision_goal);
int     gsl_sf_ellint_RD_e(double x, double y, double z, double * result, gsl_prec_t precision_goal);
double  gsl_sf_ellint_RD(double x, double y, double z, gsl_prec_t precision_goal);

int     gsl_sf_ellint_RF_impl(double x, double y, double z, double * result, gsl_prec_t precision_goal);
int     gsl_sf_ellint_RF_e(double x, double y, double z, double * result, gsl_prec_t precision_goal);
double  gsl_sf_ellint_RF(double x, double y, double z, gsl_prec_t precision_goal);

int     gsl_sf_ellint_RJ_impl(double x, double y, double z, double p, double * result, gsl_prec_t precision_goal);
int     gsl_sf_ellint_RJ_e(double x, double y, double z, double p, double * result, gsl_prec_t precision_goal);
double  gsl_sf_ellint_RJ(double x, double y, double z, double p, gsl_prec_t precision_goal);


#endif  /* !GSL_ELLINT_H_ */
