/* Author: G. Jungman
 * RCS: $Id$
 */
#ifndef GSL_ELLINT_H_
#define GSL_ELLINT_H_


/* Legendre form of complete elliptic integrals */

int gsl_sf_ellint_Kcomp_e(double k, double prec, double * result);
int gsl_sf_ellint_Ecomp_e(double k, double prec, double * result);

double gsl_sf_ellint_Kcomp(double k, double prec);
double gsl_sf_ellint_Ecomp(double k, double prec);


/* Legendre form of incomplete elliptic integrals */

int gsl_sf_ellint_F_e(double phi, double k, double prec, double * result);
int gsl_sf_ellint_E_e(double phi, double k, double prec, double * result);
int gsl_sf_ellint_P_e(double phi, double k, double n, double prec, double * result);
int gsl_sf_ellint_D_e(double phi, double k, double n, double prec, double * result);

double gsl_sf_ellint_F(double phi, double k, double prec);
double gsl_sf_ellint_E(double phi, double k, double prec);
double gsl_sf_ellint_P(double phi, double k, double n, double prec);
double gsl_sf_ellint_D(double phi, double k, double n, double prec);



/* Carlson's symmetric basis of functions */

int gsl_sf_ellint_RC_e(double x, double y, double errtol, double * result);
int gsl_sf_ellint_RD_e(double x, double y, double z, double errtol, double * result);
int gsl_sf_ellint_RF_e(double x, double y, double z, double errtol, double * result);
int gsl_sf_ellint_RJ_e(double x, double y, double z, double p, double errtol, double * result);

double gsl_sf_ellint_RC(double x, double y, double prec);
double gsl_sf_ellint_RD(double x, double y, double z, double prec);
double gsl_sf_ellint_RF(double x, double y, double z, double prec);
double gsl_sf_ellint_RJ(double x, double y, double z, double p, double prec);



int gsl_sf_ellint_F_impl(double phi, double k, double prec, double * result);
int gsl_sf_ellint_E_impl(double phi, double k, double prec, double * result);
int gsl_sf_ellint_P_impl(double phi, double k, double n, double prec, double * result);
int gsl_sf_ellint_D_impl(double phi, double k, double n, double prec, double * result);

int gsl_sf_ellint_Kcomp_impl(double k, double prec, double * result);
int gsl_sf_ellint_Ecomp_impl(double k, double prec, double * result);

int gsl_sf_ellint_RC_impl(double x, double y, double prec, double * result);
int gsl_sf_ellint_RD_impl(double x, double y, double z, double prec, double * result);
int gsl_sf_ellint_RF_impl(double x, double y, double z, double prec, double * result);
int gsl_sf_ellint_RJ_impl(double x, double y, double z, double p, double prec, double * result);


#endif  /* !GSL_ELLINT_H_ */
