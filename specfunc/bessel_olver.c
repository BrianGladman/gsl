/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>


static double olver_B0_lt1(double t, double zeta)
{
  return -5./(48.*zeta*zeta) + t*(-3 + 5.*t*t)/(24.*sqrt(zeta));
}

static double olver_B0_gt1(double t, double minus_zeta)
{
  double mz2 = minus_zeta*minus_zeta;
  return -5./(48.*mz2) + t*(3 + 5.*t*t)/(24.*sqrt(minus_zeta));
}

static double olver_A1_lt1(double t, double zeta)
{
  double rz = sqrt(zeta);
  double t2 = t*t;
  double term1 = t2*(81. - 462.*t2 + 385.*t2*t2)/1152.;
  double term2 = -455./(4608.*zeta*zeta*zeta);
  double term3 = 7.*t*(-3 + 5.*t2)/(1152.*rz*rz*rz);
  return term1 + term2 + term3;
}

static double olver_A1_gt1(double t, double minus_zeta)
{
  double rz = sqrt(minus_zeta);
  double t2 = t*t;
  double term1 = -t2*(81. + 462.*t2 + 385.*t2*t2)/1152.;
  double term2 =  455./(4608.*minus_zeta*minus_zeta*minus_zeta);
  double term3 = -7.*t*(3 + 5.*t2)/(1152.*rz*rz*rz);
  return term1 + term2 + term3;
}

static double olver_Asum_lt1(double nu, double t, double zeta)
{
  return 1. + olver_A1_lt1(t, zeta)/(nu*nu);
}

static double olver_Bsum_lt1(double nu, double t, double zeta)
{
  return olver_B0_lt1(t, zeta);
}


int gsl_sf_bessel_Jnu_asymp_Olver_impl(double nu, double x, double * result)
{
  double pre;
  double asum, bsum;
  double ai, aip;
  double z = x/nu;
  double crnu = pow(nu, 1./3.);
  double rt   = sqrt(fabs(1.-z)*(1+z));

  if(z < 1.) {
    double zeta = pow(1.5*(log((1+rt)/z) - rt), 2./3.);
    double arg  = crnu*crnu *zeta;
    pre  = sqrt(sqrt(4.*zeta/(rt*rt)));
    asum = olver_Asum_lt1(nu, 1./rt, zeta);
    bsum = olver_Bsum_lt1(nu, 1./rt, zeta);
    gsl_sf_airy_Ai_impl(arg, &ai);
    gsl_sf_airy_Ai_deriv_impl(arg, &aip);
  }
  else if(z > 1.) {
    double minus_zeta = pow(1.5*(rt - acos(1./z)), 2./3.);
    double arg  = crnu*crnu *(-minus_zeta);
    pre  = sqrt(sqrt(-4.*minus_zeta/(rt*rt)));
    asum = olver_Asum_gt1(nu, 1./rt, zeta);
    bsum = olver_Bsum_gt1(nu, 1./rt, zeta);
    gsl_sf_airy_Ai_impl(arg, &ai);
    gsl_sf_airy_Ai_deriv_impl(arg, &aip);
  }

  *result = pre * (ai*asum/crnu + aip*bsum/(nu*crnu*crnu));
  return GSL_SUCCESS;
}

