/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_gamma.h"
#include "gsl_sf_pow_int.h"
#include "gsl_sf_bessel.h"


/* sum which occurs in Taylor series for J_nu(x) * Gamma(nu+1)/((x/2)^nu)
 * [Abramowitz+Stegun, 9.1.10]
 */
int gsl_sf_bessel_Jnu_taylorsum_e(double nu, double x, int kmax, double * result)
{
  int k;
  double y = -0.25 * x*x;
  
  double kfact  = 1.;
  double nuprod = 1.;
  double yprod  = 1.;
  double ans    = 1.;
  double delta;
 
  for(k=1; k<=kmax; k++) {
    nuprod *= (nu + k);
    kfact  *= k;
    yprod  *= y;
    delta = yprod / (kfact * nuprod);
    ans += delta;
  }
  
  *result = ans;
  
  if(delta > 10. * GSL_MACH_EPS) {
    return GSL_ELOSS;
  }
  else {
    return GSL_SUCCESS;
  }
}


/* Taylor expansion for J_nu(x)
 * nu >= 0.0 and x >= 0.0
 * Return:
 *   GSL_SUCCESS on success
 *   GSL_EDOM    if nu < 0.0 or x < 0.0
 *   GSL_ELOSS   if loss of accuracy
 */
int gsl_sf_bessel_Jnu_taylor_e(double nu, double x, int kmax, double * result)
{
  if(nu < 0.0 || x < 0.0) {
    return GSL_EDOM;
  }
  else {
    double ts;
    double p   = (nu == 0. ? 1. :  pow(0.5*x, nu));
    double pre = pow(0.5*x, nu) / exp(gsl_sf_lngamma(nu+1.));
    int status = gsl_sf_bessel_Jnu_taylorsum_e(nu, x, kmax, &ts);
    *result = pre * ts;
    return status;
  }
}

int gsl_sf_bessel_Jn_taylor_e(int n, double x, int kmax, double * result)
{
  if(n < 0 || x < 0.0) {
    return GSL_EDOM;
  }  
  else {
    double ts;
    double pre = gsl_sf_pow_int(0.5*x, n) / exp(gsl_sf_lngamma(n+1.));
    int status = gsl_sf_bessel_Jnu_taylorsum_e(n, x, kmax, &ts);
    *result = pre * ts;
    return status;
  }
}

/************************************************************************
 *                                                                      *
  Asymptotic approximations 8.11.5, 8.12.5, and 8.42.7 from
  G.N.Watson, A Treatise on the Theory of Bessel Functions,
  2nd Edition (Cambridge University Press, 1944).
  Higher terms in expansion for x near l given by
  Airey in Phil. Mag. 31, 520 (1916).

  This approximation is accurate to near 0.1% at the boundaries
  between the asymptotic regions; well away from the boundaries
  the accuracy is better than 10^{-5}.
 *                                                                      *
 ************************************************************************/
double besselJ_meissel(double nu, double x)
{
  double beta = pow(nu, 0.325);
  double result;

  /* Fitted matching points.   */
  double llimit = 1.1 * beta;
  double ulimit = 1.3 * beta;

  double nu2 = nu * nu;

  if (nu < 5. && x < 1.)
    {
      /* Small argument and order. Use a Taylor expansion. */
      int k;
      double xo2 = 0.5 * x;
      double gamfactor = pow(nu,nu) * exp(-nu) * sqrt(nu * 2. * M_PI)
	* (1. + 1./(12.*nu) + 1./(288.*nu*nu));
      double prefactor = pow(xo2, nu) / gamfactor;
      double C[5];

      C[0] = 1.;
      C[1] = -C[0] / (nu+1.);
      C[2] = -C[1] / (2.*(nu+2.));
      C[3] = -C[2] / (3.*(nu+3.));
      C[4] = -C[3] / (4.*(nu+4.));
      
      result = 0.;
      for(k=0; k<5; k++)
	result += C[k] * pow(xo2, 2.*k);

      result *= prefactor;
    }
  else if(x < nu - llimit)
    {
      /* Small x region: x << l.    */
      double z = x / nu;
      double z2 = z*z;
      double rtomz2 = sqrt(1.-z2);
      double omz2_2 = (1.-z2)*(1.-z2);

      /* Calculate Meissel exponent. */
      double term1 = 1./(24.*nu) * ((2.+3.*z2)/((1.-z2)*rtomz2) -2.);
      double term2 = - z2*(4. + z2)/(16.*nu2*(1.-z2)*omz2_2);
      double V_nu = term1 + term2;
      
      /* Calculate the harmless prefactor. */
      double sterlingsum = 1. + 1./(12.*nu) + 1./(288*nu2);
      double harmless = 1. / (sqrt(rtomz2*2.*M_PI*nu) * sterlingsum);

      /* Calculate the logarithm of the nu dependent prefactor. */
      double ln_nupre = rtomz2 + log(z) - log(1. + rtomz2);

      result = harmless * exp(nu*ln_nupre - V_nu);
    } 
  else if(x < nu + ulimit)
    {         
      /* Intermediate region 1: x near nu. */
      double eps = 1.-nu/x;
      double eps_x = eps * x;
      double eps_x_2 = eps_x * eps_x;
      double xo6 = x/6.;
      double B[6];
      static double gam[6] = {2.67894, 1.35412, 1., 0.89298, 0.902745, 1.};
      static double sf[6] = {0.866025, 0.866025, 0., -0.866025, -0.866025, 0.};
      
      /* Some terms are identically zero, because sf[] can be zero.
       * Some terms do not appear in the result.
       */
      B[0] = 1.;
      B[1] = eps_x;
      /* B[2] = 0.5 * eps_x_2 - 1./20.; */
      B[3] = eps_x * (eps_x_2/6. - 1./15.);
      B[4] = eps_x_2 * (eps_x_2 - 1.)/24. + 1./280.;
      /* B[5] = eps_x * (eps_x_2*(0.5*eps_x_2 - 1.)/60. + 43./8400.); */

      result  = B[0] * gam[0] * sf[0] / pow(xo6, 1./3.);
      result += B[1] * gam[1] * sf[1] / pow(xo6, 2./3.);
      result += B[3] * gam[3] * sf[3] / pow(xo6, 4./3.);
      result += B[4] * gam[4] * sf[4] / pow(xo6, 5./3.);

      result /= (3.*M_PI);
    }
  else 
    {
      /* Region of very large argument. Use expansion
       * for x>>l, and we need not be very exacting.
       */
      double secb = x/nu;
      double sec2b= secb*secb;
      
      double cotb = 1./sqrt(sec2b-1.);      /* cotb=cot(beta) */

      double beta = acos(nu/x);
      double trigarg = nu/cotb - nu*beta - 0.25 * M_PI;
      
      double cot3b = cotb * cotb * cotb;
      double cot6b = cot3b * cot3b;

      double sum1, sum2, expterm, prefactor, trigcos;

      sum1  = 2.0 + 3.0 * sec2b;
      trigarg -= sum1 * cot3b / (24.0 * nu);

      trigcos = cos(trigarg);

      sum2 = 4.0 + sec2b;
      expterm = sum2 * sec2b * cot6b / (16.0 * nu2);

      expterm = exp(-expterm);
      prefactor = sqrt(2. * cotb / (nu * M_PI));
      
      result = prefactor * expterm * trigcos;
    }

  return  result;
}


double besselJprime_meissel(double nu, double x, double J_nu)
{
  return  besselJ_meissel(nu-1., x) - nu * J_nu / x;
}


double sphbesselj_meissel(double el, double x)
{
  return sqrt(M_PI/(2. * x)) * besselJ_meissel(el+0.5, x);
}


double sphbesseljprime_meissel(double el, double x, double j_l)
{
  return sphbesselj_meissel(el-1., x) - j_l * (el+1.)/x;
}


void asymp_sphbesselj_meissel(double l, double x,
			      double *jl, double *jlp, int dflag)
{
  double arg = x - (l+0.5)*0.5*M_PI;
  double x2 = x*x;
  double c = cos(arg);
  double s = sin(arg);

  *jl = (c + s/(8.*x))/x;
  if(dflag){
    *jlp = -(s*(8.*x2 + 2.) + 7.*x*c)/ (8.*x2*x);
  }
}



#define ACC GSL_MACH_EPS
#define RootPiOver2_  0.88622693
#define Gamma1pt5_    RootPiOver2_
void gsl_sf_bessel_j_steed(double x, int lmax, double * jl_x)
{
  if(x < ACC) {
    /* first term of Taylor series */
    int l;
    double inv_gam = 1./Gamma1pt5_;
    for(l=0; l<=lmax; l++) {
      jl_x[l] = RootPiOver2_ * pow(0.5*x, l) * inv_gam;
      inv_gam = inv_gam / (l+1.5);
    }
  }
  else {
    /* Steed/Barnett algorithm */
    double x_inv = 1./x;
    double W = 2.*x_inv;
    double F = 1.;
    double FP = (lmax+1.) * x_inv;
    double B = 2.*FP + x_inv;
    double end = B + 20000.*W;
    double D = 1./B;
    double del = -D;
    
    FP += del;
    
    /* continued fraction */
    do {
      B += W;
      D = 1./(B-D);
      del *= (B*D - 1.);
      FP += del;
      if(D < 0.) F = -F;
      if(B > end) {
	GSL_ERROR_RETURN(
		"gsl_sf_bessel_j_steed: continued fraction not converging",
		GSL_EFAILED,
		0.
		);
      }
    }
    while(fabs(del) >= fabs(FP) * ACC);
    
    FP *= F;
    
    if(lmax > 0) {
      /* downward recursion */
      double XP2 = FP;
      double PL = lmax * x_inv;
      int L  = lmax;
      int LP;
      jl_x[lmax] = F;
      for(LP = 1; LP<=lmax; LP++) {
	jl_x[L-1] = PL * jl_x[L] + XP2;
	FP = PL*jl_x[L-1] - jl_x[L];
	XP2 = FP;
	PL -= x_inv;
	--L;
      }
      F = jl_x[0];
    }
    
    /* normalization */
    W = x_inv / sqrt(FP*FP + F*F);
    jl_x[0] = W*F;
    if(lmax > 0) {
      int L;
      for(L=1; L<=lmax; L++) {
	jl_x[L] *= W;
      }
    }
  }
}
#undef ACC
#undef RootPiOver2_
#undef Gamma1pt5_
