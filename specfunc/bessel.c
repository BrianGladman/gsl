#include <math.h>
#include <gsl_errno.h>
#include "gsl_sf_bessel.h"


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


double gsl_sf_bessel_I0(double x)
{
  double ax = fabs(x);
  double y, ans;

  if(ax < 3.75) {
    
    /* Polynomial fit for small argument. */
    y = x / 3.75;
    y *= y;
    ans = 1. + y * (3.5156229 
		    + y * (3.0899424
			   + y * (1.2067492
				  + y * (0.2659732
					 + y * (0.360768e-1
						+ y * 0.45813e-2)
					 )
				  )
			   )
		    );
  }
  else {
    
    /* Prefactor with asymptotic correction. */
    y = 3.75 / ax;
    ans = exp(ax) / sqrt(ax);
    ans *= 0.39894228 
      + y * (0.1328592e-1
	     + y * (0.225319e-2
		    + y * (-0.157565e-2
			   + y * (0.916281e-2
				  + y * (-0.2057706e-1
					 + y * (0.2635537e-1
						+ y * (-0.1647633e-1
						       + y * 0.392377e-2)
						)
					 )
				  )
			   )
		    )
	     );
  }

  return ans;
}


double gsl_sf_bessel_I1(double x)
{
  double ax = fabs(x);
  double ans, y;

  if( (ax=fabs(x)) < 3.75){

    /* Polynomial approximation. */
    y = x / 3.75;
    y *= y;
    ans = ax * (0.5 + y * (0.87890594
			   + y * (0.51498869
				  + y * (0.15084934
					 + y * (0.2658733e-1
						+ y * (0.301532e-2
						       + y * 0.32411e-3)
						)
					 )
				  )
			   )
		);
  }
  else {

    /* Prefactor with asymptotic correction. */
    y = 3.75 / ax;
    ans = 0.2282967e-1 + y * (-0.2895312e-1
			      + y * (0.1787654e-1
				     - y * 0.420059e-2)
			      );
    ans = 0.39894228 + y * (-.3988024e-1
			    + y * (-0.362018e-2
				   + y * (0.163801e-2
					  + y * (-0.1031555e-1
						 + y * ans)
					  )
				   )
			    );
    ans *= (exp(ax) / sqrt(ax));
  }
  
  return x < 0. ? -ans : ans;
}


#define ACC 40.
#define BIGNO 1.e+10
#define BIGNI 1.e-10

double gsl_sf_bessel_I(int n, double x)
{
  if(n < 2){
    GSL_MESSAGE("besselI: n<2 in besselI(n,x)\n");
    return 0.;
  }
  
  if(x == 0.)
    return 0.;
  else {
    int j;
    double tox = 2./fabs(x);
    double bip = 0.;
    double ans = 0.;
    double bi = 1.;
    double bim;
    
    /* Downward recursion. */
    for(j=2*(n+(int) sqrt(ACC*n)); j>0; j--){
      bim = bip + j * tox * bi;
      bip = bi;
      bi = bim;
      
      /* Renormalize to prevent overflow. */
      if(fabs(bi) > BIGNO){
	ans *= BIGNI;
	bi *= BIGNI;
	bip *= BIGNI;
      }
      
      if(j == n) ans = bip;
    }

    ans *= gsl_sf_bessel_I0(x) / bi;

    return x < 0. && (n & 1) ? -ans : ans;
  }
}
#undef ACC
#undef BIGNO
#undef BIGNI


double gsl_sf_log_bessel_I0(double x)
{
  double ax = fabs(x);
  double y, ans;
  double poly;

  if(x <= 0.){
    GSL_MESSAGE("log_besselI0: log(I_0(x)) with x <= 0.\n");
    return 0.;
  }

  if(ax < 3.75) {
    
    /* Polynomial fit for small argument. */
    y = x / 3.75;
    y *= y;
    ans = 1. + y * (3.5156229 
		    + y * (3.0899424
			   + y * (1.2067492
				  + y * (0.2659732
					 + y * (0.360768e-1
						+ y * 0.45813e-2)
					 )
				  )
			   )
		    );
    ans = log(ans);
  }
  else {
    
    /* Prefactor with asymptotic correction. */
    y = 3.75 / ax;
    ans = ax - 0.5 * log(ax);
    poly = 0.39894228 
      + y * (0.1328592e-1
	     + y * (0.225319e-2
		    + y * (-0.157565e-2
			   + y * (0.916281e-2
				  + y * (-0.2057706e-1
					 + y * (0.2635537e-1
						+ y * (-0.1647633e-1
						       + y * 0.392377e-2)
						)
					 )
				  )
			   )
		    )
	     );
    ans += log(poly);
  }
  
  return ans;
}


double gsl_sf_log_bessel_I1(double x)
{
  double ax = fabs(x);
  double ans, y;
  double poly1, poly2;

  if( x <= 0.){
    GSL_MESSAGE("log_besselI1: log(I_1(x)), x <= 0.\n");
    return 0.;
  }

  if(ax < 3.75){

    /* Polynomial approximation. */
    y = x / 3.75;
    y *= y;
    ans = ax * (0.5 + y * (0.87890594
			   + y * (0.51498869
				  + y * (0.15084934
					 + y * (0.2658733e-1
						+ y * (0.301532e-2
						       + y * 0.32411e-3)
						)
					 )
				  )
			   )
		);
    ans = log(ans);
  }
  else {

    /* Prefactor with asymptotic correction. */
    y = 3.75 / ax;
    poly1 = 0.2282967e-1 + y * (-0.2895312e-1
			      + y * (0.1787654e-1
				     - y * 0.420059e-2)
			      );
    poly2 = 0.39894228 + y * (-.3988024e-1
			    + y * (-0.362018e-2
				   + y * (0.163801e-2
					  + y * (-0.1031555e-1
						 + y * poly1)
					  )
				   )
			    );
    ans = ax - 0.5 * log(ax) + log(poly2);
  }
  
  return x < 0. ? -ans : ans;
}


#define ACC 40.
#define BIGNO 1.e+10
#define BIGNI 1.e-10

double gsl_sf_log_bessel_I(int n, double x)
{
 if(n < 2){
    GSL_MESSAGE("log_besselI: n<2 in log_besselI(n,x)");
    return 0.;
  }
  if(x <= 0.){
    GSL_MESSAGE("log_besselI: log(I_n(x)) with x <= 0.");
    return 0.;
  }
  else {
    int j;
    double logfactor;
    double tox = 2./fabs(x);
    double bip = 0.;
    double factor = 0.;
    double bi = 1.;
    double bim;
    double ans;
    
    /* Downward recursion. */
    for(j=2*(n+(int) sqrt(ACC*n)); j>0; j--){
      bim = bip + j * tox * bi;
      bip = bi;
      bi = bim;
      
      /* Renormalize to prevent overflow. */
      if(fabs(bi) > BIGNO){
	factor *= BIGNI;
	bi  *= BIGNI;
	bip *= BIGNI;
      }
      
      if(j == n) factor = bip;
    }

    /* Cutoff the real low values. This is
     * mainly to trap for log(0).
     */
    if(factor <= 1.e-200)
      logfactor = log(1.e-200);
    else
      logfactor = log(factor);

    ans = logfactor + log_besselI0(x) - log(bi);

    return ans;
  }
}
#undef ACC
#undef BIGNO
#undef BIGNI


#define ACC 1.e-14
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
	char buff[100];
	sprintf(buff, "gsl_sf_bessel_j_steed: continued fraction not converging");
	GSL_MESSAGE(buff);
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
	FP = PL*jl_x[L-1]         - jl_x[L];
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
