/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"

#include "airy_impl.h"
#include "bessel.h"

#define CubeRoot2_ 1.25992104989487316476721060728


/* Chebyshev fit for f(x) = A_2(2/(x+1)) / (1+x) */
static double A2_gt1_data[28] = {
  0.00042740836099373760083877383820,
  0.000198015693944389527785047335096,
 -0.000054454552746900308473446906559,
 -0.0000222673228383085134964757854674,
  0.0000152232502016475221849527768055,
 -3.3243013103482365073560071893e-6,
 -5.3060213036097467515440647626e-7,
  8.0934253293797313682621766902e-7,
 -4.3970412810377017046814010185e-7,
  1.79071424519372576492608662140e-7,
 -6.2467398631792389300175040895e-8,
  1.96894735620475710289112015304e-8,
 -5.7675858040046097235060721371e-9,
  1.59728455317871235349735453543e-9,
 -4.2302024090656279644231530804e-10,
  1.08004643558298536450496126138e-10,
 -2.67437492468402217094751959285e-11,
  6.4518542339596162421910187855e-12,
 -1.52189607648056016027864364660e-12,
  3.5202190883266175011615161782e-13,
 -8.0030236559195190012512311175e-14,
  1.79175234158451644829625971592e-14,
 -3.9568041553037887210172021266e-15,
  8.6307619819767052015119074384e-16,
 -1.86166067383562915065864440616e-16,
  3.9750004549810259594153022435e-17,
 -8.4088913756029427849463017114e-18,
  1.76364042377515061860916815522e-18
};
static struct gsl_sf_ChebSeries A2_gt1_cs = {
  A2_gt1_data,
  -1, 1,
  27,
  (double *)0,
  (double *)0
};


/* Chebyshev fit for f(x) = A_3(2/(x+1)) / (1+x)^3 */
static double A3_gt1_data[24] = {
  -0.000105857038770346154791201255476,
  -0.0000205760351047268016023585098044,
   0.000038843887932590651228554092924,
  -6.4950843727032906490961018776e-6,
  -6.8550260909371744214172673062e-6,
   5.2395352749196099765976515407e-6,
  -1.69616633544227242227774179427e-6,
   1.72521003343504841706687755356e-8,
   3.3883914247036941743279709445e-7,
  -2.57613705551957452862302913685e-7,
   1.33889960634289278564724992969e-7,
  -5.7713728962127577476860031398e-8,
   2.20603248663723535274887935371e-8,
  -7.7303494297479242603927953267e-9,
   2.53299507665549622149741647084e-9,
  -7.8631261188698187150310486463e-10,
   2.33399931618284706272149336039e-10,
  -6.6701601102133563477507979861e-11,
   1.84501510088942544841583986784e-11,
  -4.9603190140521104131104222919e-12,
   1.30057299865001322997537278176e-12,
  -3.3349097589870270196152533718e-13,
   8.3823663352585629972074916939e-14,
  -2.06935497620521059555227840228e-14
};
static struct gsl_sf_ChebSeries A3_gt1_cs = {
  A3_gt1_data,
  -1, 1,
  23,
  (double *)0,
  (double *)0
};

/* Chebyshev fit for f(x) = B_1(2/(x+1)) */
static double B1_gt1_data[30] = {
  -0.00125638387443615968469140078230,
  -0.00079518496830875387886042742516,
  -0.000117550533792019687221797909112,
   0.000053177832099392027220254310301,
  -2.69942462974709280278958305596e-6,
  -3.8442655040391955801226902901e-6,
   1.93735278352654820441613088112e-6,
  -5.7567177444348444716886594860e-7,
   1.09655991539838265398335388742e-7,
  -3.9493463262644952802342029354e-10,
  -1.29925137787315793022146720225e-8,
   9.0685864569780974583429213845e-9,
  -4.9217596347029779922435454110e-9,
   2.52379362404573993294129066226e-9,
  -1.30781535782041343419561644550e-9,
   7.0384726336811679054311137949e-10,
  -3.9656903287869259099757487561e-10,
   2.33557792366149146495284096555e-10,
  -1.43039836294329736948171530669e-10,
   9.0574939745376874544020136846e-11,
  -5.8990591168315025785999930659e-11,
   3.9336908506031143837035749223e-11,
  -2.67442765030341217028365302040e-11,
   1.84575544800677159001265401512e-11,
  -1.28626375497176648076677290858e-11,
   8.9834445611871474998305912925e-12,
  -6.2114446574034950328703160154e-12,
   4.1547451016435261417400468910e-12,
  -2.55075855193449883493817495022e-12,
   1.21351807208778202835648785394e-12
};
static struct gsl_sf_ChebSeries B1_gt1_cs = {
  B1_gt1_data,
  -1, 1,
  29,
  (double *)0,
  (double *)0
};

/* Chebyshev fit for f(x) = B_2(2/(x+1)) / (1+x)^2 */
static double B2_gt1_data[28] = {
  0.000171948661718211851305571216621,
  0.000079514741279756654016155767774,
 -0.0000231980045534389356735443835521,
 -9.6608078112814154253915121642e-6,
  6.9234246115168214059342873692e-6,
 -1.37541247706222772353756248204e-6,
 -4.3232083748880289164411169258e-7,
  4.7506644879380763788488850903e-7,
 -2.28714945515401399251738734230e-7,
  7.6411016497329937280971501352e-8,
 -1.68648036990450591006180459630e-8,
 -4.7817821661397417398426456903e-12,
  2.77300632054028838759266566595e-9,
 -2.20212690117293776747358435428e-9,
  1.33978703926356648679100873427e-9,
 -7.5528675426754967540212371899e-10,
  4.2211557074388851904570120605e-10,
 -2.41136554272714597841816511490e-10,
  1.42538905397504612817193580626e-10,
 -8.7372629703949883171709170161e-11,
  5.5382078895910144286543982412e-11,
 -3.6132035104786592416603496209e-11,
  2.41397234409420703745082294584e-11,
 -1.64282187755817239822773509280e-11,
  1.13202145709942687808326952275e-11,
 -7.8352182544028993620167163596e-12,
  5.3791924356175443953371734348e-12,
 -3.5787357658563595558317302866e-12
};
static struct gsl_sf_ChebSeries B2_gt1_cs = {
  B2_gt1_data,
  -1, 1,
  27,
  (double *)0,
  (double *)0
};


static double olver_B0(double z, double abs_zeta)
{
  if(fabs(1.-z) < 10.*GSL_SQRT_MACH_EPS) {
    return CubeRoot2_ * (1./70.  + 2./225.*(1-z));
  }
  else if(z < 1.) {
    double t = 1./sqrt((1-z)*(1+z));
    return -5./(48.*abs_zeta*abs_zeta) + t*(-3 + 5.*t*t)/(24.*sqrt(abs_zeta));
  }
  else {
    double t = 1./sqrt((z-1)*(z+1));
    return -5./(48.*abs_zeta*abs_zeta) + t*( 3 + 5.*t*t)/(24.*sqrt(abs_zeta));
  }
}

static double olver_B1(double z, double abs_zeta)
{
  if(fabs(1.-z) < 10.*GSL_SQRT_MACH_EPS) {
    return 0.5*CubeRoot2_ * (-1213./511875. - 3757./1347500.*(1-z));
  }
  else if(z < 1.) {
    return 0.; /* FIXME */
  }
  else {
    return gsl_sf_cheb_eval(2./z - 1., &B1_gt1_cs);
  }
}

static double olver_B2(double z, double abs_zeta)
{
  if(fabs(1.-z) < 10.*GSL_SQRT_MACH_EPS) {
    return 0.00055221307672129279005986982501
           + 0.00089586516310476929281129228969 * (1. - z);
  }
  else if(z < 1.) {
    return 0.; /* FIXME */
  }
  else {
    double x  = 2./z - 1.;
    double f2 = (1+x)*(1+x);
    return f2 * gsl_sf_cheb_eval(x,&B2_gt1_cs);
  }
}

static double olver_B3(double z, double abs_zeta)
{
  if(fabs(1.-z) < 10.*GSL_SQRT_MACH_EPS) {
    return -0.00047461779655995980754441833105
           -0.44003318369401117358772369675 * (1. - z);
  }
  else if(z < 1.) {
    return 0.; /* FIXME */
  }
  else {
    /* FIXME */
  }
}

static double olver_A1(double z, double abs_zeta)
{
  if(fabs(1.-z) < 10.*GSL_SQRT_MACH_EPS) {
    return -1./255. - 71./38500. * (1-z);
  }
  else if(z < 1.){
    double t = 1./sqrt((1-z)*(1+z));
    double rz = sqrt(abs_zeta);
    double t2 = t*t;
    double term1 =  t2*(81. - 462.*t2 + 385.*t2*t2)/1152.;
    double term2 = -455./(4608.*abs_zeta*abs_zeta*abs_zeta);
    double term3 =  7.*t*(-3 + 5.*t2)/(1152.*rz*rz*rz);
    return term1 + term2 + term3;
  }
  else {
    double t = 1./sqrt((z-1)*(z+1));
    double rz = sqrt(abs_zeta);
    double t2 = t*t;
    double term1 = -t2*(81. + 462.*t2 + 385.*t2*t2)/1152.;
    double term2 =  455./(4608.*abs_zeta*abs_zeta*abs_zeta);
    double term3 = -7.*t*( 3 + 5.*t2)/(1152.*rz*rz*rz);
    return term1 + term2 + term3;
  }
}

static double olver_A2(double z, double abs_zeta)
{
  if(fabs(1.-z) < 10.*GSL_SQRT_MACH_EPS) {
    return 151439./218295000. + 68401./147262500.*(1-z);
  }
  else if(z < 1.){
    /* FIXME */
  }
  else {
    double x  = 2./z -1.;
    return (1+x) * gsl_sf_cheb_eval(x, &A2_gt1_cs);
  }
}

static double olver_A3(double z, double abs_zeta)
{
  if(fabs(1.-z) < 10.*GSL_SQRT_MACH_EPS) {
    return -0.00035421197145774384077112575920
           -0.000312322527890318832782774881353 * (1. - z);
  }
  else if(z < 1.){
    /* FIXME */
  }
  else {
    double x = 2./z -1.;
    double f = (1+x)*(1+x)*(1+x);
    return f * gsl_sf_cheb_eval(x, &A3_gt1_cs);
  }
}

static double olver_A4(double z, double abs_zeta)
{
  if(fabs(1.-z) < 10.*GSL_SQRT_MACH_EPS) {
    return 0.00037819419920177291402661228437
          +0.00040494390552363233477213857527 * (1.-z);
  }
  else if(z < 1.){
    /* FIXME */
  }
  else {
    /* FIXME */
  }
}

static double olver_Asum(double nu, double z, double abs_zeta)
{
  double nu2 = nu*nu;
  double A1 = olver_A1(z, abs_zeta);
  double A2 = olver_A2(z, abs_zeta);
  double A3 = olver_A3(z, abs_zeta);
  return 1. + A1/nu2 + A2/(nu2*nu2) + A3/(nu2*nu2*nu2);
}

static double olver_Bsum(double nu, double z, double abs_zeta)
{
  double nu2 = nu*nu;
  double B0 = olver_B0(z, abs_zeta);
  double B1 = olver_B1(z, abs_zeta);
  double B2 = olver_B2(z, abs_zeta);
  return B0 + B1/nu2 + B2/(nu2*nu2);
}


int gsl_sf_bessel_Jnu_asymp_Olver_impl(double nu, double x, double * result)
{
  if(x <= 0. || nu <= 0.) {
    return GSL_EDOM;
  }  
  else {
    double zeta, abs_zeta;
    double arg;
    double pre;
    double asum, bsum;
    double ai, aip;
    double z = x/nu;
    double crnu = pow(nu, 1./3.);
    double rt   = sqrt(fabs(1.-z)*(1+z));
    
    if(fabs(1-z) < GSL_SQRT_MACH_EPS) {
      /* z near 1 */
      pre  = CubeRoot2_*(0.5 + 2./5.*(1-z) + 51./175.*(1-z)*(1-z));
      zeta = CubeRoot2_*(1-z)*(1. + 3./10.*(1-z));
      abs_zeta = fabs(zeta);
    }
    else if(z < 1.) {
      /* z < 1 */
      abs_zeta = pow(1.5*(log((1+rt)/z) - rt), 2./3.);
      zeta = abs_zeta;
      pre  = sqrt(sqrt(4.*abs_zeta/(rt*rt)));
    }
    else {
      /* z > 1 */
      abs_zeta = pow(1.5*(rt - acos(1./z)), 2./3.);
      zeta = -abs_zeta;
      pre  = sqrt(sqrt(4.*abs_zeta/(rt*rt)));
    }

    asum = olver_Asum(nu, z, abs_zeta);
    bsum = olver_Bsum(nu, z, abs_zeta);
    arg  = crnu*crnu * zeta;
    gsl_sf_airy_Ai_impl(arg, &ai);
    gsl_sf_airy_Ai_deriv_impl(arg, &aip);

    *result = pre * (ai*asum/crnu + aip*bsum/(nu*crnu*crnu));
    return GSL_SUCCESS;
  }
}

int gsl_sf_bessel_Ynu_asymp_Olver_impl(double nu, double x, double * result)
{
  if(x <= 0. || nu <= 0.) {
    return GSL_EDOM;
  }  
  else {
    double zeta, abs_zeta;
    double arg;
    double pre;
    double asum, bsum;
    double bi, bip;
    double z = x/nu;
    double crnu = pow(nu, 1./3.);
    double rt   = sqrt(fabs(1.-z)*(1+z));
    
    if(fabs(1-z) < GSL_SQRT_MACH_EPS) {
      /* z near 1 */
      pre  = CubeRoot2_*(0.5 + 2./5.*(1-z) + 51./175.*(1-z)*(1-z));
      zeta = CubeRoot2_*(1-z)*(1. + 3./10.*(1-z));
      abs_zeta = fabs(zeta);
    }
    else if(z < 1.) {
      /* z < 1 */
      abs_zeta = pow(1.5*(log((1+rt)/z) - rt), 2./3.);
      zeta = abs_zeta;
      pre  = sqrt(sqrt(4.*abs_zeta/(rt*rt)));
    }
    else {
      /* z > 1 */
      abs_zeta = pow(1.5*(rt - acos(1./z)), 2./3.);
      zeta = -abs_zeta;
      pre  = sqrt(sqrt(4.*abs_zeta/(rt*rt)));
    }

    asum = olver_Asum(nu, z, abs_zeta);
    bsum = olver_Bsum(nu, z, abs_zeta);
    arg  = crnu*crnu * zeta;
    gsl_sf_airy_Bi_impl(arg, &bi);
    gsl_sf_airy_Bi_deriv_impl(arg, &bip);

    *result = -pre * (bi*asum/crnu + bip*bsum/(nu*crnu*crnu));
    return GSL_SUCCESS;
  }
}
