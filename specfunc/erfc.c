/* Author:  J. Theiler (modifications by G. Jungman)
 * RCS:     $Id$
 */
/*
 * See Hart et al, Computer Approximations, John Wiley and Sons, New York (1968)
 * (This applies only to the erfc8 stuff, which is the part
 *  of the original code that survives. I have replaced much of
 *  the other stuff with Chebyshev fits. These are simpler and
 *  more precise than the original approximations. [GJ])
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_erf.h"

#define LogRootPi_  0.57236494292470008706


static double erfc8_sum(double x)
{
    /* estimates erfc(x) valid for 8 < x < 100 */
    /* This is based on index 5725 in Hart et al */

    static double P[] = {
        2.97886562639399288862,
        7.409740605964741794425,
        6.1602098531096305440906,
        5.019049726784267463450058,
        1.275366644729965952479585264,
        0.5641895835477550741253201704
    };
    static double Q[] = {
        3.3690752069827527677,
        9.608965327192787870698,
        17.08144074746600431571095,
        12.0489519278551290360340491,
        9.396034016235054150430579648,
        2.260528520767326969591866945,
        1.0
    };
    double num=0.0, den=0.0;
    int i;

    num = P[5];
    for (i=4; i>=0; --i) {
        num = x*num + P[i];
    }
    den = Q[6];
    for (i=5; i>=0; --i) {
        den = x*den + Q[i];
    }

    return num/den;
}

#ifdef HAVE_INLINE
inline
#endif
static double erfc8(double x)
{
    double e;
    e = erfc8_sum(x);
    e *= exp(-x*x);
    return e;
}
#ifdef HAVE_INLINE
inline
#endif
static double log_erfc8(double x)
{
    double e;
    e = erfc8_sum(x);
    e = log(e) - x*x;
    return e;
}


/* Abramowitz+Stegun, 7.2.14 */
static double erfcasympsum(double x)
{
  int i;
  double e = 1.;
  double coef = 1.;
  for (i=1; i<5; ++i) {
    /* coef *= -(2*i-1)/(2*x*x); ??? [GJ] */
    coef *= -(2*i+1)/(i*(4*x*x*x*x));
    e += coef;
    /*
    if (fabs(coef) < 1.0e-15) break;
    if (fabs(coef) > 1.0e10) break;
    
    [GJ]: These tests are not useful. This function is only
    used below. Took them out; they gum up the pipeline.
    */
  }
  return e;
}


/* Abramowitz+Stegun, 7.1.5 */
static double erfseries(double x)
{
  double coef = x;
  double e    = coef;
  int k;
  for (k=1; k<30; ++k) {
    coef *= -x*x/k;
    e += coef/(2.0*k+1.0);
  }
  return 2.0 / M_SQRTPI * e;
}


/* Chebyshev fit for erfc((t+1)/2), -1 < t < 1
 */
static double erfc_xlt1_data[20] = {
  1.06073416421769980345174155056,
 -0.42582445804381043569204735291,
  0.049552626796204340403576830799,
  0.0044929348876838274955800124202,
 -0.00129194104658496953494224761120,
 -0.0000183638929214939627041697886684,
  0.0000221111470409952629153855547035,
 -5.2333748523425713467369317902e-7,
 -2.78184788833537885382530989578e-7,
  1.41158092748813114560316684249e-8,
  2.72571296330561699984539141865e-9,
 -2.06343904872070629406401492476e-10,
 -2.14273991996785367924201401812e-11,
  2.22990255539358204580285098119e-12,
  1.36250074650698280575807934155e-13,
 -1.95144010922293091898995913038e-14,
 -6.8562716923170459944280637069e-16,
  1.44506492869699938239521607493e-16,
  2.45935306460536488037576200030e-18,
 -9.2959956122052339600735932854e-19
};
static struct gsl_sf_cheb_series erfc_xlt1_cs = {
  erfc_xlt1_data,
  19,
  -1, 1,
  (double *) 0,
  (double *) 0
};

/* Chebyshev fit for erfc(x) exp(x^2), 1 < x < 5, x = 2t + 3, -1 < t < 1
 */
static double erfc_x15_data[25] = {
  0.44045832024338111077637466616,
 -0.143958836762168335790826895326,
  0.044786499817939267247056666937,
 -0.0133431242002712112036183531023,
  0.0038246827397504697676923725561,
 -0.00105869922719512654730648253009,
  0.000283859419210073742736310108006,
 -0.000073906170662206760483959432047,
  0.0000187253125214891790158729342482,
 -4.6253098116491944513129726443e-6,
  1.11558657244432857487884006422e-6,
 -2.63098662650834130067808832725e-7,
  6.0746212272455177737211940871e-8,
 -1.37460865539865444777251011793e-8,
  3.05157051905475145520096717210e-9,
 -6.6517478972031071375730772479e-10,
  1.42483346273207784489792999706e-10,
 -3.00141127395323902092018744545e-11,
  6.2217179264534809147291400125e-12,
 -1.26994639225668496876152836555e-12,
  2.55385883033257575402681845385e-13,
 -5.0625823750703869839226549977e-14,
  9.8970540947832732164126422711e-15,
 -1.90685978789192181051961024995e-15,
  3.5082664803273784924511375734e-16
};
static struct gsl_sf_cheb_series erfc_x15_cs = {
  erfc_x15_data,
  24,
  -1, 1,
  (double *) 0,
  (double *) 0
};

/* Chebyshev fit for erfc(x) x exp(x^2), 5 < x < 10, x = (5t + 15)/2, -1 < t < 1
 */
static double erfc_x510_data[20] = {
  1.11684990123545698684297865808,
  0.0037362403593819985206549275366,
 -0.00091662394804547023876361987059,
  0.000199094325044940833965078819585,
 -0.000040276384918650072591781859519,
  7.7651526469706104947712760579e-6,
 -1.44464794206689070402099225301e-6,
  2.61311930343463958393485241947e-7,
 -4.6183302663484415234530409556e-8,
  8.0025311151294360159873214434e-9,
 -1.36291114862793031395712122089e-9,
  2.28570483090160869607683087722e-10,
 -3.7802252156325180504405697456e-11,
  6.1725368387452828572991046213e-12,
 -9.9601929095531688844583059743e-13,
  1.58953143706980770269506726000e-13,
 -2.51045971047162509999527428316e-14,
  3.9260782898912581001358128756e-15,
 -6.0797061938416037439253545342e-16,
  9.1260060726479471731550747767e-17
};
static struct gsl_sf_cheb_series erfc_x510_cs = {
  erfc_x510_data,
  19,
  -1, 1,
  (double *) 0,
  (double *) 0
};



/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_erfc_asymptotic(double x)
{
  return exp(-x*x)/x * erfcasympsum(x) / M_SQRTPI;
}


double gsl_sf_log_erfc_asymptotic(double x)
{
  return log(erfcasympsum(x)/x) - x*x - LogRootPi_;
}


double gsl_sf_erfc(double x)
{
  const double ax = fabs(x);
  double e;

  if(ax <= 1.0) {
    double t = 2.0*ax - 1.0;
    e = gsl_sf_cheb_eval(&erfc_xlt1_cs, t);
  }
  else if(ax <= 5.0) {
    double t = 0.5*(ax-3.0);
    e = gsl_sf_cheb_eval(&erfc_x15_cs, t) * exp(-x*x);
  }
  else if(ax < 10.0) {
    double t = (2.0*x - 15.0)/5.0;
    double c = gsl_sf_cheb_eval(&erfc_x510_cs, t);
    e = c * exp(-x*x) / x;
  }
  else {
    e = erfc8(ax);
  }

  if(x < 0.0)
    return 2.0 - e;
  else
    return e;
}


double gsl_sf_log_erfc(double x)
{
  if (x > 8.0) {
    return log_erfc8(x);
  }
  else {
    return log(gsl_sf_erfc(x));
  }
}


double gsl_sf_erf(double x)
{
  if (fabs(x) < 1.0) 
    return erfseries(x);
  else
    return 1.0-gsl_sf_erfc(x);
}


double gsl_sf_erf_Z(double x)
{
  return exp(-x*x/2.0) / (M_SQRT2 * M_SQRTPI);
}


double gsl_sf_erf_Q(double x)
{
  return 0.5 * gsl_sf_erfc(x/M_SQRT2);
}
