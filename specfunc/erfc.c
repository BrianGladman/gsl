/* $Id$ */
/**
 *  erfc() is the complimentary error function.
 *  This algorithm is from:
 *    Hart et al, Computer Approximations, John Wiley and Sons, New York (1968)
 **/

#include <math.h>

static inline double erfc8_sum(double x)
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
    double num=0,den=0;
    double e;
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

static inline double erfc8(double x)
{
    double e;
    e = erfc8_sum(x);
    e *= exp(-x*x);
    return e;
}
static inline double log_erfc8(double x)
{
    double e;
    e = erfc8_sum(x);
    e = log(e) - x*x;
    return e;
}

static inline double erfc0(double x) {
    /* ERFC 5663, from Hart et al */
    /* this gives >9 places of absolute accuracy for 0<x<10 */
    static double P[] = {
        10.0046411662521,
        8.42655286493187,
        3.460259332062085,
        0.56235361207678689
    };
    static double Q[] = {
        10.0046411714619,
        19.71558073869718,
        15.702288085239204,
        6.0907487871984499,
        1.0
    };
    double num=0,den=0;
    double e;
    int i;

    num = P[3];
    for (i=2; i>=0; --i) {
        num = x*num + P[i];
    }
    den = Q[4];
    for (i=3; i>=0; --i) {
        den = x*den + Q[i];
    }
    
    e = exp(-x*x)*num/den;
    return e;
}

static inline double erfcasympsum(double x)
{
    double e,coef;
    int i;
    e=1;
    coef = 1;
    for (i=1; i<10; ++i) {
        coef *= -(2*i-1)/(2*x*x);
        e += coef;
        if (fabs(coef) < 1.0e-15) break;
        if (fabs(coef) > 1.0e10) break;
    }
    return e;
}

double gsl_erfc_asymptotic(double x)
{
    /* valid as x->infty */
    double e;
    e = erfcasympsum(x);
    e *= exp(-x*x)/x;
    return 0.564189583547756287*e;
}

double gsl_log_erfc_asymptotic(double x)
{
    double e;
    e = erfcasympsum(x);
    e = log(e) - x*x - log(x) - 0.57236494292470008706;
    return e;
}

double gsl_log_erfc(double x)
{
    if (x > 5) {
        return log_erfc8(x);
    } else {
        return log(erfc(x));
    }
}

double gsl_erfc(double x)
{
    double e;
    if (x < 0) {
        if (x < -9) { e = 2.0 - erfc8(-x); }
        else        { e = 2.0 - erfc0(-x); }
    } else {
        if (x>9) { e = erfc8(x); }
        else     { e = erfc0(x); }
    }
    return e;
}



static inline double erfseries(double x)
{
    int k;
    double coef,e;

    coef = x;
    e = coef;
    for (k=1; k<50; ++k) {
        coef *= -x*x/k;
        e += coef/(2*k+1);
        if (fabs(coef) < 1.0e-15) break; /* we are accurate enough */
        if (fabs(coef) > 1.0e10) break; /* we are overflowing! */
    }
    return 2*0.564189583547756287*e;
}

double gsl_erf(double x)
{
    if (fabs(x) < 1.0) 
        return erfseries(x);
    else
        return 1.0-erfc(x);
}


double gsl_Z(double x)
{
    const double oneover_sqrt_twopi = .39894228040143267794;
    exp(-x*x/2.0)*oneover_sqrt_twopi;
}

double gsl_Q(double x)
{
    /* Abramowitz+Stegun, 26.2.17 */
    const double p=0.2316419;
    const double b[] = {
          0.319381530,
         -0.356563782,
          1.781477937,
         -1.821255978,
          1.330274429
    };
    double t,e;
    
    t = 1.0/(1.0+p*fabs(x));
    e = gsl_Z(x)*t*(b[0] + t*(b[1] + t*(b[2] + t*(b[3] + t*b[4]))));
    if (x < 0) e = 1.0-e;
    return e;
}
