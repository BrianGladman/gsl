/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_synchrotron.h"


/* based on SYNCH1(), SYNCH2() of MISCFUN, Allan J. MacLeod  [TOMS 757] */

/* SYNCH1(x) = x * Integral{x to inf} K(5/3)(t) dt */

/* SYNCH2(x) = x * K(2/3)(x)                       */


static double async1_data[14] = {
  30.36468298250107627340,
  17.07939527740839457449,
   4.56013213354507288887,
   0.54928124673041997963,
   0.3729760750693011724e-01,
   0.161362430201041242e-02,
   0.4819167721203707e-04,
   0.105124252889384e-05,
   0.1746385046697e-07,
   0.22815486544e-09,
   0.240443082e-11,
   0.2086588e-13,
   0.15167e-15,
   0.94e-18
};
static struct gsl_sf_cheb_series async1_cs = {
  async1_data,
  13,
  -1.0, 1.0,
  (double *)0,
  (double *)0
};

static double async2_data[12] = {
  0.44907216235326608443,
  0.8983536779941872179e-01,
  0.810445737721512894e-02,
  0.42617169910891619e-03,
  0.1476096312707460e-04,
  0.36286336153998e-06,
  0.666348074984e-08,
  0.9490771655e-10,
  0.107912491e-11,
  0.1002201e-13,
  0.7745e-16,
  0.51e-18
};
static struct gsl_sf_cheb_series async2_cs = {
  async2_data,
  11,
  -1.0, 1.0,
  (double *)0,
  (double *)0
};

static double asynca_data[25] = {
  2.13293051613550009848,
  0.7413528649542002401e-01,
  0.869680999099641978e-02,
  0.117038262487756921e-02,
  0.16451057986191915e-03,
  0.2402010214206403e-04,
  0.358277563893885e-05,
  0.54477476269837e-06,
  0.8388028561957e-07,
  0.1306988268416e-07,
  0.205309907144e-08,
  0.32518753688e-09,
  0.5179140412e-10,
  0.830029881e-11,
  0.133527277e-11,
  0.21591498e-12,
  0.3499673e-13,
  0.569942e-14,
  0.92906e-15,
  0.15222e-15,
  0.2491e-16,
  0.411e-17,
  0.67e-18,
  0.11e-18,
  0.2e-19,
};
static struct gsl_sf_cheb_series asynca_cs = {
  asynca_data,
  24,
  -1.0, 1.0,
  (double *)0,
  (double *)0
};

static double asyn21_data[15] = {
  38.61783992384308548014,
  23.03771559496373459697,
  5.38024998683357059676,
  0.61567938069957107760,
  0.4066880046688955843e-01,
  0.172962745526484141e-02,
  0.5106125883657699e-04,
  0.110459595022012e-05,
  0.1823553020649e-07,
  0.23707698034e-09,
  0.248872963e-11,
  0.2152868e-13,
  0.15607e-15,
  0.96e-18,
  0.1e-19
};
static struct gsl_sf_cheb_series asyn21_cs = {
  asyn21_data,
  14,
  -1.0, 1.0,
  (double *)0,
  (double *)0
};

static double asyn22_data[14] = {
   7.90631482706608042875,
   3.13534636128534256841,
   0.48548794774537145380,
   0.3948166758272372337e-01,
   0.196616223348088022e-02,
   0.6590789322930420e-04,
   0.158575613498559e-05,
   0.2868653011233e-07,
   0.40412023595e-09,
   0.455684443e-11,
   0.4204590e-13,
   0.32326e-15,
   0.210e-17,
   0.1e-19
};
static struct gsl_sf_cheb_series asyn22_cs = {
  asyn22_data,
  13,
  -1.0, 1.0,
  (double *)0,
  (double *)0
};

static double asyn2a_data[19] = {
  2.02033709417071360032,
  0.1095623712180740443e-01,
  0.85423847301146755e-03,
  0.7234302421328222e-04,
  0.631244279626992e-05,
  0.56481931411744e-06,
  0.5128324801375e-07,
  0.471965329145e-08,
  0.43807442143e-09,
  0.4102681493e-10,
  0.386230721e-11,
  0.36613228e-12,
  0.3480232e-13,
  0.333010e-14,
  0.3185e-15,
  0.3074e-16,
  0.295e-17,
  0.29e-18,
  0.3e-19
};
static struct gsl_sf_cheb_series asyn2a_cs = {
  asyn2a_data,
  18,
  -1.0, 1.0,
  (double *)0,
  (double *)0
};


int gsl_sf_synchrotron_1_impl(const double x, double * result)
{
  const double CONLOW = 2.14952824153447863671;
  const double PIBRT3 = 1.81379936423421785059;
  const double LNRTP2 = 0.22579135264472743236;
  const double xlow = 2.0*M_SQRT2 * GSL_SQRT_MACH_EPS;
  const double xhi1 = -8.0*GSL_LOG_DBL_MIN / 7.0;

  if(x < 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x < xlow) {
    *result = CONLOW * pow(x, 1.0/3.0);
    return GSL_SUCCESS;
  }
  else if(x <= 4.0) {
    double t = x*x/8.0 - 1.0;
    double cheb1 = gsl_sf_cheb_eval(&async1_cs, t);
    double cheb2 = gsl_sf_cheb_eval(&async2_cs, t);
    *result = pow(x,1.0/3.0) * cheb1 - pow(x,11.0/3.0) * cheb2 - PIBRT3 * x;
    return GSL_SUCCESS;
  }
  else if(x < xhi1) {
    double t = (12.0 - x) / (x + 4.0);
    double cheb1 = gsl_sf_cheb_eval(&asynca_cs, t);
    double y = LNRTP2 - x + log(sqrt(x) * cheb1);
    if(y < GSL_LOG_DBL_MIN) {
      *result = 0.0;
      return GSL_EUNDRFLW;
    }
    else {
      *result = exp(y);
      return GSL_SUCCESS;
    }
  }
  else {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }
}

int gsl_sf_synchrotron_2_impl(const double x, double * result)
{
  const double CONLOW = 1.07476412076723931836;
  const double LNRTP2 = 0.22579135264472743236;
  const double xlow = 2.0*M_SQRT2*GSL_SQRT_MACH_EPS;
  const double xhi1 = -8.0*GSL_LOG_DBL_MIN / 7.0;
  
  if(x < 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x < xlow) {
    *result = CONLOW * pow(x, 1.0/3.0);
    return GSL_SUCCESS;
  }
  else if(x <= 4.0) {
    double t = x*x/8.0 - 1.0;
    double cheb1 = gsl_sf_cheb_eval(&asyn21_cs, t);
    double cheb2 = gsl_sf_cheb_eval(&asyn22_cs, t);
    *result = pow(x, 1.0/3.0) * cheb1 - pow(x, 5.0/3.0) * cheb2;
    return GSL_SUCCESS;
  }
  else if(x < xhi1) {
    double t = (10.0 - x) / (x + 2.0);
    double cheb1 = gsl_sf_cheb_eval(&asyn2a_cs, t);
    double y = LNRTP2 - x + log(sqrt(x) * cheb1);
    if(y < GSL_LOG_DBL_MIN) {
      *result = 0.0;
      return GSL_EUNDRFLW;
    }
    else {
      *result = exp(y);
      return GSL_SUCCESS;
    }
  }
  else {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }
};
