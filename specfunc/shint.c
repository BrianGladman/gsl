#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_expint.h"


extern int gsl_sf_expint_Ei_impl(double, double *);
extern int gsl_sf_expint_E1_impl(double, double *);



/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

/* based on SLATEC shi.f, W. Fullerton

 series for shi  on the interval  0.00000e+00 to  1.40625e-01
					with weighted error   4.67e-20
					 log weighted error  19.33
			       significant figures required  17.07
				    decimal places required  19.75
*/
static double shi_data[7] = {
   0.0078372685688900950695,
   0.0039227664934234563973,
   0.0000041346787887617267,
   0.0000000024707480372883,
   0.0000000000009379295591,
   0.0000000000000002451817,
   0.0000000000000000000467
};
static struct gsl_sf_ChebSeries shi_cs = {
  shi_data,
  6,
  -1, 1
};

int gsl_sf_Shi_impl(const double x, double * result)
{
  const double xsml = GSL_SQRT_MACH_EPS;  /* sqrt (r1mach(3)) */
  const double ax = fabs(x);

  if(x < xsml) {
    *result = x;
    return GSL_SUCCESS;
  }
  else if(ax <= 0.375) {
    *result = x * (1.0 + gsl_sf_cheb_eval(128.*x*x/9.-1., &shi_cs));
    return GSL_SUCCESS;
  }
  else {
    double Ei, E1;
    int status_Ei = gsl_sf_expint_Ei_impl(x, &Ei);
    int status_E1 = gsl_sf_expint_E1_impl(x, &E1);
    *result = 0.5*(Ei + E1);
    if(status_Ei == GSL_EUNDRFLW && status_E1 == GSL_EUNDRFLW) {
      return GSL_EUNDRFLW;
    }
    else {
      return GSL_SUCCESS;
    }
  }
}

int gsl_sf_Chi_impl(const double x, double * result)
{
  double Ei, E1;
  int status_Ei = gsl_sf_expint_Ei_impl(x, &Ei);
  int status_E1 = gsl_sf_expint_Ei_impl(x, &E1);
  if(status_Ei == GSL_EDOM || status_E1 == GSL_EDOM) {
    return GSL_EDOM;
  }
  else if(status_Ei == GSL_EUNDRFLW && status_E1 == GSL_EUNDRFLW) {
    *result = 0.;
    return GSL_EUNDRFLW;
  }
  else {
    *result = 0.5 * (Ei - E1);
    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_Shi_e(double x, double * result)
{
  int status = gsl_sf_Shi_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_Shi_e", status);
  }
  return status;
}

int gsl_sf_Chi_e(double x, double * result)
{
  int status = gsl_sf_Chi_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_Chi_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_Shi(double x)
{
  double y;
  int status = gsl_sf_Shi_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_Shi");
  }
  return y;
}

double gsl_sf_Chi(double x)
{
  double y;
  int status = gsl_sf_Chi_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_Chi");
  }
  return y;
}
