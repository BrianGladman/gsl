/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_log.h"


#define locMAX(a,b) ((a) > (b) ? (a) : (b))
#define locMIN(a,b) ((a) < (b) ? (a) : (b))



/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* Chebyshev expansion for log(2 + t/2), -1<t<1
 */
static double lopx_data[30] = {
  -0.138672928390147820418835791217,
   0.53589838486224541294510731699,
  -0.071796769724490825890214633977,
   0.0128252576445603980588699182746,
  -0.00257738807143578123150243783557,
   0.00055248724185826110548585010907,
  -0.000123365425237016197284541806513,
   0.0000283334280567261724856789729930,
  -6.6429292707794557364700505816e-6,
   1.58219336309548624080986748608e-6,
  -3.8155269052018587773948237390e-7,
   9.2942486631641722026407570943e-8,
  -2.28285422158647493139243885126e-8,
   5.6463594933422326752171005708e-9,
  -1.40487050452993623885308455025e-9,
   3.5133832268182697313084600379e-10,
  -8.8257018593156940003341080926e-11,
   2.22573146902729145998374498345e-11,
  -5.6325056360026163263137556112e-12,
   1.42979242408150509933291942482e-12,
  -3.6395613910731953507806644858e-13,
   9.2877860525401087567579002972e-14,
  -2.37553409987900202968605039263e-14,
   6.0884755132754112143765304910e-15,
  -1.56342687634151828778127905972e-15,
   4.0216169888456051529183271296e-16,
  -1.03612365066938659240288862503e-16,
   2.67270435429422491702929511820e-17,
  -6.8766567654198801943924765289e-18,
   1.66708629835556875107467462252e-18
};
static struct gsl_sf_cheb_series lopx_cs = {
  lopx_data,
  29,
  -1, 1,
  (double *)0,
  (double *)0
};


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_log_impl(const double x, double * result)
{
  if(x <= 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else {
    *result = log(x);
    return GSL_SUCCESS;
  }
}

int gsl_sf_log_abs_impl(const double x, double * result)
{
  if(x == 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else {
    *result = log(fabs(x));
    return GSL_SUCCESS;
  }
}

int gsl_sf_complex_log_impl(const double zr, const double zi, double * lnr, double * theta)
{
  if(zr != 0.0 || zi != 0.0) {
    double ax = fabs(zr);
    double ay = fabs(zi);
    double min = locMIN(ax, ay);
    double max = locMAX(ax, ay);
    *lnr   = log(max) + 0.5 * log(1.0 + (min/max)*(min/max));
    *theta = atan2(zi, zr);
    return GSL_SUCCESS;
  }
  else {
    return GSL_EDOM;
  }
}

int gsl_sf_log_1plusx_impl(const double x, double * result)
{
  if(x <= -1.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(fabs(x) < GSL_ROOT4_MACH_EPS) {
    *result = x * (1.0  - x/2.0 + x*x/3.0 - x*x*x/4.0);
    return GSL_SUCCESS;
  }
  else if(fabs(x) <  0.5) {
    double t = 2.0*x;
    *result = gsl_sf_cheb_eval(&lopx_cs, t);
    return GSL_SUCCESS;
  }
  else {
    *result = log(1.0 + x);
    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Error Handling Versions *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_log_e(const double x, double * result)
{
  int status = gsl_sf_log_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_log_e", status);
  }
  return status;
}

int gsl_sf_log_abs_e(const double x, double * result)
{
  int status = gsl_sf_log_abs_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_log_abs_e", status);
  }
  return status;
}

int gsl_sf_complex_log_e(const double zr, const double zi, double * lnr, double * theta)
{
  int status = gsl_sf_complex_log_impl(zr, zi, lnr, theta);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_complex_log_e", status);
  }
  return status;
}

int gsl_sf_log_1plusx_e(const double x, double * result)
{
  int status = gsl_sf_log_1plusx_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_log_1plusx_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_log(const double x)
{
  double y;
  int status = gsl_sf_log_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_log", status);
  }
  return y;
}

double gsl_sf_log_abs(const double x)
{
  double y;
  int status = gsl_sf_log_abs_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_log_abs", status);
  }
  return y;
}

double gsl_sf_log_1plusx(const double x)
{
  double y;
  int status = gsl_sf_log_1plusx_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_log_1plusx", status);
  }
  return y;
}
