/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_bessel.h"


/* correction terms to leading McMahon expansion
 * [Abramowitz+Stegun 9.5.12]
 */
static int
mcmahon_correction(const double mu, const double beta)
{
  const double eb   = 8.0*beta;
  const double ebsq = eb*eb;
  const double n2 = 4.0/3.0    * (7.0*mu - 31.0);
  const double n3 = 32.0/15.0  * (mu*(83.0*mu - 982.0) + 3779.0);
  const double n4 = 64.0/105.0 * (mu*(mu*(6949.0*mu - 153855.0) + 1585743.0) - 6277237.0);
  const double term1 = -(mu-1.0) / eb;
  const double term2 =  term1 * n2 / (ebsq);
  const double term3 =  term1 * n3 / (ebsq*ebsq);
  const double term4 =  term1 * n4 / (ebsq*ebsq*ebsq);
}


int
gsl_sf_bessel_zero_J0_impl(int s, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else {
    /* See [F. Lether, J. Comp. Appl .Math. 67, 167 (1996)]. */

    const static double P[] = { 1567450796.0/12539606369.0,
                                8903660.0/2365861.0,
                                10747040.0/536751.0,
                                17590991.0/1696654.0
                              };
    const static double Q[] = { 1.0,
                                29354255.0/954518.0,
                                76900001.0/431847.0,
                                67237052.0/442411.0
                              };

    const double beta = (s - 0.25) * M_PI;
    const double bi2  = 1.0/(beta*beta);
    const double R33num = P[0] + bi2 * (P[1] + bi2 * (P[2] + P[3] * bi2));
    const double R33den = Q[0] + bi2 * (Q[1] + bi2 * (Q[2] + Q[3] * bi2));
    const double R33 = R33num/R33den;
    result->val = beta + R33/beta;
    result->err = fabs(3.0e-15 * result->val);
    return GSL_SUCCESS;
  }
}


int
gsl_sf_bessel_zero_J0_e(int s, gsl_sf_result * result)
{
  int status = gsl_sf_bessel_zero_J0_impl(s, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_zero_e", status);
  }
  return status;
}
