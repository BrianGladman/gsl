#include <config.h>
#include <float.h>
#include <math.h>
#include <gsl_integration.h>

#include <stdio.h>

#include "err.h"

void
gsl_integration_qk (const int n,
		    const double xgk[], const double wg[], 
		    const double wgk[],
		    double fv1[], double fv2[],
		    const gsl_function *f,
		    double a, double b,
		    double * result, double * abserr,
		    double * resabs, double * resasc)
{

  const double center = 0.5 * (a + b);
  const double half_length = 0.5 * (b - a);
  const double abs_half_length = fabs (half_length);
  const double f_center = GSL_FN_EVAL(f, center);

  double result_gauss = 0;
  double result_kronrod = f_center * wgk[n - 1];

  double result_abs = fabs (result_kronrod);
  double result_asc = 0 ;
  double mean = 0, err = 0 ;
  
  int j ;

  if (n % 2 == 0) 
    {
      result_gauss = f_center * wg[n/2 - 1];	
    }

  for (j = 0; j < (n-1)/2 ; j++)
    {
      const int jtw = j * 2 + 1;   /* j=1,2,3 jtw=2,4,6 */
      const double abscissa = half_length * xgk[jtw];
      const double fval1 = GSL_FN_EVAL(f, center - abscissa);
      const double fval2 = GSL_FN_EVAL(f, center + abscissa);
      const double fsum = fval1 + fval2;
      fv1[jtw] = fval1;
      fv2[jtw] = fval2;
      result_gauss +=  wg[j] * fsum;
      result_kronrod +=  wgk[jtw] * fsum;
      result_abs += wgk[jtw] * (fabs (fval1) + fabs (fval2));
    }

  for (j = 0; j < n/2; j++)
    {
      int jtwm1 = j * 2 ;
      const double abscissa = half_length * xgk[jtwm1];
      const double fval1 = GSL_FN_EVAL(f, center - abscissa);
      const double fval2 = GSL_FN_EVAL(f, center + abscissa);
      fv1[jtwm1] = fval1;
      fv2[jtwm1] = fval2;
      result_kronrod += wgk[jtwm1] * (fval1 + fval2);
      result_abs +=  wgk[jtwm1] * (fabs (fval1) + fabs (fval2));
    } ;

  mean = result_kronrod * 0.5;

  result_asc = wgk[n-1] * fabs (f_center - mean);

  for (j = 0; j < n - 1 ; j++)
    {
      result_asc += wgk[j] * (fabs (fv1[j] - mean) + fabs (fv2[j] - mean));
    }

  /* scale by the width of the integration region */

  err = (result_kronrod - result_gauss) * half_length ;

  result_kronrod *= half_length ;
  result_abs *= abs_half_length ;
  result_asc *= abs_half_length ;

  *result = result_kronrod ;
  *resabs = result_abs ;
  *resasc = result_asc ;
  *abserr = rescale_error (err, result_abs, result_asc) ;

}
