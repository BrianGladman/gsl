#include <float.h>
#include <math.h>
#include <gsl_integration.h>

#include <stdio.h>

int
gsl_integration_qk (const int n,
			 const double xgk[], const double wg[], 
			 const double wgk[],
			 double fv1[], double fv2[],
			 double (*f) (double x),
			 double a, double b,
			 double * result, double * abserr,
			 double * resabs, double * resasc)
{

  const double center = 0.5 * (a + b);
  const double half_length = 0.5 * (b - a);
  const double fc = (*f) (center);

  double result_gauss = 0;
  double result_kronrod = fc * wgk[n - 1];

  double result_abs = fabs (result_kronrod);
  double result_asc = 0 ;
  double mean = 0 ;
  
  int j ;

  if (n % 2 == 0) 
    {
      result_gauss = fc * wg[n/2 - 1];	
    }

  for (j = 0; j < (n-1)/2 ; j++)
    {
      const int jtw = j * 2 + 1;   /* j=1,2,3 jtw=2,4,6 */
      const double abscissa = half_length * xgk[jtw];
      const double fval1 = (*f) (center - abscissa);
      const double fval2 = (*f) (center + abscissa);
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
      const double fval1 = (*f)(center - abscissa);
      const double fval2 = (*f)(center + abscissa);
      fv1[jtwm1] = fval1;
      fv2[jtwm1] = fval2;
      result_kronrod += wgk[jtwm1] * (fval1 + fval2);
      result_abs +=  wgk[jtwm1] * (fabs (fval1) + fabs (fval2));
    } ;

  mean = result_kronrod * 0.5;

  result_asc = wgk[n-1] * fabs (fc - mean);

  for (j = 0; j < n - 1 ; j++)
    {
      result_asc += wgk[j] * (fabs (fv1[j] - mean) + fabs (fv2[j] - mean));
    }

  *result = result_kronrod * half_length;

  {
      const double abs_half_length = fabs (half_length);
      result_abs *= abs_half_length;
      result_asc *= abs_half_length;
  }

  *resabs = result_abs ;
  *resasc = result_asc ;

  {
    double err = fabs ((result_kronrod - result_gauss) * half_length);

    printf("result_gauss   = %.18g\n",result_gauss) ;
    printf("result_kronrod = %.18g\n",result_kronrod) ;
    printf("err = %g\n",err) ;

    if (result_asc != 0 && err != 0)
      {
	double scale = pow((200 * err / result_asc), 1.5) ;

	if (scale < 1)
	  {
	    err = result_asc * scale ;
	  }
	else 
	  {
	    err = result_asc ;
	  }
      }
    if (result_abs > DBL_MIN / (50 * DBL_EPSILON))
      {
	double min_err = (DBL_EPSILON * 50) * result_abs ;
	if (min_err > err) 
	  {
	    printf("min_err = %g, err = %g\n",min_err,err) ;
	    err = min_err ;
	  }
      }
    
    printf("err = %g\n",err) ;

    *abserr = err ;
  }
  return 0 ;
}
