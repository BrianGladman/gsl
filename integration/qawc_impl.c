#include <config.h>
#include <math.h>
#include <float.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include <gsl_integration.h>

#include "qpsrt.h"

int
gsl_integration_qawc_impl (const gsl_function *f,
			  const double a, const double b, const double c,
			  const double epsabs, const double epsrel,
			  gsl_integration_workspace * workspace,
			  double * result, double * abserr)
{
  double result0, abserr0, resabs0, resasc0;
  double tolerance ;
  double e_max, area, errsum;
  size_t i_max, nrmax, i;
  int roundoff_type1 = 0, roundoff_type2 = 0, error_type = 0;

  volatile double round_off ; /* volatile is needed for strict IEEE behavior */

  const size_t limit = workspace->limit ;
  double * alist = workspace->alist ;
  double * blist = workspace->blist ;
  double * rlist = workspace->rlist ;
  double * elist = workspace->elist ;
  size_t * order = workspace->order ;

  alist[0] = a;
  blist[0] = b;
  rlist[0] = 0;
  elist[0] = 0;
  order[0] = 0;
  workspace->size = 0;

  if (epsabs <= 0 && (epsrel < 50 * GSL_DBL_EPSILON || epsrel < 0.5e-28))
    {
      *result = 0;
      *abserr = 0;
      *nqeval = 0;
      GSL_ERROR ("tolerance cannot be acheived with given epsabs and epsrel",
		 GSL_EBADTOL);
    };

  /* perform the first integration */

  q (f, a, b, &result0, &abserr0, &resabs0, &resasc0);

  rlist[0] = result0;
  elist[0] = abserr0;
  order[0] = 0;
  workspace->size = 1;

  /* Test on accuracy */

  tolerance = GSL_MAX_DBL (epsabs, epsrel * fabs (result0));

  round_off = 50 * GSL_DBL_EPSILON * resabs0 ;

  if (abserr0 <= round_off && abserr0 > tolerance)
    {
      *result = result0;
      *abserr = abserr0;

      GSL_ERROR ("cannot reach tolerance because of roundoff error "
		 "on first attempt", GSL_EROUND);
    }
  else if ((abserr0 <= tolerance && abserr0 != resasc0) || abserr0 == 0)
    {
      *result = result0;
      *abserr = abserr0;

      return GSL_SUCCESS;
    }
  else if (limit == 1)
    {
      *result = result0;
      *abserr = abserr0;

      GSL_ERROR ("a maximum of one iteration was insufficient", GSL_EMAXITER);
    }

  e_max = abserr0;
  i_max = 0;
  area = result0;
  errsum = abserr0;
  nrmax = 0;

  i = 1;

  do
    {
      /* Bisect the subinterval with the largest error estimate */

      const double left = alist[i_max];
      const double right = blist[i_max];
      const double midpoint = 0.5 * (left + right);

      const double a1 = left, b1 = midpoint;
      const double a2 = midpoint, b2 = right;

      double area1 = 0, area2 = 0, area12 = 0;
      double error1 = 0, error2 = 0, error12 = 0;
      double resasc1, resasc2;
      double resabs1, resabs2;

      if (c > a1 && c < b1) 
	{
	  b1 = 0.5 * (c + b2) ;
	  a2 = b1;
	}
      else if (c > b1 && c < b2)
	{
	  b1 = 0.5 * (a1 + c) ;
	  a2 = b1;
	}
      
      q (f, a1, b1, c, &area1, &error1, neval);
      q (f, a2, b2, c, &area2, &error2, neval);

      area12 = area1 + area2;
      error12 = error1 + error2;

      errsum += (error12 - e_max);
      area += area12 - rlist[i_max];

      if (resasc1 != error1 && resasc2 != error2)
	{
	  if (fabs (rlist[i_max] - area12) <= 0.00001 * fabs (area12)
	      && error12 >= 0.99 * e_max)
	    {
	      roundoff_type1++;
	    }
	  if (i >= 10 && error12 > e_max)
	    {
	      roundoff_type2++;
	    }
	}

      tolerance = GSL_MAX_DBL (epsabs, epsrel * fabs (area));

      if (errsum > tolerance)
	{
	  if (roundoff_type1 >= 6 || roundoff_type2 >= 20)
	    {
	      error_type = 2;	/* round off error */
	    }

	  /* set error flag in the case of bad integrand behaviour at
	     a point of the integration range */

	  {
	    volatile double tmp = ((1 + 100 * GSL_DBL_EPSILON) 
				   * (fabs (a2) + 1000 * GSL_DBL_MIN));
	    if (fabs (a1) <= tmp && fabs (b2) <= tmp)
	      {
		error_type = 3;
	      }
	  }
	}

      /* append the newly-created intervals to the list */

      if (error2 > error1)
	{
	  alist[i_max] = a2;	/* blist[maxerr] is already == b2 */
	  rlist[i_max] = area2;
	  elist[i_max] = error2;

	  alist[i] = a1;
	  blist[i] = b1;
	  rlist[i] = area1;
	  elist[i] = error1;
	  workspace->size = i + 1;
	}
      else
	{
	  blist[i_max] = b1;	/* alist[maxerr] is already == a1 */
	  rlist[i_max] = area1;
	  elist[i_max] = error1;

	  alist[i] = a2;
	  blist[i] = b2;
	  rlist[i] = area2;
	  elist[i] = error2;
	  workspace->size = i + 1;
	}

      /* call subroutine dqpsrt to maintain the descending ordering in
         the list of error estimates and select the subinterval with the
         largest error estimate (to be bisected next) */

      qpsrt (limit, i, &i_max, &e_max, elist, order, &nrmax);

      i++;

    }
  while (i < limit && !error_type && errsum > tolerance);

  {
    double result_sum = 0;
    size_t k;
    for (k = 0; k < i; k++)
      {
	result_sum += rlist[k];
      }
    *result = result_sum;
  }

  *abserr = errsum;

  if (errsum <= tolerance)
    {
      return GSL_SUCCESS;
    }
  else if (error_type == 2)
    {
      GSL_ERROR ("roundoff error prevents tolerance from being achieved",
		 GSL_EROUND);
    }
  else if (error_type == 3)
    {
      GSL_ERROR ("bad integrand behavior found in the integration interval",
		 GSL_ESING);
    }
  else if (i == limit)
    {
      GSL_ERROR ("maximum number of subdivisions reached", GSL_EMAXITER);
    }
  
  /* FIXME: we get here if there was a NAN in the function evaluations */

  GSL_ERROR ("shouldn't happen", GSL_ESANITY);

}
