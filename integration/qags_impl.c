#include <config.h>
#include <math.h>
#include <float.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include <gsl_integration.h>

#include "qpsrt.h"

int
gsl_integration_qags_impl (const gsl_function * f,
			   const double a, const double b,
			   const double epsabs, const double epsrel,
			   gsl_integration_workspace * workspace,
			   double *result, double *abserr,
			   size_t * last, size_t * nqeval,
			   gsl_integration_rule_t * const q)
{
  double q_result, q_abserr, q_defabs, q_resabs;
  double tolerance, maxerr_value, area, errsum, last_maxerr_value;
  double small = 0, ertest = 0;
  double error_over_large_intervals = 0;
  double reseps = 0, abseps = 0, correc = 0;
  size_t maxerr_index, nrmax = 0, i = 0, nres = 0, numrl2 = 1, ktmin = 0;
  int roundoff_type1 = 0, roundoff_type2 = 0, roundoff_type3 = 0;
  int error_type = 0, error_type2 = 0;

  const size_t limit = workspace->limit;
  double *alist = workspace->alist;
  double *blist = workspace->blist;
  double *rlist = workspace->rlist;
  double *elist = workspace->elist;
  size_t *iord = workspace->iord;

  int positive_integrand = 0;
  int extrapolate = 0;
  int disallow_extrapolation = 0;

  double res3la[3], rlist2[52];

  alist[0] = a;
  blist[0] = b;
  rlist[0] = 0;
  elist[0] = 0;
  iord[0] = 0;

  /* Test on accuracy */

  if (epsabs <= 0 && (epsrel < 50 * GSL_DBL_EPSILON || epsrel < 0.5e-28))
    {
      *result = 0;
      *abserr = 0;
      *nqeval = 0;
      GSL_ERROR ("tolerance cannot be acheived with given epsabs and epsrel",
		 GSL_EBADTOL);
    }

  /* Perform the first integration */

  q (f, a, b, &q_result, &q_abserr, &q_defabs, &q_resabs);

  rlist[0] = q_result;
  elist[0] = q_abserr;
  iord[0] = 0;

  tolerance = GSL_MAX_DBL (epsabs, epsrel * fabs (q_result));

  if (q_abserr <= 100 * GSL_DBL_EPSILON * q_defabs && q_abserr > tolerance)
    {
      *result = q_result;
      *abserr = q_abserr;
      *nqeval = 1;
      *last = 0;

      GSL_ERROR ("cannot reach tolerance because of roundoff error"
		 "on first attempt", GSL_EROUND);
    }
  else if ((q_abserr <= tolerance && q_abserr != q_resabs) || q_abserr == 0)
    {
      *result = q_result;
      *abserr = q_abserr;
      *nqeval = 1;
      *last = 0;
      return GSL_SUCCESS;
    }
  else if (limit == 1)
    {
      *result = q_result;
      *abserr = q_abserr;
      *nqeval = 1;
      *last = 0;
      GSL_ERROR ("a maximum of one iteration was insufficient", GSL_EMAXITER);
    }

  /* Initialization */

  rlist2[0] = q_result;
  maxerr_value = q_abserr;
  maxerr_index = 0;
  area = q_result;
  errsum = q_abserr;
  *abserr = GSL_DBL_MAX;

  /* Compare the integral of f(x) with the integral of |f(x)|
     to determine if f(x) covers both positive and negative values */

  if (fabs (q_result) >= (1 - 50 * GSL_DBL_EPSILON) * q_defabs)
    {
      positive_integrand = 1;
    }
  else
    {
      positive_integrand = 0;
    }

  do
    {
      /* Bisect the subinterval with the nrmax-th largest error estimate */

      const double left = alist[maxerr_index];
      const double right = blist[maxerr_index];
      const double midpoint = 0.5 * (left + right);

      const double a1 = left, b1 = midpoint;
      const double a2 = midpoint, b2 = right;

      double area1 = 0, area2 = 0, area12 = 0;
      double error1 = 0, error2 = 0, error12 = 0;
      double defab1, defab2;
      double resabs = 0;

      i++;

      q (f, a1, b1, &area1, &error1, &resabs, &defab1);
      q (f, a2, b2, &area2, &error2, &resabs, &defab2);

      area12 = area1 + area2;
      error12 = error1 + error2;
      last_maxerr_value = maxerr_value;

      /* Improve previous approximations to the integral and test for
         accuracy.

         We write these expressions in the same way as the original
         QUADPACK code so that the rounding errors are the same, which
         makes testing easier. */

      errsum = errsum + error12 - maxerr_value;
      area = area + area12 - rlist[maxerr_index];

      tolerance = GSL_MAX_DBL (epsabs, epsrel * fabs (area));

      if (defab1 != error1 && defab2 != error2)
	{
	  if (fabs (rlist[maxerr_index] - area12) <= 1e-5 * fabs (area12)
	      && error12 >= 0.99 * maxerr_value)
	    {
	      if (!extrapolate)
		{
		  roundoff_type1++;
		}
	      else
		{
		  roundoff_type2++;
		}
	    }
	  if (i > 10 && error12 > maxerr_value)
	    {
	      roundoff_type3++;
	    }
	}

      /* Test for roundoff and eventually set error flag */

      if (roundoff_type1 + roundoff_type2 >= 10 || roundoff_type3 >= 20)
	{
	  error_type = 2;	/* round off error */
	}

      if (roundoff_type2 >= 5)
	{
	  error_type2 = 1;
	}

      /* set error flag in the case of bad integrand behaviour at
         a point of the integration range */

      {
	double tmp = ((1 + 100 * GSL_DBL_EPSILON)
		      * (fabs (a2) + 1000 * GSL_DBL_MIN));
	if (fabs (a1) <= tmp && fabs (b2) <= tmp)
	  {
	    error_type = 4;
	  }
      }

      /* append the newly-created intervals to the list */

      if (error2 > error1)
	{
	  alist[maxerr_index] = a2;	/* already done blist[maxerr] = b2 */
	  rlist[maxerr_index] = area2;
	  elist[maxerr_index] = error2;

	  alist[i] = a1;
	  blist[i] = b1;
	  rlist[i] = area1;
	  elist[i] = error1;
	}
      else
	{
	  blist[maxerr_index] = b1;	/* alist[maxerr] is already == a1 */
	  rlist[maxerr_index] = area1;
	  elist[maxerr_index] = error1;

	  alist[i] = a2;
	  blist[i] = b2;
	  rlist[i] = area2;
	  elist[i] = error2;
	}

      /* call subroutine dqpsrt to maintain the descending ordering in
         the list of error estimates and select the subinterval with
         the nrmax-th largest error estimate (to be bisected next) */

      qpsrt (limit, i, &maxerr_index, &maxerr_value, elist, iord, &nrmax);

      if (errsum <= tolerance)
	{
	  goto compute_result;
	}

      if (error_type)
	{
	  break;
	}

      if (i >= limit - 1)
	{
	  error_type = 1;
	  break;
	}

      if (i == 1)		/* set up variables on first iteration */
	{
	  small = fabs (b - a) * 0.375;
	  error_over_large_intervals = errsum;
	  ertest = tolerance;
	  rlist2[1] = area;
	  continue;
	}

      if (disallow_extrapolation)
	{
	  continue;
	}

      error_over_large_intervals += -last_maxerr_value;

      if (fabs (b1 - a1) > small)
	{
	  error_over_large_intervals += error12;
	}

      if (!extrapolate)
	{
	  /* test whether the interval to be bisected next is the
	     smallest interval. */
	  if (fabs (blist[maxerr_index] - alist[maxerr_index]) > small)
	    {
	      continue;
	    }
	  extrapolate = 1;
	  nrmax = 1;
	}

      if (!error_type2 && error_over_large_intervals > ertest)
	{
	  /* The smallest interval has the largest error.  Before
	     bisecting decrease the sum of the errors over the larger
	     intervals (error_over_large_intervals) and perform
	     extrapolation. */

	  int k, flag = 0;
	  int id = nrmax;
	  int jupbnd;
	  if (i > (1 + limit / 2))
	    {
	      jupbnd = limit + 1 - i;
	    }
	  else
	    {
	      jupbnd = i;
	    }

	  for (k = id; k <= jupbnd && !flag; k++)
	    {

	      maxerr_index = iord[nrmax];
	      maxerr_value = elist[maxerr_index];
	      if (fabs (blist[maxerr_index] - alist[maxerr_index]) > small)
		{
		  flag = 1;
		  break;
		}
	      nrmax++;
	    }
	  if (flag)
	    continue;
	}

      /* Perform extrapolation */

      numrl2++;
      rlist2[numrl2] = area;

      gsl_integration_qelg (&numrl2, rlist2, &reseps, &abseps, res3la, &nres);

      ktmin++;

      if (ktmin > 5 && *abserr < 0.001 * errsum)
	{
	  error_type = 5;
	}

      if (abseps < *abserr)
	{
	  ktmin = 0;
	  *abserr = abseps;
	  *result = reseps;
	  correc = error_over_large_intervals;
	  ertest = GSL_MAX_DBL (epsabs, epsrel * fabs (reseps));
	  if (*abserr <= ertest)
	    break;
	}

      /* Prepare bisection of the smallest interval. */

      if (numrl2 == 0)
	{
	  disallow_extrapolation = 1;
	}

      if (error_type == 5)
	{
	  break;
	}

      maxerr_index = iord[0];
      maxerr_value = elist[maxerr_index];
      nrmax = 0;
      extrapolate = 0;
      small *= 0.5;
      error_over_large_intervals = errsum;

    }
  while (i < limit);

  if (*abserr == GSL_DBL_MAX)
    goto compute_result;

  if (error_type || error_type2)
    {
      if (error_type2)
	{
	  *abserr += correc;
	}

      if (error_type == 0)
	error_type = 3;

      if (result != 0 && area != 0)
	{
	  if (*abserr / fabs (*result) > errsum / fabs (area))
	    goto compute_result;
	}
      else if (*abserr > errsum)
	{
	  goto compute_result;
	}
      else if (area == 0)
	{
	  goto return_error;
	}
    }

  /*  Test on divergence. */

  {
    double max_area = GSL_MAX_DBL (fabs (*result), fabs (area));

    if (!positive_integrand && max_area < 0.01 * q_defabs)
      goto return_error;
  }

  {
    double ratio = *result / area;

    if (ratio < 0.01 || ratio > 100 || errsum > fabs (area))
      error_type = 6;
  }

  goto return_error;

compute_result:

  {
    /* Compute global integral sum. */

    double result_sum = 0;
    size_t k;
    for (k = 0; k <= i; k++)
      {
	result_sum += rlist[k];
      }
    *result = result_sum;
  }

  *abserr = errsum;

return_error:

  *nqeval = 2 * i + 1;
  *last = i + 1;

  if (error_type > 2)
    error_type--;

  if (error_type == 1)
    {
      GSL_ERROR ("number of iterations was insufficient", GSL_EMAXITER);
    }
  else if (error_type == 2)
    {
      GSL_ERROR ("cannot reach tolerance because of roundoff error",
		 GSL_EROUND);
    }
  else if (error_type == 3)
    {
      GSL_ERROR ("bad integrand behavior found in the integration interval",
		 GSL_ESING);
    }
  else if (error_type == 4)
    {
      GSL_ERROR ("roundoff error detected in the extrapolation table",
		 GSL_EROUND);
    }
  else if (error_type == 5)
    {
      GSL_ERROR ("integral is divergent, or slowly convergent",
		 GSL_EDIVERGE);
    }

  return GSL_SUCCESS;

}
