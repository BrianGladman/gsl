#include <config.h>
#include <math.h>
#include <float.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include <gsl_integration.h>

#include "qpsrt.h"

int
gsl_integration_qagp_impl (const gsl_function * f,
			   const double *pts, const size_t npts,
			   const double epsabs, const double epsrel,
			   gsl_integration_workspace * workspace,
			   gsl_integration_workspace_pts * workspace_pts,
			   double *result, double *abserr,
			   size_t * last, size_t * nqeval,
			   gsl_integration_rule_t * const q)
{
  double q_result, q_abserr, q_resabs;
  double tolerance, maxerr_value, area, errsum, last_maxerr_value;
  double ertest = 0;
  double error_over_large_intervals = 0;
  double reseps = 0, abseps = 0, correc = 0;
  size_t maxerr_index, nrmax = 0, i = 0, nres = 0, numrl2, ktmin = 0;
  size_t maximum_level = 1;
  int roundoff_type1 = 0, roundoff_type2 = 0, roundoff_type3 = 0;
  int error_type = 0, error_type2 = 0;

  const size_t nint = npts - 1;	/* number of intervals */

  const size_t limit = workspace->limit;
  double *alist = workspace->alist;
  double *blist = workspace->blist;
  double *rlist = workspace->rlist;
  double *elist = workspace->elist;
  size_t *iord = workspace->iord;

  unsigned int *level = workspace_pts->level;
  unsigned int *ndin = workspace_pts->ndin;

  int positive_integrand = 0;
  int extrapolate = 0;
  int disallow_extrapolation = 0;

  double res3la[3], rlist2[52];

  /* Test on validity of parameters */

  if (npts > workspace_pts->npts)
    {
      *result = 0;
      *abserr = 0;
      *nqeval = 0;
      GSL_ERROR ("npts exceeds size of workspace_pts", GSL_EINVAL);
    }

  if (epsabs <= 0 && (epsrel < 50 * GSL_DBL_EPSILON || epsrel < 0.5e-28))
    {
      *result = 0;
      *abserr = 0;
      *nqeval = 0;
      GSL_ERROR ("tolerance cannot be acheived with given epsabs and epsrel",
		 GSL_EBADTOL);
    }

  /* Check that the integration range and break points are an
     ascending sequence */

  for (i = 0; i < nint; i++)
    {
      if (pts[i + 1] < pts[i])
	{
	  GSL_ERROR ("points are not in an ascending sequence", GSL_EINVAL);
	}
    }

  /* Perform the first integration */

  q_result = 0;
  q_abserr = 0;
  q_resabs = 0;

  for (i = 0; i < nint; i++)
    {
      double area1, error1, resabs1, resasc1;
      const double a1 = pts[i];
      const double b1 = pts[i + 1];

      q (f, a1, b1, &area1, &error1, &resabs1, &resasc1);

      q_result = q_result + area1;
      q_abserr = q_abserr + error1;
      q_resabs = q_resabs + resabs1;

      alist[i] = a1;
      blist[i] = b1;
      rlist[i] = area1;
      elist[i] = error1;
      iord[i] = i;
      level[i] = 0;

      if (error1 == resasc1 && error1 != 0.0)
	{
	  ndin[i] = 1;
	}
      else
	{
	  ndin[i] = 0;
	}
    }

  /* Compute the initial error estimate */

  errsum = 0.0;

  for (i = 0; i < nint; i++)
    {
      if (ndin[i])
	{
	  elist[i] = q_abserr;
	}
      errsum = errsum + elist[i];
    }

  /* Sort results into order of decreasing error via the indirection
     array iord[] */

  for (i = 0; i < nint; i++)
    {
      size_t i1 = iord[i];
      size_t e1 = elist[i1];
      size_t i_max = i1;
      size_t j;

      for (j = i + 1; j < nint; j++)
	{
	  size_t i2 = iord[j];
	  size_t e2 = elist[i2];

	  if (e2 >= e1)
	    {
	      i_max = i2;
	      e1 = e2;
	    }
	}

      if (i_max != i1)
	{
	  iord[i] = iord[i_max];
	  iord[i_max] = i1;
	}
    }

  /* Test on accuracy */

  tolerance = GSL_MAX_DBL (epsabs, epsrel * fabs (q_result));

  if (q_abserr <= 100 * GSL_DBL_EPSILON * q_resabs && q_abserr > tolerance)
    {
      *result = q_result;
      *abserr = q_abserr;
      *nqeval = nint;
      *last = 0;

      GSL_ERROR ("cannot reach tolerance because of roundoff error"
		 "on first attempt", GSL_EROUND);
    }
  else if ((q_abserr <= tolerance && q_abserr != q_resabs) || q_abserr == 0)
    {
      *result = q_result;
      *abserr = q_abserr;
      *nqeval = nint;
      *last = 0;
      return GSL_SUCCESS;
    }
  else if (limit == 1)
    {
      *result = q_result;
      *abserr = q_abserr;
      *nqeval = nint;
      *last = 0;
      GSL_ERROR ("a maximum of one iteration was insufficient", GSL_EMAXITER);
    }

  /* Initialization */

  numrl2 = 0;
  rlist2[0] = q_result;

  maxerr_index = iord[0];
  maxerr_value = elist[maxerr_index];

  area = q_result;

  *abserr = GSL_DBL_MAX;

  error_over_large_intervals = errsum;
  ertest = tolerance;

  /* Compare the integral of f(x) with the integral of |f(x)|
     to determine if f(x) covers both positive and negative values */

  if (fabs (q_result) >= (1 - 50 * GSL_DBL_EPSILON) * q_resabs)
    {
      positive_integrand = 1;
    }
  else
    {
      positive_integrand = 0;
    }

  i = nint - 1;	

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
      double resasc1, resasc2;
      double resabs1, resabs2;

      const size_t current_level = level[maxerr_index] + 1;

      i++;

      q (f, a1, b1, &area1, &error1, &resabs1, &resasc1);
      q (f, a2, b2, &area2, &error2, &resabs2, &resasc2);

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

      if (resasc1 != error1 && resasc2 != error2)
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
	  alist[maxerr_index] = a2;	/* blist[maxerr] is already == b2 */
	  rlist[maxerr_index] = area2;
	  elist[maxerr_index] = error2;
	  level[maxerr_index] = current_level;

	  alist[i] = a1;
	  blist[i] = b1;
	  rlist[i] = area1;
	  elist[i] = error1;
	  level[i] = current_level;
	}
      else
	{
	  blist[maxerr_index] = b1;	/* alist[maxerr] is already == a1 */
	  rlist[maxerr_index] = area1;
	  elist[maxerr_index] = error1;
	  level[maxerr_index] = current_level;

	  alist[i] = a2;
	  blist[i] = b2;
	  rlist[i] = area2;
	  elist[i] = error2;
	  level[i] = current_level;
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

      if (disallow_extrapolation)
	{
	  continue;
	}

      error_over_large_intervals += -last_maxerr_value;

      if (current_level < maximum_level)
	{
	  error_over_large_intervals += error12;
	}

      if (!extrapolate)
	{
	  /* test whether the interval to be bisected next is the
	     smallest interval. */
	  if (level[maxerr_index] < maximum_level)
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
	      if (level[maxerr_index] < maximum_level)
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

      if (numrl2 < 2) 
	{
	  goto skip_extrapolation;
	} 

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

    skip_extrapolation:
      maxerr_index = iord[0];
      maxerr_value = elist[maxerr_index];
      nrmax = 0;
      extrapolate = 0;
      maximum_level++;
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

    if (!positive_integrand && max_area < 0.01 * q_resabs)
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

  *nqeval = 2 * (i + 1 - nint) + nint;
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
