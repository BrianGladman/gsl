/* Author:  G. Jungman
   RCS:     $Id$ */

#include <config.h>
#include <gsl_math.h>
#include <gsl_test.h>
#include <gsl_errno.h>
#include <gsl_sum.h>

int
gsl_sum_levin_u_trunc_accel (const double *array,
			     const size_t array_size,
			     double *q_num,
			     double *q_den,
			     double *sum_accel,
			     size_t * n_used,
			     double *sum_plain,
			     double *precision)
{
  return gsl_sum_levin_u_trunc_accel_minmax (array, array_size,
					     0, array_size - 1,
					     q_num, q_den,
					     sum_accel, n_used,
					     sum_plain, precision);
}


int
gsl_sum_levin_u_trunc_accel_minmax (const double *array,
				    const size_t array_size,
				    const size_t min_terms,
				    const size_t max_terms,
				    double *q_num,
				    double *q_den,
				    double *sum_accel,
				    size_t * n_used,
				    double *sum_plain,
				    double *precision)
{
  if (array_size == 0)
    {
      *sum_accel = 0.0;
      *sum_plain = 0.0;
      *n_used    = 0;
      return GSL_SUCCESS;
    }
  else if (array_size == 1)
    {
      *sum_accel = array[0];
      *sum_plain = array[0];
      *n_used    = 1;
      return GSL_SUCCESS;
    }
  else
    {
      const double SMALL = 0.01;
      const size_t nmax = GSL_MAX (max_terms, array_size) - 1;
      double trunc_n = 0.0, trunc_nm1 = 0.0;
      double result_n = 0.0, result_nm1 = 0.0;
      size_t n;
      int better = 0;
      int before = 0;
      int converging = 0;
      double least_trunc = GSL_DBL_MAX;
      double result_least_trunc;

      /* Calculate specified minimum number of terms. No convergence
         tests are made, and no truncation information is stored. */

      for (n = 0; n < min_terms; n++)
	{
	  const double t = array[n];

	  result_nm1 = result_n;
	  gsl_sum_levin_u_trunc_step (t, n, q_num, q_den, 
				      &result_n, sum_plain);
	}

      /* Assume the result after the minimum calculation is the best. */

      result_least_trunc = result_n;

      /* Calculate up to maximum number of terms. Check truncation
         condition. */

      for (; n <= nmax; n++)
	{
	  const double t = array[n];

	  result_nm1 = result_n;
	  gsl_sum_levin_u_trunc_step (t, n, q_num, q_den, 
				      &result_n, sum_plain);

	  trunc_nm1 = trunc_n;
	  trunc_n = fabs (result_n - result_nm1);

	  /* Determine if we are in the convergence region. */

	  better = (trunc_n < trunc_nm1 || trunc_n < SMALL * fabs (result_n));
	  converging = converging || (better && before);
	  before = better;

	  if (converging)
	    {
	      if (trunc_n < least_trunc)
		{
		  /* Found a low truncation point in the convergence
		     region. Save it. */

		  least_trunc = trunc_n;
		  result_least_trunc = result_n;
		}

	      if (fabs (trunc_n / result_n) < 10.0 * GSL_MACH_EPS)
		break;
	    }
	}

      if (converging)
	{
	  /* Stopped in the convergence region. Return result and
	     error estimate. */

	  *sum_accel = result_least_trunc;
	  *precision = fabs (least_trunc / *sum_accel);
	  *n_used = n ;
	  return GSL_SUCCESS;
	}
      else
	{
	  /* Never reached the convergence region. Use the last
	     calculated values. */

	  *sum_accel = result_n;
	  *precision = fabs (trunc_n / result_n);
	  *n_used = n ;
	  return GSL_SUCCESS;
	}
    }
}

int
gsl_sum_levin_u_trunc_step (const double term,
			    const size_t n,
			    double *q_num,
			    double *q_den,
			    double *sum_accel,
			    double *sum_plain)
{
  if (term == 0.0)
    {
      /* This is actually harmless when treated in this way. A term
         which is exactly zero is simply ignored; the state is not
         changed. We return GSL_EZERODIV as an indicator that this
         occured. */

      return GSL_EZERODIV;
    }
  else if (n == 0)
    {
      *sum_accel = term;
      *sum_plain = term;
      q_den[0] = 1.0 / term;
      q_num[0] = 1.0;
      return GSL_SUCCESS;
    }
  else
    {
      double factor = 1.0;
      double ratio = (double) n / (n + 1.0);
      int j;

      *sum_plain += term;
      q_den[n] = 1.0 / (term * (n + 1.0) * (n + 1.0));
      q_num[n] = *sum_plain * q_den[n];

      for (j = n - 1; j >= 0; j--)
	{
	  double c = factor * (j + 1) / (n + 1);
	  factor *= ratio;
	  q_den[j] = q_den[j + 1] - c * q_den[j];
	  q_num[j] = q_num[j + 1] - c * q_num[j];
	}

      *sum_accel = q_num[0] / q_den[0];
      return GSL_SUCCESS;
    }
}
