#include <gsl_math.h>
#include <gsl_test.h>
#include <gsl_errno.h>
#include <gsl_sum.h>

int
gsl_sum_levin_u_with_derivs_accel (const double *array, 
				   const size_t array_size,
				   double *q_num,
				   double *q_den,
				   double *dq_num,
				   double *dq_den,
				   double *dsum,		       
				   double *sum_accel,
				   double *sum_plain,
				   double *precision)
{
  return gsl_sum_levin_u_with_derivs_accel_minmax (array, array_size,
						   0, array_size - 1,
						   q_num, q_den,
						   dq_num, dq_den, dsum,
						   sum_accel, 
						   sum_plain, 
						   precision);
}

int
gsl_sum_levin_u_with_derivs_accel_minmax (const double *array, 
					  const size_t array_size,
					  const size_t min_terms,
					  const size_t max_terms,
					  double *q_num,
					  double *q_den,
					  double *dq_num,
					  double *dq_den,
					  double *dsum,
					  double *sum_accel,
					  double *sum_plain,
					  double *precision)
{
  if (array_size == 0)
    {
      *sum_accel = 0.0;
      *sum_plain = 0.0;
      return GSL_SUCCESS;
    }
  else if (array_size == 1)
    {
      *sum_accel = array[0];
      *sum_plain = array[0];
      return GSL_SUCCESS;
    }
  else
    {
      const double SMALL = 0.01;
      const size_t nmax = MAX (max_terms, array_size) - 1;
      double noise_n = 0.0, noise_nm1 = 0.0;
      double trunc_n = 0.0, trunc_nm1 = 0.0;
      double result_n = 0.0, result_nm1 = 0.0;
      double variance = 0;
      size_t n;
      unsigned int i;
      int better = 0;
      int before = 0;
      int converging = 0;
      double least_trunc = DBL_MAX;
      double least_trunc_noise = DBL_MAX;
      double least_trunc_result;

      /* Calculate specified minimum number of terms.  No convergence
	 tests are made, and no truncation information is stored.  */

      for (n = 0; n < min_terms; n++)
	{
	  const double t = array[n];
	  result_nm1 = result_n;
	  gsl_sum_levin_u_with_derivs_step (t, n, nmax, q_num, q_den, 
					    dq_num, dq_den, dsum,
					    &result_n,sum_plain);
	}

      least_trunc_result = result_n ;

      variance = 0 ;
      for (i = 0 ; i <= n; i++) 
	{
	  double dn = dsum[i] * GSL_MACH_EPS * array[i];
	  variance += dn * dn ;
	}
      noise_n = sqrt(variance) ;

      /* Calculate up to maximum number of terms.  Check truncation
         condition.  */

      for (; n <= nmax; n++)
	{
	  const double t = array[n];

	  result_nm1 = result_n;
	  gsl_sum_levin_u_with_derivs_step (t, n, nmax, q_num, q_den, 
					    dq_num, dq_den, dsum,
					    &result_n, sum_plain);

	  trunc_nm1 = trunc_n;
	  trunc_n = fabs (result_n - result_nm1);

	  noise_nm1 = noise_n;
	  variance = 0 ;
	  
	  for (i = 0 ; i <= n; i++) 
	    {
	      double dn = dsum[i] * GSL_MACH_EPS * array[i];
	      variance += dn * dn ;
	    }
	  
	  noise_n = sqrt(variance) ;

	  /* Determine if we are in the convergence region.  */

	  better = (trunc_n < trunc_nm1 || trunc_n < SMALL * fabs (result_n));
	  converging = converging || (better && before);
	  before = better;

	  printf("trun = %g  noise = %g\n", trunc_n, noise_n) ;

	  if (converging)
	    {

	      if (trunc_n < least_trunc)
		{
		  /* Found a low truncation point in
		   * the convergence region. Save it.
		   */
		  least_trunc_result = result_n;
		  least_trunc = trunc_n;
		  least_trunc_noise = noise_n;
		}

	      if (trunc_n < 3.0 * noise_n)
		{
		  printf("trun = %g < 3* noise = %g\n", trunc_n, noise_n) ;
		  break;
		}

	      if (trunc_n < 10.0 * GSL_MACH_EPS * fabs(result_n))
		break;
	    }

	}

      if (converging)
	{
	  /* Stopped in the convergence region.  Return result and
	     error estimate.  */

	  *sum_accel = least_trunc_result;
	  *precision = MAX(least_trunc,least_trunc_noise) / fabs(*sum_accel);
	  return GSL_SUCCESS;
	}
      else
	{
	  /* Never reached the convergence region.  Use the last
	     calculated values.  */

	  *sum_accel = result_n;
	  *precision = MAX(trunc_n, noise_n) / fabs(result_n);
	  return GSL_SUCCESS;
	}
    }
}


int
gsl_sum_levin_u_with_derivs_step (const double term,
				  const size_t n,
				  const size_t nmax,
				  double *q_num,
				  double *q_den,
				  double *dq_num,
				  double *dq_den,
				  double *dsum,
				  double *sum_accel,
				  double *sum_plain)
{

#define IN(i,j) ((i)*nmax + j)


  if (n == 0)
    {
      *sum_accel = term;
      *sum_plain = term;

      q_den[0] = 1.0 / term;
      q_num[0] = 1.0;

      dq_den[IN(0,0)] = -1.0 / (term * term);
      dq_num[IN(0,0)] = 0.0;

      dsum[0] = 1.0;

      return GSL_SUCCESS;
    }
  else
    {
      double factor = 1.0;
      double ratio = (double) n / (n + 1.0);
      unsigned int i;
      int j;

      *sum_plain += term;

      q_den[n] = 1.0 / (term * (n + 1.0) * (n + 1.0));
      q_num[n] = *sum_plain * q_den[n];

      for (i = 0; i < n; i++)
	{
	  dq_den[IN(i,n)] = 0;
	  dq_num[IN(i,n)] = q_den[i] ;
	}

      dq_den[IN(n,n)] = -q_den[n]/ term;
      dq_num[IN(n,n)] = q_den[n] + (*sum_plain) * (dq_den[IN(n,n)]);

      for (j = n - 1; j >= 0; j--)
	{
	  double c = factor * (j + 1) / (n + 1);
	  factor *= ratio;
	  q_den[j] = q_den[j + 1] - c * q_den[j];
	  q_num[j] = q_num[j + 1] - c * q_num[j];

	  for (i = 0; i < n; i++)
	    {
	      dq_den[IN(i,j)] = dq_den[IN(i,j+1)] - c * dq_den[IN(i,j)]; 
	      dq_num[IN(i,j)] = dq_num[IN(i,j+1)] - c * dq_num[IN(i,j)]; 
	    } 
	  
	  dq_den[IN(n,j)] = dq_den[IN(n,j+1)] ;
	  dq_num[IN(n,j)] = dq_num[IN(n,j+1)] ;
	}

      *sum_accel = q_num[0] / q_den[0];

      for (i = 0; i <= n; i++)
	{
	  dsum[i] = ((dq_num[IN(i,0)] - q_num[0] * dq_den[IN(i,0)] / q_den[0]) 
		     / q_den[0]);
	}

      return GSL_SUCCESS;
    }
}

