#include <math.h>
#include <float.h>
#include <gsl_errno.h>
#include <gsl_integration.h>

int
gsl_integration_qage (double (*f)(double x),
		      double a, double b,
		      double epsabs, double epsrel,
		      int key,
		      size_t limit,
		      double alist[], double blist[], double rlist[], 
		      double elist[], size_t iord[], size_t * last,
		      double * result, double * abserr, double * neval)
{
  gsl_integration_rule_t * integration_rule ;

  if (key < GSL_INTEG_GAUSS15)
    {
      key = GSL_INTEG_GAUSS15 ;
    } 
  else if (key > GSL_INTEG_GAUSS61) 
    {
      key = GSL_INTEG_GAUSS61 ;
    }

  switch (key) 
    {
    case GSL_INTEG_GAUSS15:
      integration_rule = &gsl_integration_qk15 ;
      break ;
    case GSL_INTEG_GAUSS21:
      integration_rule = &gsl_integration_qk21 ;
      break ;
    case GSL_INTEG_GAUSS31:
      integration_rule = &gsl_integration_qk31 ; 
      break ;
    case GSL_INTEG_GAUSS41:
      integration_rule = &gsl_integration_qk41 ;
      break ;      
    case GSL_INTEG_GAUSS51:
      integration_rule = &gsl_integration_qk51 ;
      break ;      
    case GSL_INTEG_GAUSS61:
      integration_rule = &gsl_integration_qk61 ;
      break ;      
    }

  gsl_integration_qage_impl (f, a, b, epsabs, epsrel, limit,
			     alist, blist, rlist, elist, iord, last,
			     result, abserr, neval, integration_rule) ;
}

int
gsl_integration_qage_impl (double (*f)(double x),
			   const double a, const double b,
			   const double epsabs, const double epsrel,
			   const size_t limit,
			   double alist[], double blist[], double rlist[], 
			   double elist[], size_t iord[], size_t * last,
			   double * result, double * abserr, double * nqeval,
			   const gsl_integration_rule_t * const q)
{
  double q_result, q_abserr, q_defabs, q_resabs ;
  double tolerance,  maxerr_value, area, errsum ;
  size_t maxerr_index, nrmax, i ;
  int roundoff_type1 = 0, roundoff_type2 = 0, error_type = 0 ;

  alist[0] = a ;
  blist[0] = b ;
  rlist[0] = 0 ;
  elist[0] = 0 ;
  iord[0] = 0 ;
  
  if (epsabs <= 0 && (epsrel < 50 * DBL_EPSILON || epsrel < 0.5e-28))
    {
      * result = 0;
      * abserr = 0;
      * nqeval = 0;
      GSL_ERROR ("tolerance cannot be acheived with given epsabs and epsrel",
		 GSL_EBADTOL);
    };

  /* perform the first integration */

  q (f, a, b, &q_result, &q_abserr, &q_defabs, &q_resabs) ;

  rlist[0] = q_result ;
  elist[0] = q_abserr ;
  iord[0] = 0 ;

  /* Test on accuracy */
  {
    const double absolute = epsabs ;
    const double relative = epsrel * fabs(q_result) ;
    
    if (absolute > relative) 
      {
	tolerance = absolute ;
      }
    else 
      {
	tolerance = relative ;
      }
  }

  if (q_abserr <= 50 * DBL_EPSILON * q_defabs && q_abserr > tolerance)
    {
      *result = q_result ;
      *abserr = q_abserr ;
      *nqeval = 1 ;
      *last = 0 ;
	
      GSL_ERROR("cannot reach tolerance because of roundoff error", 
		GSL_EROUND) ;
    } 
  else if ((q_abserr <= tolerance && q_abserr != q_resabs) || q_abserr == 0)
    {
      *result = q_result ;
      *abserr = q_abserr ;
      *nqeval = 1 ;
      *last = 0 ;
      return GSL_SUCCESS ;
    }
  else if (limit == 1) 
    {
      *result = q_result ;
      *abserr = q_abserr ;
      *nqeval = 1 ;
      *last = 0 ;
      GSL_ERROR("maximum of one iteration was insufficient", GSL_EMAXITER) ;
    }

  maxerr_value = q_abserr ;
  maxerr_index = 0 ;
  area = q_result ;
  errsum = q_abserr ;
  nrmax = 1 ;

  for (i = 1; i < limit; i++)
    {
      /* Bisect the subinterval with the largest error estimate */

      const double left = alist[maxerr_index] ;
      const double right = blist[maxerr_index] ;
      const double midpoint = 0.5 * (left + right) ; 

      const double a1 = left ; 
      const double b1 = midpoint ;
      
      const double a2 = midpoint ; 
      const double b2 = right ;

      double area1, area2, area12 ;
      double error1, error2, error12 ;
      double defab1, defab2 ;
      double resabs ;

      q (f, a1, b1, &area1, &error1, &resabs, &defab1) ;
      q (f, a2, b2, &area2, &error2, &resabs, &defab2) ;
      
      area12 = area1 + area2;
      error12 = error1 + error2;
      errsum += (error12 - maxerr_value) ;
      area += area12 - rlist[maxerr_index] ;
      
      if (defab1 != error1 && defab2 != error2) 
	{
	  if (fabs(rlist[maxerr_index] - area12) < 0.0001 * abs(area12)
	      && error12 > 0.99 * maxerr_value)
	    {
	      roundoff_type1++ ;
	    }
	  if (i > 10 && error12 > maxerr_value)
	    {
	      roundoff_type2++ ;
	    }
	}
      
      {
	const double absolute = epsabs ;
	const double relative = epsrel * fabs(area) ;
      
	if (absolute > relative) 
	  {
	    tolerance = absolute ;
	  }
	else 
	  {
	    tolerance = relative ;
	  }
      }

      if (errsum > tolerance)
	{
	  if (roundoff_type1 >= 6 || roundoff_type2 >=20)
	    {
	      error_type = 2 ;
	    } 

	  /* set error flag in the case that the number of
             subintervals equals limit */
	  
	  if (i == limit)
	    {
	      error_type = 1;
	    }

	  /* set error flag in the case of bad integrand behaviour at
	     a point of the integration range */

	  {
	    double tmp = 1 + 100 * DBL_EPSILON *(fabs(a2) + 1000 * DBL_MIN) ;
	    if (fabs(a1) <= tmp || fabs(a2) <= tmp)
	      {
		error_type = 3;
	      }
	  }
	}
      /* append the newly-created intervals to the list */

      if (error2 > error1) 
	{
	  alist[maxerr_index] = a2 ; /* already done blist[maxerr] = b2 */
	  rlist[maxerr_index] = area2 ;
	  elist[maxerr_index] = error2 ;

	  alist[i] = a1 ;
	  blist[i] = b1 ;
	  rlist[i] = area1 ;
	  elist[i] = error1 ;
	}
      else
	{
	  blist[maxerr_index] = b1 ;  /* alist[maxerr] is already == a1 */
	  rlist[maxerr_index] = area1 ;
	  elist[maxerr_index] = error1;

	  alist[i] = a2 ;
	  blist[i] = b2 ;
	  rlist[i] = area2 ;
	  elist[i] = error2;
	}

      /* call subroutine dqpsrt to maintain the descending ordering in
	 the list of error estimates and select the subinterval with the
	 largest error estimate (to be bisected next) */
      
      dqpsrt(limit,i,&maxerr_index,&maxerr_value,elist,iord,&nrmax) ;
      
      if (error_type != 0 || errsum <= tolerance)
	last ;
    }
  
  {
    double result_sum = 0 ;
    size_t k ;
    for (k = 0; k <= i; k++)
      {
	result_sum += rlist[k] ;
      }
    * result = result_sum ;
  }
  
  *abserr = errsum ;
  
}
  



