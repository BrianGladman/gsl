#include <math.h>
#include <float.h>
#include <gsl_errno.h>
#include <gsl_integration.h>

#include "qpsrt.h"
#include "max.h"

int
qagse_impl (double (*f)(double x), double a, double b, 
	    double epsabs, double epsrel, size_t limit,
	    double * result, double * abserr, 
	    double alist[], double blist[], double rlist[], double elist[],
	    size_t iord[], size_t * last, size_t * nqeval,
	    gsl_integration_rule_t * const q)
{
  double q_result, q_abserr, q_defabs, q_resabs ;
  double tolerance,  maxerr_value, area, errsum, erlast ;
  double small, ertest, errbnd ;
  double error_over_large_intervals ;
  double reseps = 0, abseps = 0, correc = 0 ;
  size_t maxerr_index, nrmax, i, nres = 0, numrl2 = 1, ktmin = 0 ;
  int roundoff_type1 = 0, roundoff_type2 = 0, roundoff_type3 ;
  int error_type = 0, error_type2 = 0 ;

  int positive_integrand = 0 ;
  int extrapolate = 0 ;
  int disallow_extrapolation = 0 ;

  double res3la[3], rlist2[52] ;

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

  /* Perform the first integration */

  q (f, a, b, &q_result, &q_abserr, &q_defabs, &q_resabs) ;

  rlist[0] = q_result ;
  elist[0] = q_abserr ;
  iord[0] = 0 ;

  tolerance = max(epsabs, epsrel * fabs(q_result)) ;

  if (q_abserr <= 100 * DBL_EPSILON * q_defabs && q_abserr > tolerance)
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

  rlist2[0] = q_result ;
  maxerr_value = q_abserr ;
  maxerr_index = 0 ;
  area = q_result ;
  errsum = q_abserr ;
  tolerance = DBL_MAX ;

  nrmax = 0 ;

  roundoff_type1 = 0 ;
  roundoff_type2 = 0 ;
  roundoff_type3 = 0 ;

  if(fabs(q_result) >= (1 - 50 * DBL_EPSILON)*q_defabs) 
    {
      positive_integrand = 1 ;
    }
  else
    {
      positive_integrand = 0 ;
    }

  i = 1; 

  do
    {
      /* Bisect the subinterval with the nrmax-th largest error estimate */

      const double left = alist[maxerr_index] ;
      const double right = blist[maxerr_index] ;
      const double midpoint = 0.5 * (left + right) ; 
      
      const double a1 = left, b1 = midpoint ;
      const double a2 = midpoint, b2 = right ;

      double area1 = 0, area2 = 0, area12 = 0 ;
      double error1 = 0, error2 = 0, error12 = 0 ;
      double defab1, defab2 ;
      double resabs = 0 ;

      printf("i=%d, limit=%d\n",i,limit) ;

      q (f, a1, b1, &area1, &error1, &resabs, &defab1) ;
      q (f, a2, b2, &area2, &error2, &resabs, &defab2) ;

      area12 = area1 + area2;
      error12 = error1 + error2;
      erlast = maxerr_value ;

      printf("a1 = %g b1 = %g a2 = %g b2 = %g\n", a1,b1,a2,b2) ;
      printf("area1 = %g area = %g\n", area1, area2) ;
      printf("error1 = %g error2 = %g\n", error1, error2) ;
      printf("error12 = %g maxerr_value = %g\n", error12, maxerr_value) ;
      printf("delta = %g\n", error12 - maxerr_value) ;
      printf("errsum before = %g\n", errsum) ;
      errsum += (error12 - maxerr_value) ;
      printf("errsum after = %g\n", errsum) ;
      area += area12 - rlist[maxerr_index] ;

      if (defab1 != error1 && defab2 != error2) 
	{
	  if (fabs(rlist[maxerr_index] - area12) <= 0.00001 * fabs(area12)
	      && error12 >= 0.99 * maxerr_value)
	    {
	      if (!extrapolate) 
		{
		  roundoff_type1++ ;
		}
	      else
		{
		  roundoff_type2++ ;
		}
	    }
	  if (i >= 10 && error12 > maxerr_value)
	    {
	      roundoff_type3++ ;
	    }
	}
  
      
      if (roundoff_type1 + roundoff_type2 >= 10 || roundoff_type3 >= 20)
	{
	  error_type = 2 ; /* round off error */
	} 

      if (roundoff_type2 >= 5) 
	{
	  error_type2 = 1 ;
	}

      /* set error flag in the case of bad integrand behaviour at
	 a point of the integration range */
      
      {
	double tmp = (1 + 100 * DBL_EPSILON)*(fabs(a2) + 1000 * DBL_MIN) ;
	if (fabs(a1) <= tmp && fabs(b2) <= tmp)
	  {
	    error_type = 4;
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
	 the list of error estimates and select the subinterval with
	 the nrmax-th largest error estimate (to be bisected next) */
      
      qpsrt(limit,i,&maxerr_index,&maxerr_value,elist,iord,&nrmax) ;

      tolerance = max (epsabs, epsrel * fabs(area)) ;

      if (errsum <= tolerance)
	{
	  continue ;
	}

      if(error_type) continue; 

      if(i == 1) /* set up variables on first iteration */
	{
	  small = fabs(b-a)*0.375 ;
	  error_over_large_intervals = errsum ;
	  ertest = errbnd ;
	  rlist2[1] = area ;
	  continue ;
	}
      
      if (disallow_extrapolation) 
	continue ;
            
      error_over_large_intervals += - erlast ;

      if(fabs(b1-a1) > small) 
	{
	  error_over_large_intervals += error12 ;
	}
      
      if(!extrapolate) 
	{
	  /* test whether the interval to be bisected next is the
	     smallest interval. */
	  if(fabs(blist[maxerr_index]-alist[maxerr_index]) > small)
	    continue ;
	  extrapolate = 1 ;
	  nrmax = 1 ;
	}

      if(error_type2 != 1 && error_over_large_intervals > ertest)
	{
	  
	  /* The smallest interval has the largest error.  Before
	     bisecting decrease the sum of the errors over the larger
	     intervals (error_over_large_intervals) and perform
	     extrapolation. */
	  int k ;
	  int id = nrmax ;
	  int jupbnd ;
	  if(i > (2+limit/2)) { 
	    jupbnd = limit+3-i ;
	  } else {
	    jupbnd = i ;
	  }
	  for (k = id ; k < jupbnd ; k++) 
	    {
	      maxerr_index = iord[nrmax];
	      maxerr_value = elist[maxerr_index] ;
	      if(fabs(blist[maxerr_index]-alist[maxerr_index]) > small) 
		continue ;
	      nrmax++ ;
	    }
	}
      
      numrl2++ ; rlist2[numrl2] = area ; 

      qelg(&numrl2, rlist2, &reseps, &abseps, res3la, &nres) ;
      ktmin++ ;
      
      if(ktmin > 5 && *abserr < 0.001*errsum)
	error_type = 5 ;
      
      if(abseps < *abserr) 
	{
	  ktmin = 0 ;
	  *abserr = abseps ;
	  *result = reseps ;
	  correc = error_over_large_intervals ;
	  ertest = max(epsabs,epsrel*fabs(reseps)) ;
	  if(*abserr <= ertest) 
	    goto done ;
	}
      else 
	{
	  /* Prepare bisection of the smallest interval. */
	  
	  if(numrl2 == 0) 
	    disallow_extrapolation = 1 ;

	  if(error_type == 5)
	    goto done ;
	  maxerr_index = iord[0] ;
	  maxerr_value = elist[maxerr_index] ;
	  nrmax = 1 ;
	  extrapolate = 0 ;
	  small *= 0.5 ;
	  error_over_large_intervals = errsum ;
	}

    }  while (i < limit && !error_type && errsum > tolerance) ;
  

done:

  if(*abserr == DBL_MAX) 
    goto 115 ;
  if(error_type == 0 && error_type2 = 0) 
    goto 110 ;
  if(error_type2) 
    *abserr += correc ;
  if(error_type == 0) 
    error_type = 3 ;
  if(result != 0 && area != 0) 
    goto 105 ;
  if(*abserr > errsum) 
    goto 115 ;
  if(area == 0) 
    goto 130 ;
  go to 110 ;

  105 if(*abserr/fabs(result) > errsum/fabs(area)) 
    go to 115 ;

  /*  Test on divergence. */

  110 if(!positive_integrand && max(fabs(result),fabs(area)) < 0.01 * defabs)
    go to 130 ;

  if((result/area) < 0.01 || (result/area) > 100 || errsum > fabs(area))
    error_type = 6 ;
  go to 130 ;

  115 {
    /* Compute global integral sum. */
    
    double result_sum = 0 ;
    size_t k ;
    for (k = 0; k <= i; k++)
      {
	result_sum += rlist[k] ;
      }
    *result = result_sum ;
  }
  



