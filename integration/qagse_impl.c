#include <math.h>
#include <float.h>
#include <gsl_errno.h>
#include <gsl_integration.h>

#include "qpsrt.h"
#include "max.h"

int
gsl_integration_qagse_impl (double (*f)(double x), 
			    const double a, const double b, 
			    const double epsabs, const double epsrel, 
			    const size_t limit,
			    double * result, double * abserr, 
			    double alist[], double blist[], 
			    double rlist[], double elist[],
			    size_t iord[], size_t * last, size_t * nqeval,
			    gsl_integration_rule_t * const q)
{
  double q_result, q_abserr, q_defabs, q_resabs ;
  double tolerance,  maxerr_value, area, errsum, last_maxerr_value ;
  double small = 0, ertest = 0 ;
  double error_over_large_intervals = 0 ;
  double reseps = 0, abseps = 0, correc = 0 ;
  size_t maxerr_index, nrmax = 0, i = 0, nres = 0, numrl2 = 1, ktmin = 0 ;
  int roundoff_type1 = 0, roundoff_type2 = 0, roundoff_type3 = 0 ;
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
  
  /* Test on accuracy */

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

  printf("a = %g b = %g result = %g err = %g\n", a,b,q_result,q_abserr) ;

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

  /* Initialization */

  rlist2[0] = q_result ;
  maxerr_value = q_abserr ;
  maxerr_index = 0 ;
  area = q_result ;
  errsum = q_abserr ;
  *abserr = DBL_MAX ;

  /* Compare the integral of f(x) with the integral of |f(x)|
     to determine if f(x) covers both positive and negative values */

  if(fabs(q_result) >= (1 - 50 * DBL_EPSILON)*q_defabs) 
    {
      positive_integrand = 1 ;
    }
  else
    {
      positive_integrand = 0 ;
    }

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

      i++ ;
      printf("i=%d, limit=%d\n",i,limit) ;

      q (f, a1, b1, &area1, &error1, &resabs, &defab1) ;
      q (f, a2, b2, &area2, &error2, &resabs, &defab2) ;

      area12 = area1 + area2;
      error12 = error1 + error2;
      last_maxerr_value = maxerr_value ;

      /* Improve previous approximations to the integral and test for
         accuracy */

      printf("a1 = %g b1 = %g a2 = %g b2 = %g\n", a1,b1,a2,b2) ;
      printf("area1 = %g area = %g\n", area1, area2) ;
      printf("error1 = %g error2 = %g\n", error1, error2) ;
      printf("error12 = %g maxerr_value = %g\n", error12, maxerr_value) ;
      printf("delta = %g\n", error12 - maxerr_value) ;
      printf("errsum before = %g\n", errsum) ;
      errsum += (error12 - maxerr_value) ;
      printf("errsum after = %g\n", errsum) ;
      area += area12 - rlist[maxerr_index] ;

      tolerance = max (epsabs, epsrel * fabs(area)) ;
      printf("result = %.15f abserr = %.16g\n", area, errsum) ;
      printf("errsum = %g  vs tolerance = %g\n", errsum, tolerance) ;

      printf("defab1 = %g, error1 = %g\n", defab1, error1) ;
      printf("defab2 = %g, error2 = %g\n", defab2, error2) ;

      if (defab1 == error1 || defab2 == error2) 
	printf(" JUMP\n") ;

      if (defab1 != error1 && defab2 != error2) 
	{
	  printf("defab1 == error1 || defab2 == error2\n") ;
	  printf("fabs(rlist[maxerr_index] - area12) = %g\n",
		 fabs(rlist[maxerr_index] - area12)) ;
	  printf("1e-5 * fabs(area12 = %g\n", 1e-5 * fabs(area12) ) ;
	  printf("error12 = %g\n", error12) ;
	  printf("0.99 * maxerr_value = %g\n", 0.99 * maxerr_value) ;
	  if (fabs(rlist[maxerr_index] - area12) <= 1e-5 * fabs(area12)
	      && error12 >= 0.99 * maxerr_value)
	    {
	      if (!extrapolate) 
		{
		  roundoff_type1++ ;
		  printf("roundoff_type1 increased to %d\n", roundoff_type1);
		}
	      else
		{
		  roundoff_type2++ ;
		  printf("roundoff_type2 increased to %d\n", roundoff_type2);
		}
	    }
	  if (i > 10 && error12 > maxerr_value)
	    {
	      roundoff_type3++ ;
	      printf("roundoff_type3 increased to %d\n", roundoff_type3);
	    }
	}
  
      /* Test for roundoff and eventually set error flag */

      
      if (roundoff_type1 + roundoff_type2 >= 10 || roundoff_type3 >= 20)
	{
	  error_type = 2 ; /* round off error */
	} 

      if (roundoff_type2 >= 5) 
	{
	  printf("roundoff_type2 >=5, setting error_type2\n") ;
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

      
      if (errsum <= tolerance)
	{
	  printf("BREAK\n") ;
	  break ;
	}

      if(error_type) 
	break; 

      if (i == limit)
	break ;

      if(i == 1) /* set up variables on first iteration */
	{
	  small = fabs(b-a)*0.375 ;
	  error_over_large_intervals = errsum ;
	  ertest = tolerance ;
	  rlist2[1] = area ;
	  continue ;
	}
      
      if (disallow_extrapolation)
	{
	  printf("extrapolation is OFF, continuing...\n") ;
	  continue ;
	} ;

      error_over_large_intervals += - last_maxerr_value ;

      if(fabs(b1-a1) > small) 
	{
	  error_over_large_intervals += error12 ;
	}
      
      if(!extrapolate) 
	{
	  /* test whether the interval to be bisected next is the
	     smallest interval. */
	  if(fabs(blist[maxerr_index]-alist[maxerr_index]) > small)
	    {
	      printf("NOT extrapolating, continuing\n") ;
	      continue ;
	    }
	  extrapolate = 1 ;
	  nrmax = 1 ;
	}

      printf("Checking smallest interval %g vs %g\n",
	     error_over_large_intervals, ertest) ;

      if(!error_type2 && error_over_large_intervals > ertest)
	{
	  /* The smallest interval has the largest error.  Before
	     bisecting decrease the sum of the errors over the larger
	     intervals (error_over_large_intervals) and perform
	     extrapolation. */
	  int k, flag = 0 ;
	  int id = nrmax ;
	  int jupbnd ;
	  if(i > (2+limit/2)) {   /* FIXME */
	    jupbnd = limit+3-i ;
	  } else {
	    jupbnd = i ;
	  }
	  printf("Loop will run from %d to %d - 1\n",id, jupbnd) ;
	  for (k = id ; k < jupbnd && !flag; k++) 
	    {
	      printf("nrmax=%d\n",nrmax) ;
	      maxerr_index = iord[nrmax];
	      maxerr_value = elist[maxerr_index] ;
	      if(fabs(blist[maxerr_index]-alist[maxerr_index]) > small) 
		{
		  flag = 1 ; 
		  break ;
		}
	      nrmax++ ;
	    }
	  if (flag) continue ;
	}
      
      /* Perform extrapolation */

      numrl2++ ; 
      rlist2[numrl2] = area ; 
      
      printf("Before qelg numrl2 = %d\n",numrl2) ;
      printf("rlist2[0,1,2] = %g %g %g\n",rlist2[0],rlist2[1],rlist2[2]) ;
      gsl_integration_qelg(&numrl2, rlist2, &reseps, &abseps, res3la, &nres) ;
      printf("After qelg numrl2 = %d\n",numrl2) ;

      ktmin++ ;

      printf("ktmin = %d  *abserr = %g, 0.001*errsum = %g\n",
	     ktmin,*abserr,0.001*errsum) ;
      
      if(ktmin > 5 && *abserr < 0.001 * errsum)
	error_type = 5 ;
      
      if(abseps < *abserr) 
	{
	  ktmin = 0 ;
	  printf("setting abserr <- abseps = %g\n",abseps) ;
	  *abserr = abseps ;
	  *result = reseps ;
	  correc = error_over_large_intervals ;
	  ertest = max(epsabs,epsrel*fabs(reseps)) ;
	  if(*abserr <= ertest) 
	    break ;
	}

      /* Prepare bisection of the smallest interval. */
      
      if(numrl2 == 0) 
	{
	  printf("Turning off extrapolation because numrl2 == 0\n") ;
	  disallow_extrapolation = 1 ;
	}

      if(error_type == 5)
	break ;

      maxerr_index = iord[0] ;
      maxerr_value = elist[maxerr_index] ;
      nrmax = 0 ;
      extrapolate = 0 ;
      small *= 0.5 ;
      error_over_large_intervals = errsum ;

    }  while (i < limit && !error_type && errsum > tolerance) ;
  
  printf("error_type = %d, error_type2 = %d\n", error_type, error_type2) ;

  if(*abserr == DBL_MAX) 
    goto compute_result ;
  if(error_type == 0 && error_type2 == 0) 
    goto test_divergence ;
  if(error_type2) {
    printf("adding correction of %g\n", correc) ;
    *abserr += correc ;
  }
  if(error_type == 0) 
    error_type = 3 ;
  if(result != 0 && area != 0) 
    goto check_error ;
  if(*abserr > errsum) 
    goto compute_result ;
  if(area == 0) 
    goto fixmefixme ;

  goto test_divergence ;

check_error:
  if(*abserr/fabs(*result) > errsum/fabs(area)) 
    goto compute_result ;

  /*  Test on divergence. */

test_divergence:
  printf("testing divergence\n") ;
  if(!positive_integrand && max(fabs(*result),fabs(area)) < 0.01 * q_defabs)
    goto fixmefixme ;

  if((*result/area) < 0.01 || (*result/area) > 100 || errsum > fabs(area))
    error_type = 6 ;
  goto fixmefixme ;

compute_result:
  printf("computing result\n") ;
  {
    /* Compute global integral sum. */
    
    double result_sum = 0 ;
    size_t k ;
    for (k = 0; k <= i; k++)
      {
	result_sum += rlist[k] ;
      }
    *result = result_sum ;
  }

  *abserr = errsum ;

fixmefixme: 

  *nqeval = 2*i + 1;

  return 0 ;

}


