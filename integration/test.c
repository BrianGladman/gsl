#include <stdio.h>
#include <math.h>
#include <gsl_integration.h>
#include <gsl_test.h>

#include "tests.h"

double alpha = 0 ;

int main (void)
{
  size_t neval ;

  /* Test a smooth function. This will find any discrepancies in the
     Gauss-Kronrod calculation */

  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result = 7.716049357767090777E-02;
    double exp_abserr = 2.990224871000550874E-06;
    double exp_resabs = 7.716049357767090777E-02;
    double exp_resasc = 4.434273814139995384E-02;

    alpha = 2.6 ;
    gsl_integration_qk15 (book1, 0.0, 1.0, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk15(book1) smooth result") ;
    gsl_test_rel(abserr,exp_abserr,1e-15,"qk15(book1) smooth abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk15(book1) smooth resabs") ;    
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk15(book1) smooth resasc") ;
  }


  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result = 7.716049379303084599E-02;
    double exp_abserr = 9.424302194248481445E-08;
    double exp_resabs = 7.716049379303084599E-02;
    double exp_resasc = 4.434311425038358484E-02;

    alpha = 2.6 ;
    gsl_integration_qk21 (book1, 0.0, 1.0, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-11,"qk21(book1) smooth result") ;
    gsl_test_rel(abserr,exp_abserr,1e-15,"qk21(book1) smooth abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk21(book1) smooth resabs") ;    
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk21(book1) smooth resasc") ;
  }


  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result = 7.716049382494900855E-02;
    double exp_abserr = 1.713503193600029893E-09;
    double exp_resabs = 7.716049382494900855E-02;
    double exp_resasc = 4.427995051868838933E-02;

    alpha = 2.6 ;
    gsl_integration_qk31 (book1, 0.0, 1.0, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk31(book1) smooth result") ;
    gsl_test_rel(abserr,exp_abserr,1e-15,"qk31(book1) smooth abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk31(book1) smooth resabs") ;    
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk31(book1) smooth resasc") ;
  }


  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result = 7.716049382681375302E-02;
    double exp_abserr = 9.576386660975511224E-11;
    double exp_resabs = 7.716049382681375302E-02;
    double exp_resasc = 4.421521169637691873E-02;
    
    alpha = 2.6 ;
    gsl_integration_qk41 (book1, 0.0, 1.0, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk41(book1) smooth result") ;
    gsl_test_rel(abserr,exp_abserr,1e-7,"qk41(book1) smooth abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk41(book1) smooth resabs") ;    
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk41(book1) smooth resasc") ;
  }

  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result = 7.716049382708510540E-02;
    double exp_abserr = 1.002079980317363772E-11;
    double exp_resabs = 7.716049382708510540E-02;
    double exp_resasc = 4.416474291216854892E-02;

    alpha = 2.6 ;
    gsl_integration_qk51 (book1, 0.0, 1.0, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk51(book1) smooth result") ;
    gsl_test_rel(abserr,exp_abserr,1e-15,"qk51(book1) smooth abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk51(book1) smooth resabs") ;    
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk51(book1) smooth resasc") ;
  }


  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result = 7.716049382713800753E-02;
    double exp_abserr = 1.566060362296155616E-12;
    double exp_resabs = 7.716049382713800753E-02;
    double exp_resasc = 4.419287685934316506E-02;
    
    alpha = 2.6 ;
    gsl_integration_qk61 (book1, 0.0, 1.0, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk61(book1) smooth result") ;
    gsl_test_rel(abserr,exp_abserr,1e-5,"qk61(book1) smooth abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk61(book1) smooth resabs") ;    
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk61(book1) smooth resasc") ;
  }

  /* Now test a function with a singularity. This will find any discrepancies
     in the abserr calculation */

  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result = 1.555688196612745777E+01;
    double exp_abserr = 2.350164577239293706E+01;
    double exp_resabs = 1.555688196612745777E+01;
    double exp_resasc = 2.350164577239293706E+01;

    alpha = -0.9 ;
    gsl_integration_qk15 (book1, 0.0, 1.0, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk15(book1) singular result") ;
    gsl_test_rel(abserr,exp_abserr,1e-15,"qk15(book1) singular abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk15(book1) singular resabs") ;    
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk15(book1) singular resasc") ;
  }


  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result = 1.799045317938126232E+01;
    double exp_abserr = 2.782360287710622515E+01;
    double exp_resabs = 1.799045317938126232E+01;
    double exp_resasc = 2.782360287710622515E+01;

    alpha = -0.9 ;
    gsl_integration_qk21 (book1, 0.0, 1.0, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk21(book1) singular result") ;
    gsl_test_rel(abserr,exp_abserr,1e-15,"qk21(book1) singular abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk21(book1) singular resabs") ;    
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk21(book1) singular resasc") ;
  }


  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result = 2.081873305159121657E+01;
    double exp_abserr = 3.296500137482590276E+01;
    double exp_resabs = 2.081873305159121301E+01;
    double exp_resasc = 3.296500137482590276E+01;

    alpha = -0.9 ;
    gsl_integration_qk31 (book1, 0.0, 1.0, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk31(book1) singular result") ;
    gsl_test_rel(abserr,exp_abserr,1e-15,"qk31(book1) singular abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk31(book1) singular resabs") ;    
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk31(book1) singular resasc") ;
  }


  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result = 2.288677623903126701E+01;
    double exp_abserr = 3.671538820274916048E+01;
    double exp_resabs = 2.288677623903126701E+01;
    double exp_resasc = 3.671538820274916048E+01;
    
    alpha = -0.9 ;
    gsl_integration_qk41 (book1, 0.0, 1.0, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk41(book1) singular result") ;
    gsl_test_rel(abserr,exp_abserr,1e-15,"qk41(book1) singular abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk41(book1) singular resabs") ;    
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk41(book1) singular resasc") ;
  }

  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result = 2.449953612016972215E+01;
    double exp_abserr = 3.967771249391228849E+01;
    double exp_resabs = 2.449953612016972215E+01;
    double exp_resasc = 3.967771249391228849E+01;

    alpha = -0.9 ;
    gsl_integration_qk51 (book1, 0.0, 1.0, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk51(book1) singular result") ;
    gsl_test_rel(abserr,exp_abserr,1e-15,"qk51(book1) singular abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk51(book1) singular resabs") ;    
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk51(book1) singular resasc") ;
  }


  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result = 2.583030240976628988E+01;
    double exp_abserr = 4.213750493076978643E+01;
    double exp_resabs = 2.583030240976628988E+01;
    double exp_resasc = 4.213750493076978643E+01;

    alpha = -0.9 ;
    gsl_integration_qk61 (book1, 0.0, 1.0, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk61(book1) singular result") ;
    gsl_test_rel(abserr,exp_abserr,1e-15,"qk61(book1) singular abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk61(book1) singular resabs") ;    
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk61(book1) singular resasc") ;
  }



#ifdef JUNK
 gsl_integration_qng(f,0.0,1.0,0.0,1e-5,&result,&abserr,&neval) ;
  
  printf("qng: result = %.18g, abserr = %.18g, neval = %d\n",
	 result, abserr, neval) ;


  {
    double alist[1000], blist[1000], rlist[1000], elist[1000];
    size_t iord[1000] ;
    size_t last;
    result = 0 ; abserr=0; neval=0  ;
    gsl_integration_qage(f, 0.0, 1.0, 0, 1e-10, 6, 1000,
			 alist, blist, rlist, elist, iord, &last,
			 &result, &abserr, &neval) ;
    printf("qage: result = %.18g, abserr = %.18g, neval = %d\n",
	   result, abserr, neval) ;
  }

  {
    double alist[1000], blist[1000], rlist[1000], elist[1000];
    size_t iord[1000] ;
    size_t last;
    result = 0 ; abserr=0; neval=0  ;
    status = gsl_integration_qagse(f, 0.0, 1.0, 1e-20, 1e-17, 10,
			 alist, blist, rlist, elist, iord, &last,
			 &result, &abserr, &neval) ;
    printf("qagse: result = %.18g, abserr = %.18g, neval = %d\n",
	   result, abserr, neval) ;
    printf("status=%d\n",status) ;
  }

  {
    double alist[1000], blist[1000], rlist[1000], elist[1000];
    size_t iord[1000] ;
    size_t last;
    result = 0 ; abserr=0; neval=0  ;
    status = gsl_integration_qagse(g, 0.0, 1.0, 0, 1e-10, 10,
			 alist, blist, rlist, elist, iord, &last,
			 &result, &abserr, &neval) ;
    printf("qagse: result = %.18g, abserr = %.18g, neval = %d\n",
	   result, abserr, neval) ;
    printf("status=%d\n",status) ;
  }


  {
    double alist[1000], blist[1000], rlist[1000], elist[1000];
    size_t iord[1000] ;
    size_t last;
    result = 0 ; abserr=0; neval=0  ;
    status = gsl_integration_qagse(stepfn, 0.0, 1.0, 0, 1e-10, 3,
				   alist, blist, rlist, elist, iord, &last,
				   &result, &abserr, &neval) ;
    printf("qagse: result = %.18g, abserr = %.18g, neval = %d\n",
	   result, abserr, neval) ;
    printf("status=%d\n",status) ;
  }
#endif

  return gsl_test_summary() ;
} 

