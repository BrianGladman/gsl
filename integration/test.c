#include <stdio.h>
#include <math.h>
#include <gsl_integration.h>
#include <gsl_errno.h>
#include <gsl_test.h>

#include "tests.h"

double alpha = 0 ;

void my_error_handler (const char *reason, const char *file,
		       int line, int err);

int main (void)
{

  gsl_set_error_handler (&my_error_handler);

  /* Test the basic rules with a smooth positive function. This will
     find some types of discrepancies in the Gauss-Kronrod calculation. */

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
    gsl_test_rel(result,exp_result,1e-15,"qk21(book1) smooth result") ;
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
    /* FIXME: seems to be a loss of precision on the error here */
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
    /* FIXME: seems to be a loss of precision on the error here */
    gsl_test_rel(abserr,exp_abserr,1e-5,"qk61(book1) smooth abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk61(book1) smooth resabs") ;    
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk61(book1) smooth resasc") ;
  }

  /* Now test the basic rules with a positive function that has a
     singularity. This will find some types of discrepancies in the
     abserr calculation. */

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

  /* Test the basic rules with a smooth oscillating function, over a
     "random" range. This will find different types of discrepancies
     in the Gauss-Kronrod calculation. */

  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result =-7.238969575483799046E-01;
    double exp_abserr = 8.760080200939757174E-06;
    double exp_resabs = 1.165564172429140788E+00;
    double exp_resasc = 9.334560307787327371E-01;
    
    alpha = 1.3 ;
    gsl_integration_qk15 (book3, 0.3, 2.71, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk15(book3) oscill result") ;
    gsl_test_rel(abserr,exp_abserr,1e-15,"qk15(book3) oscill abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk15(book3) oscill resabs") ;
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk15(book3) oscill resasc") ;
  }

  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result =-7.238969575482959717E-01;
    double exp_abserr = 7.999213141433641888E-11;
    double exp_resabs = 1.150829032708484023E+00;
    double exp_resasc = 9.297591249133687619E-01;
    
    alpha = 1.3 ;
    gsl_integration_qk21 (book3, 0.3, 2.71, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk21(book3) oscill result") ;
    gsl_test_rel(abserr,exp_abserr,1e-15,"qk21(book3) oscill abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk21(book3) oscill resabs") ;
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk21(book3) oscill resasc") ;
  }

  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result =-7.238969575482959717E-01;
    double exp_abserr = 1.285805464427459261E-14;
    double exp_resabs = 1.158150602093290571E+00;
    double exp_resasc = 9.277828092501518853E-01;
    
    alpha = 1.3 ;
    gsl_integration_qk31 (book3, 0.3, 2.71, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk31(book3) oscill result") ;
    gsl_test_rel(abserr,exp_abserr,1e-15,"qk31(book3) oscill abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk31(book3) oscill resabs") ;
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk31(book3) oscill resasc") ;
  }

  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result =-7.238969575482959717E-01;
    double exp_abserr = 1.286535726271015626E-14;
    double exp_resabs = 1.158808363486595328E+00;
    double exp_resasc = 9.264382258645686985E-01;
    
    alpha = 1.3 ;
    gsl_integration_qk41 (book3, 0.3, 2.71, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk41(book3) oscill result") ;
    gsl_test_rel(abserr,exp_abserr,1e-15,"qk41(book3) oscill abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk41(book3) oscill resabs") ;
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk41(book3) oscill resasc") ;
  }

  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result =-7.238969575482961938E-01;
    double exp_abserr = 1.285290995039385778E-14;
    double exp_resabs = 1.157687209264406381E+00;
    double exp_resasc = 9.264666884071264263E-01;
    
    alpha = 1.3 ;
    gsl_integration_qk51 (book3, 0.3, 2.71, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk51(book3) oscill result") ;
    gsl_test_rel(abserr,exp_abserr,1e-15,"qk51(book3) oscill abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk51(book3) oscill resabs") ;
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk51(book3) oscill resasc") ;
  }

  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result =-7.238969575482959717E-01;
    double exp_abserr = 1.286438572027470736E-14;
    double exp_resabs = 1.158720854723590099E+00;
    double exp_resasc = 9.270469641771273972E-01;

    alpha = 1.3 ;
    gsl_integration_qk61 (book3, 0.3, 2.71, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk61(book3) oscill result") ;
    gsl_test_rel(abserr,exp_abserr,1e-15,"qk61(book3) oscill abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk61(book3) oscill resabs") ;
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk61(book3) oscill resasc") ;
  }

  /* Test the non-adaptive gaussian integrator QNG */

  {
    int status = 0; size_t neval = 0 ;
    double result = 0, abserr = 0 ;
    double exp_result = 7.716049379303084599E-02;
    double exp_abserr = 9.424302194248481445E-08;
    int exp_neval  =  21;
    int exp_ier    =   0;

    alpha = 2.6 ;
    status = gsl_integration_qng (book1, 0.0, 1.0, 1e-1, 0.0,
				  &result, &abserr, &neval) ;
    gsl_test_rel(result,exp_result,1e-15,"qng(book1) smooth result") ;
    gsl_test_rel(abserr,exp_abserr,1e-15,"qng(book1) smooth abserr") ;
    gsl_test_int((int)neval,exp_neval,"qng(book1) smooth neval") ;  
    gsl_test_int(status,exp_ier,"qng(book1) smooth status") ;
  }

  {
    int status = 0; size_t neval = 0 ;
    double result = 0, abserr = 0 ;

    double exp_result = 7.716049382706505200E-02;
    double exp_abserr = 2.666891413688592065E-12;
    int exp_neval  =  43;
    int exp_ier    =   0;

    alpha = 2.6 ;
    status = gsl_integration_qng (book1, 0.0, 1.0, 0.0, 1e-9,
				  &result, &abserr, &neval) ;
    gsl_test_rel(result,exp_result,1e-15,"qng(book1) 43pt result") ;
    gsl_test_rel(abserr,exp_abserr,1e-15,"qng(book1) 43pt abserr") ;
    gsl_test_int((int)neval,exp_neval,"qng(book1) 43pt neval") ;  
    gsl_test_int(status,exp_ier,"qng(book1) 43pt status") ;
  }

  {
    int status = 0; size_t neval = 0 ;
    double result = 0, abserr = 0 ;

    double exp_result = 7.716049382716028138E-02;
    double exp_abserr = 8.566535680046932640E-16;
    int exp_neval  =  87;
    int exp_ier    =   0;

    alpha = 2.6 ;
    status = gsl_integration_qng (book1, 0.0, 1.0, 0.0, 1e-13,
				  &result, &abserr, &neval) ;
    gsl_test_rel(result,exp_result,1e-15,"qng(book1) 87pt result") ;
    gsl_test_rel(abserr,exp_abserr,1e-15,"qng(book1) 87pt abserr") ;
    gsl_test_int((int)neval,exp_neval,"qng(book1) 87pt neval") ;  
    gsl_test_int(status,exp_ier,"qng(book1) 87pt status") ;
  }

  {
    int status = 0; size_t neval = 0 ;
    double result = 0, abserr = 0 ;

    double exp_result = 3.222948711817264211E+01;
    double exp_abserr = 2.782360287710622515E+01;
    int exp_neval  =  87;
    int exp_ier    =  GSL_ETOL;

    alpha = -0.9 ;
    status = gsl_integration_qng (book1, 0.0, 1.0, 0.0, 1e-3,
				  &result, &abserr, &neval) ;
    gsl_test_rel(result,exp_result,1e-15,"qng(book1) beyond 87pt result") ;
    gsl_test_rel(abserr,exp_abserr,1e-15,"qng(book1) beyond 87pt abserr") ;
    gsl_test_int((int)neval,exp_neval,"qng(book1) beyond 87pt neval") ;  
    gsl_test_int(status,exp_ier,"qng(book1) beyond 87pt status") ;
  }

  {
    int status; size_t neval = 0 ;
    double result = 0, abserr = 0 ;
    double exp_result =-7.238969575482961938E-01;
    double exp_abserr = 1.277676889520056369E-14;
    int exp_neval  =  43;
    int exp_ier    =   0;

    alpha = 1.3 ;
    status = gsl_integration_qng (book3, 0.3, 2.71, 0.0, 1e-12,
				  &result, &abserr, &neval) ;
    gsl_test_rel(result,exp_result,1e-15,"qnq(book3) oscill result") ;
    gsl_test_rel(abserr,exp_abserr,1e-15,"qng(book3) oscill abserr") ;
    gsl_test_int((int)neval,exp_neval,"qng(book3) oscill neval") ;
    gsl_test_int(status,exp_ier,"qng(book3) oscill status") ;
  }




#ifdef JUNK
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

void
my_error_handler (const char *reason, const char *file, int line, int err)
{
  if (0) printf ("(caught [%s:%d: %s (%d)])\n", file, line, reason, err) ;
}
