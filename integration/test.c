#include <stdio.h>
#include <math.h>
#include <gsl_integration.h>
#include <gsl_errno.h>
#include <gsl_test.h>
#include <gsl_ieee_utils.h>

#include "tests.h"

double alpha = 0 ;

void my_error_handler (const char *reason, const char *file,
		       int line, int err);

int main (void)
{
  gsl_ieee_env_setup () ;

  gsl_set_error_handler (&my_error_handler); 

  /* Test the basic Gauss-Kronrod rules with a smooth positive function. */

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
    /* FIXME: we have some roundoff error on this one, but it seems
       reasonable (it's suprising that the others are exact) */
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
    /* FIXME: we have some roundoff error on this one, but it seems
       reasonable (it's suprising that the others are exact) */
    gsl_test_rel(abserr,exp_abserr,1e-5,"qk61(book1) smooth abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk61(book1) smooth resabs") ;    
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk61(book1) smooth resasc") ;
  }

  /* Now test the basic rules with a positive function that has a
     singularity. This should give large values of abserr which would
     find discrepancies in the abserr calculation. */

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

  /* Test the basic Gauss-Kronrod rules with a smooth oscillating
     function, over an unsymmetric range. This should find any
     discrepancies in the abscissae. */

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
    gsl_test_rel(result,exp_result,1e-15,"qng(book1) smooth 43pt result") ;
    gsl_test_rel(abserr,exp_abserr,1e-15,"qng(book1) smooth 43pt abserr") ;
    gsl_test_int((int)neval,exp_neval,"qng(book1) smooth 43pt neval") ;  
    gsl_test_int(status,exp_ier,"qng(book1) smooth 43pt status") ;
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
    gsl_test_rel(result,exp_result,1e-15,"qng(book1) 87pt smooth result") ;
    gsl_test_rel(abserr,exp_abserr,1e-15,"qng(book1) 87pt smooth abserr") ;
    gsl_test_int((int)neval,exp_neval,"qng(book1) 87pt smooth neval") ;  
    gsl_test_int(status,exp_ier,"qng(book1) 87pt smooth status") ;
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
    gsl_test_rel(result,exp_result,1e-15,"qng(book1) sing beyond 87pt result");
    gsl_test_rel(abserr,exp_abserr,1e-15,"qng(book1) sing beyond 87pt abserr");
    gsl_test_int((int)neval,exp_neval,"qng(book1) sing beyond 87pt neval") ;  
    gsl_test_int(status,exp_ier,"qng(book1) sing beyond 87pt status") ;
  }

  /* Test the adaptive integrator QAGE */

  {
    int status = 0, i; size_t last = 0,  neval = 0;
    double result = 0, abserr=0;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000) ;

    double exp_result = 7.716049382715854665E-02 ;
    double exp_abserr = 6.679384885865053037E-12 ;
    int exp_neval  =     165;
    int exp_ier    =       0;
    int exp_last   =       6;

    double a[6] = { 0, 0.5, 0.25, 0.125, 0.0625, 0.03125 } ;
    double b[6] = { 0.03125, 1, 0.5, 0.25, 0.125, 0.0625 } ;
    double r[6] = { 3.966769831709074375E-06, 5.491842501998222409E-02,
		    1.909827770934243926E-02, 2.776531175604360531E-03,
		    3.280661030752063693E-04, 3.522704932261797744E-05 } ;
    double e[6] = { 6.678528276336181873E-12, 6.097169993333454062E-16,
		    2.120334764359736934E-16, 3.082568839745514608E-17,
		    3.642265412331439511E-18, 3.910988124757650942E-19 } ;
    int iord[6] = { 1, 2, 3, 4, 5, 6 } ;

    alpha = 2.6 ;
    status = gsl_integration_qage (book1, 0.0, 1.0, 0.0, 1e-10, 
				   GSL_INTEG_GAUSS15, w, &last, 
				   &result, &abserr, &neval) ;

    gsl_test_rel(result,exp_result,1e-15,"qage(book1) smooth result") ;
    gsl_test_rel(abserr,exp_abserr,1e-6,"qage(book1) smooth abserr") ;
    gsl_test_int((int)neval,exp_neval,"qage(book1) smooth neval") ;  
    gsl_test_int((int)last,exp_last,"qage(book1) smooth last") ;  
    gsl_test_int(status,exp_ier,"qage(book1) smooth status") ;

    for (i = 0; i < 6 ; i++) 
	gsl_test_rel(w->alist[i],a[i],1e-15,"qage(book1) smooth alist") ;

    for (i = 0; i < 6 ; i++) 
	gsl_test_rel(w->blist[i],b[i],1e-15,"qage(book1) smooth blist") ;

    for (i = 0; i < 6 ; i++) 
	gsl_test_rel(w->rlist[i],r[i],1e-15,"qage(book1) smooth rlist") ;

    for (i = 0; i < 6 ; i++) 
	gsl_test_rel(w->elist[i],e[i],1e-6,"qage(book1) smooth elist") ;

    for (i = 0; i < 6 ; i++) 
	gsl_test_int((int)w->iord[i],iord[i]-1,"qage(book1) smooth iord") ;

    gsl_integration_workspace_free (w) ;

  }

  /* Test the same function using an absolute error bound and the
     21-point rule */

  {
    int status = 0, i; size_t last = 0,  neval = 0;
    double result = 0, abserr=0;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000) ;

    double exp_result = 7.716049382716050342E-02 ;
    double exp_abserr = 2.227969521869139532E-15 ;
    int exp_neval  =     315;
    int exp_ier    =       0;
    int exp_last   =       8;

    double a[8] = { 0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625,
		    0.0078125 } ;
    double b[8] = { 0.0078125, 1, 0.5, 0.25, 0.125, 0.0625, 0.03125,
		    0.015625 } ;
    double r[8] = { 3.696942726831556522E-08, 5.491842501998223103E-02,
		    1.909827770934243579E-02, 2.776531175604360097E-03,
		    3.280661030752062609E-04, 3.522704932261797744E-05,
		    3.579060884684503576E-06, 3.507395216921808047E-07 } ;
    double e[8] = { 1.371316364034059572E-15, 6.097169993333454062E-16,
		    2.120334764359736441E-16, 3.082568839745514608E-17,
		    3.642265412331439511E-18, 3.910988124757650460E-19,
		    3.973555800712018091E-20, 3.893990926286736620E-21 } ;
    int iord[8] = { 1, 2, 3, 4, 5, 6, 7, 8 } ;

    alpha = 2.6 ;
    status = gsl_integration_qage (book1, 0.0, 1.0, 1e-14, 0.0, 
				   GSL_INTEG_GAUSS21, w, &last, 
				   &result, &abserr, &neval) ;

    gsl_test_rel(result,exp_result,1e-15,"qage(book1,21pt) smooth result") ;
    gsl_test_rel(abserr,exp_abserr,1e-6,"qage(book1,21pt) smooth abserr") ;
    gsl_test_int((int)neval,exp_neval,"qage(book1,21pt) smooth neval") ;  
    gsl_test_int((int)last,exp_last,"qage(book1,21pt) smooth last") ;  
    gsl_test_int(status,exp_ier,"qage(book1,21pt) smooth status") ;

    for (i = 0; i < 8 ; i++) 
	gsl_test_rel(w->alist[i],a[i],1e-15,"qage(book1,21pt) smooth alist") ;

    for (i = 0; i < 8 ; i++) 
	gsl_test_rel(w->blist[i],b[i],1e-15,"qage(book1,21pt) smooth blist") ;

    for (i = 0; i < 8 ; i++) 
	gsl_test_rel(w->rlist[i],r[i],1e-15,"qage(book1,21pt) smooth rlist") ;

    for (i = 0; i < 8 ; i++) 
	gsl_test_rel(w->elist[i],e[i],1e-6,"qage(book1,21pt) smooth elist") ;

    for (i = 0; i < 8 ; i++) 
	gsl_test_int((int)w->iord[i],iord[i]-1,"qage(book1,21pt) smooth iord");

    gsl_integration_workspace_free (w) ;

  }

  /* Adaptive integration of an oscillatory function which terminates because
     of roundoff error, uses the 31-pt rule */

  {
    int status = 0, i; size_t last = 0,  neval = 0;
    double result = 0, abserr=0;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000) ;

    double exp_result = -7.238969575482960828E-01 ;
    double exp_abserr = 1.285717162868479170E-14 ;
    int exp_neval  =     403;
    int exp_ier    =     GSL_EROUND;
    int exp_last   =     7;

    double a[7] = { 1.504999999999999893,
		    1.203749999999999876,
		    2.107499999999999929,
		    2.999999999999999889E-01,
		    1.806249999999999911,
		    9.024999999999999689E-01,
		    1.655624999999999902 } ;
		    
    double b[7] = { 1.655624999999999902,
		    1.504999999999999893,
		    2.709999999999999964,
		    9.024999999999999689E-01,
		    2.107499999999999929,
		    1.203749999999999876,
		    1.806249999999999911 } ;
    
    double r[7] = { -1.169563901051563493E-01,
		    -2.210377785682320351E-01,
		    -3.014214297713070645E-02,
		    1.112947029779970676E-01,
		    -1.938177241117024219E-01,
		    -1.593317266741237914E-01,
		    -1.139058980899478740E-01 } ;

    double e[7] = { 1.298476771717864271E-15,
		    2.454012310784481733E-15,
		    1.819536846275027357E-15,
		    2.099789678945222176E-15,
		    2.151808998892583982E-15,
		    1.768937515068944431E-15,
		    1.264609507000667754E-15 } ;

    int iord[7] = { 2, 5, 4, 3, 6, 1, 7 } ;

    double exp_result2 = -7.238969575482960828E-01 ;
    double exp_abserr2 = 1.2858054644274593e-14 ;
    int exp_neval2  =     31;
    int exp_ier2    =     GSL_EROUND;
    int exp_last2   =     0;

    alpha = 1.3 ;
    status = gsl_integration_qage (book3, 0.3, 2.71, 1e-14, 0.0, 
				   GSL_INTEG_GAUSS31, w, &last, 
				   &result, &abserr, &neval) ;

    if ((int)neval == exp_neval) 
      {
	/* These are the results of running with an extended precision
	   fpu such as the pentium (the extra precision leads to a
	   different branch being taken in the code) */

	gsl_test_rel(result,exp_result,1e-15,"qage(book3,31pt) extended oscill result");
	gsl_test_rel(abserr,exp_abserr,1e-6,"qage(book3,31pt) extended oscill abserr") ;
	gsl_test_int((int)neval,exp_neval,"qage(book3,31pt) extended oscill neval") ;  
	gsl_test_int((int)last,exp_last,"qage(book3,31pt) extended oscill last") ;  
	gsl_test_int(status,exp_ier,"qage(book3,31pt) extended oscill status") ;
	
	for (i = 0; i < 7 ; i++) 
	  gsl_test_rel(w->alist[i],a[i],1e-15,"qage(book3,31pt) extended oscill alist");
	
	for (i = 0; i < 7 ; i++) 
	  gsl_test_rel(w->blist[i],b[i],1e-15,"qage(book3,31pt) extended oscill blist");
	
	for (i = 0; i < 7 ; i++) 
	  gsl_test_rel(w->rlist[i],r[i],1e-15,"qage(book3,31pt) extended oscill rlist");
	
	for (i = 0; i < 7 ; i++) 
	  gsl_test_rel(w->elist[i],e[i],1e-6,"qage(book3,31pt) extended oscill elist");
	
	for (i = 0; i < 7 ; i++) 
	  gsl_test_int((int)w->iord[i],iord[i]-1,"qage(book3,31pt) extended oscill iord");
      } 
    else
      {
	gsl_test_rel(result,exp_result2,1e-15,"qage(book3,31pt) oscill result");
	gsl_test_rel(abserr,exp_abserr2,1e-6,"qage(book3,31pt) oscill abserr");
	gsl_test_int((int)neval,exp_neval2,"qage(book3,31pt) oscill neval") ;  
	gsl_test_int((int)last,exp_last2,"qage(book3,31pt) oscill last") ;  
	gsl_test_int(status,exp_ier2,"qage(book3,31pt) oscill status") ;
      }

    gsl_integration_workspace_free (w) ;

  }

  /* Check the singularity detection (singularity at x=-0.1 in this example) */

  {
    int status = 0, i; size_t last = 0,  neval = 0;
    double result = 0, abserr=0;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000) ;

    double exp_result =  2.493668291024357500E+15 ;
    double exp_abserr =  4.062997483865515000E+15;
    int exp_neval  =     5151;
    int exp_ier    =     GSL_ESING;
    int exp_last   =     51;

    double a[10] = { -1.000000000000014211E-01,
		     0.000000000000000000,
		     -1.000000000000000000,
		     -5.000000000000000000E-01,
		     -2.500000000000000000E-01,
		     -6.250000000000000000E-02,
		     -9.375000000000000000E-02,
		     -1.250000000000000000E-01,
		     -1.093750000000000000E-01,
		     -9.765625000000000000E-02 } ;
    
    double b[10] = { -9.999999999999964473E-02,
		     1.000000000000000000,
		     -5.000000000000000000E-01,
		     -2.500000000000000000E-01,
		     -1.250000000000000000E-01,
		     0.000000000000000000,
		     -6.250000000000000000E-02,
		     -1.093750000000000000E-01,
		     -1.015625000000000000E-01,
		     -9.375000000000000000E-02 } ;
    
    double r[10] = { 2.458497689591106500E+15,
		     9.090909090909089774E-02,
		     1.388888888888888985E-02,
		     4.166666666666666435E-02,
		     3.333333333333332593E-01,
		     1.666666666666666574E-01,
		     1.333333333333333703,
		     6.666666666666665186E-01,
		     5.333333333333337478,
		     2.666666666666667407 } ;

    double e[10] = { 4.060351180997878500E+15,
		    1.009293658750142202E-15,
		    1.541976423090495387E-16,
		    4.625929269271484928E-16,
		    3.700743415417187942E-15,
		    1.850371707708593971E-15,
		    1.480297366166876124E-14,
		    7.401486830834375884E-15,
		    5.921189464667505756E-14,
		    2.960594732333751616E-14 } ;

    int iord[10] = { 1, 51, 49, 50, 47, 48, 45, 43, 46, 44 } ;

    alpha = 1.0 ;
    status = gsl_integration_qage (book16, -1.0, 1.0, 1e-14, 0.0, 
				   GSL_INTEG_GAUSS51, w, &last, 
				   &result, &abserr, &neval) ;

    gsl_test_rel(result,exp_result,1e-15,"qage(book16,51pt) sing result") ;
    gsl_test_rel(abserr,exp_abserr,1e-6,"qage(book16,51pt) sing abserr") ;
    gsl_test_int((int)neval,exp_neval,"qage(book16,51pt) sing neval") ;  
    gsl_test_int((int)last,exp_last,"qage(book16,51pt) sing last") ;  
    gsl_test_int(status,exp_ier,"qage(book16,51pt) sing status") ;

    for (i = 0; i < 10 ; i++) 
	gsl_test_rel(w->alist[i],a[i],1e-15,"qage(book16,51pt) sing alist") ;

    for (i = 0; i < 10 ; i++) 
	gsl_test_rel(w->blist[i],b[i],1e-15,"qage(book16,51pt) sing blist") ;

    for (i = 0; i < 10 ; i++) 
	gsl_test_rel(w->rlist[i],r[i],1e-15,"qage(book16,51pt) sing rlist") ;

    for (i = 0; i < 10 ; i++) 
	gsl_test_rel(w->elist[i],e[i],1e-6,"qage(book16,51pt) sing elist") ;

    for (i = 0; i < 10 ; i++) 
	gsl_test_int((int)w->iord[i],iord[i]-1,"qage(book16,51pt) iord") ;

    gsl_integration_workspace_free (w) ;

  }

  /* Check for hitting the iteration limit */

  {
    int status = 0, i; size_t last = 0,  neval = 0;
    double result = 0, abserr=0;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (3) ;

    double exp_result =  9.565151449233894709 ;
    double exp_abserr =  1.570369823891028460E+01;
    int exp_neval  =     305;
    int exp_ier    =     GSL_EMAXITER;
    int exp_last   =     3;

    double a[3] = { -5.000000000000000000E-01,
		    0.000000000000000000,
		    -1.000000000000000000 } ;
    
    double b[3] = { 0.000000000000000000,
                    1.000000000000000000,
                    -5.000000000000000000E-01 } ;
    
    double r[3] = { 9.460353469435913709,
                    9.090909090909091161E-02,
                    1.388888888888888812E-02 } ;
    
    double e[3] = { 1.570369823891028460E+01,
                    1.009293658750142399E-15,
                    1.541976423090495140E-16 } ;
    
    int iord[3] = { 1, 2, 3 } ;
    
    alpha = 1.0 ;
    
    status = gsl_integration_qage (book16, -1.0, 1.0, 1e-14, 0.0, 
				   GSL_INTEG_GAUSS61, w, &last, 
				   &result, &abserr, &neval) ;

    gsl_test_rel(result,exp_result,1e-15,"qage(book16,61pt) limit result") ;
    gsl_test_rel(abserr,exp_abserr,1e-6,"qage(book16,61pt) limit abserr") ;
    gsl_test_int((int)neval,exp_neval,"qage(book16,61pt) limit neval") ;  
    gsl_test_int((int)last,exp_last,"qage(book16,61pt) limit last") ;  
    gsl_test_int(status,exp_ier,"qage(book16,61pt) limit status") ;

    for (i = 0; i < 3 ; i++) 
	gsl_test_rel(w->alist[i],a[i],1e-15,"qage(book16,61pt) limit alist") ;

    for (i = 0; i < 3 ; i++) 
	gsl_test_rel(w->blist[i],b[i],1e-15,"qage(book16,61pt) limit blist") ;

    for (i = 0; i < 3 ; i++) 
	gsl_test_rel(w->rlist[i],r[i],1e-15,"qage(book16,61pt) limit rlist") ;

    for (i = 0; i < 3 ; i++) 
	gsl_test_rel(w->elist[i],e[i],1e-6,"qage(book16,61pt) limit elist") ;

    for (i = 0; i < 3 ; i++) 
	gsl_test_int((int)w->iord[i],iord[i]-1,"qage(book16,61pt) limit iord");

    gsl_integration_workspace_free (w) ;

  }

  return gsl_test_summary() ;
} 

void
my_error_handler (const char *reason, const char *file, int line, int err)
{
  if (0) printf ("(caught [%s:%d: %s (%d)])\n", file, line, reason, err) ;
}
