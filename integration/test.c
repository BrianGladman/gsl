#include <config.h>
#include <stdio.h>
#include <math.h>
#include <gsl_integration.h>
#include <gsl_errno.h>
#include <gsl_test.h>
#include <gsl_ieee_utils.h>

#include "tests.h"

void my_error_handler (const char *reason, const char *file,
		       int line, int err);

int main (void)
{
  gsl_ieee_env_setup ();
  gsl_set_error_handler (&my_error_handler); 

  /* Test the basic Gauss-Kronrod rules with a smooth positive function. */

  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result = 7.716049357767090777E-02;
    double exp_abserr = 2.990224871000550874E-06;
    double exp_resabs = 7.716049357767090777E-02;
    double exp_resasc = 4.434273814139995384E-02;

    double alpha = 2.6 ;
    gsl_function f = { &f1, &alpha } ;

    gsl_integration_qk15 (&f, 0.0, 1.0, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk15(f1) smooth result") ;
    gsl_test_rel(abserr,exp_abserr,1e-7,"qk15(f1) smooth abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk15(f1) smooth resabs") ;    
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk15(f1) smooth resasc") ;
  }

  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result = 7.716049379303084599E-02;
    double exp_abserr = 9.424302194248481445E-08;
    double exp_resabs = 7.716049379303084599E-02;
    double exp_resasc = 4.434311425038358484E-02;

    double alpha = 2.6 ;
    gsl_function f = { &f1, &alpha } ;

    gsl_integration_qk21 (&f, 0.0, 1.0, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk21(f1) smooth result") ;
    gsl_test_rel(abserr,exp_abserr,1e-7,"qk21(f1) smooth abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk21(f1) smooth resabs") ;    
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk21(f1) smooth resasc") ;
  }

  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result = 7.716049382494900855E-02;
    double exp_abserr = 1.713503193600029893E-09;
    double exp_resabs = 7.716049382494900855E-02;
    double exp_resasc = 4.427995051868838933E-02;

    double alpha = 2.6 ;
    gsl_function f = { &f1, &alpha } ;

    gsl_integration_qk31 (&f, 0.0, 1.0, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk31(f1) smooth result") ;
    gsl_test_rel(abserr,exp_abserr,1e-7,"qk31(f1) smooth abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk31(f1) smooth resabs") ;    
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk31(f1) smooth resasc") ;
  }

  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result = 7.716049382681375302E-02;
    double exp_abserr = 9.576386660975511224E-11;
    double exp_resabs = 7.716049382681375302E-02;
    double exp_resasc = 4.421521169637691873E-02;

    double alpha = 2.6 ;
    gsl_function f = { &f1, &alpha } ;

    gsl_integration_qk41 (&f, 0.0, 1.0, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk41(f1) smooth result") ;
    gsl_test_rel(abserr,exp_abserr,1e-7,"qk41(f1) smooth abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk41(f1) smooth resabs") ;    
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk41(f1) smooth resasc") ;
  }

  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result = 7.716049382708510540E-02;
    double exp_abserr = 1.002079980317363772E-11;
    double exp_resabs = 7.716049382708510540E-02;
    double exp_resasc = 4.416474291216854892E-02;

    double alpha = 2.6 ;
    gsl_function f = { &f1, &alpha } ;

    gsl_integration_qk51 (&f, 0.0, 1.0, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk51(f1) smooth result") ;
    gsl_test_rel(abserr,exp_abserr,1e-5,"qk51(f1) smooth abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk51(f1) smooth resabs") ;    
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk51(f1) smooth resasc") ;
  }

  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result = 7.716049382713800753E-02;
    double exp_abserr = 1.566060362296155616E-12;
    double exp_resabs = 7.716049382713800753E-02;
    double exp_resasc = 4.419287685934316506E-02;

    double alpha = 2.6 ;
    gsl_function f = { &f1, &alpha } ;

    gsl_integration_qk61 (&f, 0.0, 1.0, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk61(f1) smooth result") ;
    gsl_test_rel(abserr,exp_abserr,1e-5,"qk61(f1) smooth abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk61(f1) smooth resabs") ;    
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk61(f1) smooth resasc") ;
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

    double alpha = -0.9 ;
    gsl_function f = { &f1, &alpha } ;

    gsl_integration_qk15 (&f, 0.0, 1.0, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk15(f1) singular result") ;
    gsl_test_rel(abserr,exp_abserr,1e-7,"qk15(f1) singular abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk15(f1) singular resabs") ;    
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk15(f1) singular resasc") ;
  }

  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result = 1.799045317938126232E+01;
    double exp_abserr = 2.782360287710622515E+01;
    double exp_resabs = 1.799045317938126232E+01;
    double exp_resasc = 2.782360287710622515E+01;

    double alpha = -0.9 ;
    gsl_function f = { &f1, &alpha } ;

    gsl_integration_qk21 (&f, 0.0, 1.0, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk21(f1) singular result") ;
    gsl_test_rel(abserr,exp_abserr,1e-7,"qk21(f1) singular abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk21(f1) singular resabs") ;    
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk21(f1) singular resasc") ;
  }

  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result = 2.081873305159121657E+01;
    double exp_abserr = 3.296500137482590276E+01;
    double exp_resabs = 2.081873305159121301E+01;
    double exp_resasc = 3.296500137482590276E+01;

    double alpha = -0.9 ;
    gsl_function f = { &f1, &alpha } ;

    gsl_integration_qk31 (&f, 0.0, 1.0, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk31(f1) singular result") ;
    gsl_test_rel(abserr,exp_abserr,1e-7,"qk31(f1) singular abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk31(f1) singular resabs") ;    
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk31(f1) singular resasc") ;
  }

  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result = 2.288677623903126701E+01;
    double exp_abserr = 3.671538820274916048E+01;
    double exp_resabs = 2.288677623903126701E+01;
    double exp_resasc = 3.671538820274916048E+01;

    double alpha = -0.9 ;
    gsl_function f = { &f1, &alpha } ;

    gsl_integration_qk41 (&f, 0.0, 1.0, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk41(f1) singular result") ;
    gsl_test_rel(abserr,exp_abserr,1e-7,"qk41(f1) singular abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk41(f1) singular resabs") ;    
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk41(f1) singular resasc") ;
  }

  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result = 2.449953612016972215E+01;
    double exp_abserr = 3.967771249391228849E+01;
    double exp_resabs = 2.449953612016972215E+01;
    double exp_resasc = 3.967771249391228849E+01;

    double alpha = -0.9 ;
    gsl_function f = { &f1, &alpha } ;

    gsl_integration_qk51 (&f, 0.0, 1.0, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk51(f1) singular result") ;
    gsl_test_rel(abserr,exp_abserr,1e-7,"qk51(f1) singular abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk51(f1) singular resabs") ;    
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk51(f1) singular resasc") ;
  }

  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result = 2.583030240976628988E+01;
    double exp_abserr = 4.213750493076978643E+01;
    double exp_resabs = 2.583030240976628988E+01;
    double exp_resasc = 4.213750493076978643E+01;

    double alpha = -0.9 ;
    gsl_function f = { &f1, &alpha } ;

    gsl_integration_qk61 (&f, 0.0, 1.0, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk61(f1) singular result") ;
    gsl_test_rel(abserr,exp_abserr,1e-7,"qk61(f1) singular abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk61(f1) singular resabs") ;    
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk61(f1) singular resasc") ;
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

    double alpha = 1.3 ;
    gsl_function f = { &f3, &alpha } ;

    gsl_integration_qk15 (&f, 0.3, 2.71, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk15(f3) oscill result") ;
    gsl_test_rel(abserr,exp_abserr,1e-7,"qk15(f3) oscill abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk15(f3) oscill resabs") ;
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk15(f3) oscill resasc") ;
  }

  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result =-7.238969575482959717E-01;
    double exp_abserr = 7.999213141433641888E-11;
    double exp_resabs = 1.150829032708484023E+00;
    double exp_resasc = 9.297591249133687619E-01;

    double alpha = 1.3 ;
    gsl_function f = { &f3, &alpha } ;
    
    gsl_integration_qk21 (&f, 0.3, 2.71, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk21(f3) oscill result") ;
    gsl_test_rel(abserr,exp_abserr,1e-5,"qk21(f3) oscill abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk21(f3) oscill resabs") ;
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk21(f3) oscill resasc") ;
  }

  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result =-7.238969575482959717E-01;
    double exp_abserr = 1.285805464427459261E-14;
    double exp_resabs = 1.158150602093290571E+00;
    double exp_resasc = 9.277828092501518853E-01;

    double alpha = 1.3 ;
    gsl_function f = { &f3, &alpha } ;

    gsl_integration_qk31 (&f, 0.3, 2.71, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk31(f3) oscill result") ;
    gsl_test_rel(abserr,exp_abserr,1e-7,"qk31(f3) oscill abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk31(f3) oscill resabs") ;
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk31(f3) oscill resasc") ;
  }

  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result =-7.238969575482959717E-01;
    double exp_abserr = 1.286535726271015626E-14;
    double exp_resabs = 1.158808363486595328E+00;
    double exp_resasc = 9.264382258645686985E-01;

    double alpha = 1.3 ;
    gsl_function f = { &f3, &alpha } ;

    gsl_integration_qk41 (&f, 0.3, 2.71, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk41(f3) oscill result") ;
    gsl_test_rel(abserr,exp_abserr,1e-7,"qk41(f3) oscill abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk41(f3) oscill resabs") ;
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk41(f3) oscill resasc") ;
  }

  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result =-7.238969575482961938E-01;
    double exp_abserr = 1.285290995039385778E-14;
    double exp_resabs = 1.157687209264406381E+00;
    double exp_resasc = 9.264666884071264263E-01;

    double alpha = 1.3 ;
    gsl_function f = { &f3, &alpha } ;

    gsl_integration_qk51 (&f, 0.3, 2.71, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk51(f3) oscill result") ;
    gsl_test_rel(abserr,exp_abserr,1e-7,"qk51(f3) oscill abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk51(f3) oscill resabs") ;
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk51(f3) oscill resasc") ;
  }

  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    double exp_result =-7.238969575482959717E-01;
    double exp_abserr = 1.286438572027470736E-14;
    double exp_resabs = 1.158720854723590099E+00;
    double exp_resasc = 9.270469641771273972E-01;

    double alpha = 1.3 ;
    gsl_function f = { &f3, &alpha } ;

    gsl_integration_qk61 (&f, 0.3, 2.71, 
			  &result, &abserr, &resabs, &resasc) ;
    gsl_test_rel(result,exp_result,1e-15,"qk61(f3) oscill result") ;
    gsl_test_rel(abserr,exp_abserr,1e-7,"qk61(f3) oscill abserr") ;
    gsl_test_rel(resabs,exp_resabs,1e-15,"qk61(f3) oscill resabs") ;
    gsl_test_rel(resasc,exp_resasc,1e-15,"qk61(f3) oscill resasc") ;
  }

  /* Test the non-adaptive gaussian integrator QNG */

  {
    int status = 0; size_t neval = 0 ;
    double result = 0, abserr = 0 ;
    double exp_result = 7.716049379303083211E-02;
    double exp_abserr = 9.424302199601294244E-08;
    int exp_neval  =  21;
    int exp_ier    =   0;

    double alpha = 2.6 ;
    gsl_function f = { &f1, &alpha } ;

    status = gsl_integration_qng (&f, 0.0, 1.0, 1e-1, 0.0,
				  &result, &abserr, &neval) ;
    gsl_test_rel(result,exp_result,1e-15,"qng(f1) smooth result") ;
    gsl_test_rel(abserr,exp_abserr,1e-7,"qng(f1) smooth abserr") ;
    gsl_test_int((int)neval,exp_neval,"qng(f1) smooth neval") ;  
    gsl_test_int(status,exp_ier,"qng(f1) smooth status") ;
  }

  {
    int status = 0; size_t neval = 0 ;
    double result = 0, abserr = 0 ;

    double exp_result = 7.716049382706505200E-02;
    double exp_abserr = 2.666893044866214501E-12;
    int exp_neval  =  43;
    int exp_ier    =   0;

    double alpha = 2.6 ;
    gsl_function f = { &f1, &alpha } ;

    status = gsl_integration_qng (&f, 0.0, 1.0, 0.0, 1e-9,
				  &result, &abserr, &neval) ;
    gsl_test_rel(result,exp_result,1e-15,"qng(f1) smooth 43pt result") ;
    gsl_test_rel(abserr,exp_abserr,1e-6,"qng(f1) smooth 43pt abserr") ;
    gsl_test_int((int)neval,exp_neval,"qng(f1) smooth 43pt neval") ;  
    gsl_test_int(status,exp_ier,"qng(f1) smooth 43pt status") ;
  }

  {
    int status; size_t neval = 0 ;
    double result = 0, abserr = 0 ;
    double exp_result =-7.238969575482961938E-01;
    double exp_abserr = 1.277676889520056369E-14;
    int exp_neval  =  43;
    int exp_ier    =   0;

    double alpha = 1.3 ;
    gsl_function f = { &f3, &alpha } ;

    status = gsl_integration_qng (&f, 0.3, 2.71, 0.0, 1e-12,
				  &result, &abserr, &neval) ;
    gsl_test_rel(result,exp_result,1e-15,"qnq(f3) oscill result") ;
    gsl_test_rel(abserr,exp_abserr,1e-7,"qng(f3) oscill abserr") ;
    gsl_test_int((int)neval,exp_neval,"qng(f3) oscill neval") ;
    gsl_test_int(status,exp_ier,"qng(f3) oscill status") ;
  }

  {
    int status = 0; size_t neval = 0 ;
    double result = 0, abserr = 0 ;

    double exp_result = 7.716049382716029525E-02;
    double exp_abserr = 8.566535680046930668E-16;
    int exp_neval  =  87;
    int exp_ier    =   0;

    double alpha = 2.6 ;
    gsl_function f = { &f1, &alpha } ;

    status = gsl_integration_qng (&f, 0.0, 1.0, 0.0, 1e-13,
				  &result, &abserr, &neval) ;
    gsl_test_rel(result,exp_result,1e-15,"qng(f1) 87pt smooth result") ;
    gsl_test_rel(abserr,exp_abserr,1e-7,"qng(f1) 87pt smooth abserr") ;
    gsl_test_int((int)neval,exp_neval,"qng(f1) 87pt smooth neval") ;  
    gsl_test_int(status,exp_ier,"qng(f1) 87pt smooth status") ;
  }

  {
    int status = 0; size_t neval = 0 ;
    double result = 0, abserr = 0 ;

    double exp_result = 3.222948711817264211E+01;
    double exp_abserr = 2.782360287710622870E+01;
    int exp_neval  =  87;
    int exp_ier    =  GSL_ETOL;

    double alpha = -0.9 ;
    gsl_function f = { &f1, &alpha } ;

    status = gsl_integration_qng (&f, 0.0, 1.0, 0.0, 1e-3,
				  &result, &abserr, &neval) ;
    gsl_test_rel(result,exp_result,1e-15,"qng(f1) sing beyond 87pt result");
    gsl_test_rel(abserr,exp_abserr,1e-7,"qng(f1) sing beyond 87pt abserr");
    gsl_test_int((int)neval,exp_neval,"qng(f1) sing beyond 87pt neval") ;  
    gsl_test_int(status,exp_ier,"qng(f1) sing beyond 87pt status") ;
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

    double alpha = 2.6 ;
    gsl_function f = { &f1, &alpha } ;

    status = gsl_integration_qage (&f, 0.0, 1.0, 0.0, 1e-10, 
				   GSL_INTEG_GAUSS15, w, &last, 
				   &result, &abserr, &neval) ;

    gsl_test_rel(result,exp_result,1e-15,"qage(f1) smooth result") ;
    gsl_test_rel(abserr,exp_abserr,1e-6,"qage(f1) smooth abserr") ;
    gsl_test_int((int)neval,exp_neval,"qage(f1) smooth neval") ;  
    gsl_test_int((int)last,exp_last,"qage(f1) smooth last") ;  
    gsl_test_int(status,exp_ier,"qage(f1) smooth status") ;

    for (i = 0; i < 6 ; i++) 
	gsl_test_rel(w->alist[i],a[i],1e-15,"qage(f1) smooth alist") ;

    for (i = 0; i < 6 ; i++) 
	gsl_test_rel(w->blist[i],b[i],1e-15,"qage(f1) smooth blist") ;

    for (i = 0; i < 6 ; i++) 
	gsl_test_rel(w->rlist[i],r[i],1e-15,"qage(f1) smooth rlist") ;

    for (i = 0; i < 6 ; i++) 
	gsl_test_rel(w->elist[i],e[i],1e-6,"qage(f1) smooth elist") ;

    for (i = 0; i < 6 ; i++) 
	gsl_test_int((int)w->iord[i],iord[i]-1,"qage(f1) smooth iord") ;

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

    double alpha = 2.6 ;
    gsl_function f = { &f1, &alpha } ;

    status = gsl_integration_qage (&f, 0.0, 1.0, 1e-14, 0.0, 
				   GSL_INTEG_GAUSS21, w, &last, 
				   &result, &abserr, &neval) ;

    gsl_test_rel(result,exp_result,1e-15,"qage(f1,21pt) smooth result") ;
    gsl_test_rel(abserr,exp_abserr,1e-6,"qage(f1,21pt) smooth abserr") ;
    gsl_test_int((int)neval,exp_neval,"qage(f1,21pt) smooth neval") ;  
    gsl_test_int((int)last,exp_last,"qage(f1,21pt) smooth last") ;  
    gsl_test_int(status,exp_ier,"qage(f1,21pt) smooth status") ;

    for (i = 0; i < 8 ; i++) 
	gsl_test_rel(w->alist[i],a[i],1e-15,"qage(f1,21pt) smooth alist") ;

    for (i = 0; i < 8 ; i++) 
	gsl_test_rel(w->blist[i],b[i],1e-15,"qage(f1,21pt) smooth blist") ;

    for (i = 0; i < 8 ; i++) 
	gsl_test_rel(w->rlist[i],r[i],1e-15,"qage(f1,21pt) smooth rlist") ;

    for (i = 0; i < 8 ; i++) 
	gsl_test_rel(w->elist[i],e[i],1e-6,"qage(f1,21pt) smooth elist") ;

    for (i = 0; i < 8 ; i++) 
	gsl_test_int((int)w->iord[i],iord[i]-1,"qage(f1,21pt) smooth iord");

    gsl_integration_workspace_free (w) ;

  }

  /* Adaptive integration of an oscillatory function which terminates because
     of roundoff error, uses the 31-pt rule */

  {
    int status = 0; size_t last = 0,  neval = 0;
    double result = 0, abserr=0;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000) ;

    double exp_result = -7.238969575482960828E-01 ;
    double exp_abserr =  1.2858054644274593e-14 ;
    int exp_neval   =     31;
    int exp_ier     =     GSL_EROUND;
    int exp_last    =     0;

    double alpha = 1.3 ;
    gsl_function f = { &f3, &alpha } ;

    status = gsl_integration_qage (&f, 0.3, 2.71, 1e-14, 0.0, 
				   GSL_INTEG_GAUSS31, w, &last, 
				   &result, &abserr, &neval) ;

    gsl_test_rel(result,exp_result,1e-15,"qage(f3,31pt) oscill result");
    gsl_test_rel(abserr,exp_abserr,1e-6,"qage(f3,31pt) oscill abserr");
    gsl_test_int((int)neval,exp_neval,"qage(f3,31pt) oscill neval") ;  
    gsl_test_int((int)last,exp_last,"qage(f3,31pt) oscill last") ;  
    gsl_test_int(status,exp_ier,"qage(f3,31pt) oscill status") ;

    gsl_integration_workspace_free (w) ;

  }

  /* Check the singularity detection (singularity at x=-0.1 in this example) */

  {
    int status = 0; size_t last = 0,  neval = 0;
    double result = 0, abserr=0;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000) ;

    int exp_neval  =     5151;
    int exp_ier    =     GSL_ESING;
    int exp_last   =     51;

    double alpha = 2.0 ;
    gsl_function f = { &f16, &alpha } ;

    status = gsl_integration_qage (&f, -1.0, 1.0, 1e-14, 0.0, 
				   GSL_INTEG_GAUSS51, w, &last, 
				   &result, &abserr, &neval) ;

    gsl_test_int((int)neval,exp_neval,"qage(f16,51pt) sing neval") ;  
    gsl_test_int((int)last,exp_last,"qage(f16,51pt) sing last") ;  
    gsl_test_int(status,exp_ier,"qage(f16,51pt) sing status") ;

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

    double alpha = 1.0 ;
    gsl_function f = { &f16, &alpha } ;
    
    status = gsl_integration_qage (&f, -1.0, 1.0, 1e-14, 0.0, 
				   GSL_INTEG_GAUSS61, w, &last, 
				   &result, &abserr, &neval) ;

    gsl_test_rel(result,exp_result,1e-15,"qage(f16,61pt) limit result") ;
    gsl_test_rel(abserr,exp_abserr,1e-6,"qage(f16,61pt) limit abserr") ;
    gsl_test_int((int)neval,exp_neval,"qage(f16,61pt) limit neval") ;  
    gsl_test_int((int)last,exp_last,"qage(f16,61pt) limit last") ;  
    gsl_test_int(status,exp_ier,"qage(f16,61pt) limit status") ;

    for (i = 0; i < 3 ; i++) 
	gsl_test_rel(w->alist[i],a[i],1e-15,"qage(f16,61pt) limit alist") ;

    for (i = 0; i < 3 ; i++) 
	gsl_test_rel(w->blist[i],b[i],1e-15,"qage(f16,61pt) limit blist") ;

    for (i = 0; i < 3 ; i++) 
	gsl_test_rel(w->rlist[i],r[i],1e-15,"qage(f16,61pt) limit rlist") ;

    for (i = 0; i < 3 ; i++) 
	gsl_test_rel(w->elist[i],e[i],1e-6,"qage(f16,61pt) limit elist") ;

    for (i = 0; i < 3 ; i++) 
	gsl_test_int((int)w->iord[i],iord[i]-1,"qage(f16,61pt) limit iord");

    gsl_integration_workspace_free (w) ;

  }

  /* Test the adaptive integrator with extrapolation QAGSE */

  {
    int status = 0, i; size_t last = 0,  neval = 0;
    double result = 0, abserr=0;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000) ;

    double exp_result = 7.716049382715789440E-02 ;
    double exp_abserr = 2.216394961010438404E-12 ;
    int exp_neval  =     189;
    int exp_ier    =       0;
    int exp_last   =       5;

    double a[5] = { 0, 0.5, 0.25, 0.125, 0.0625 } ;
    double b[5] = { 0.0625, 1, 0.5, 0.25, 0.125 } ;
    double r[5] = { 3.919381915366914693E-05,
		    5.491842501998223103E-02,
		    1.909827770934243579E-02,
		    2.776531175604360097E-03,
		    3.280661030752062609E-04 } ;
    double e[5] = { 2.215538742580964735E-12,
		    6.097169993333454062E-16,
		    2.120334764359736441E-16,
		    3.082568839745514608E-17,
		    3.642265412331439511E-18 } ;
    int iord[5] = { 1, 2, 3, 4, 5 } ;

    double alpha = 2.6 ;
    gsl_function f = { &f1, &alpha } ;

    status = gsl_integration_qagse (&f, 0.0, 1.0, 0.0, 1e-10, 
				    w, &last, 
				    &result, &abserr, &neval) ;

    gsl_test_rel(result,exp_result,1e-15,"qagse(f1) smooth result") ;
    gsl_test_rel(abserr,exp_abserr,1e-6,"qagse(f1) smooth abserr") ;
    gsl_test_int((int)neval,exp_neval,"qagse(f1) smooth neval") ;  
    gsl_test_int((int)last,exp_last,"qagse(f1) smooth last") ;  
    gsl_test_int(status,exp_ier,"qagse(f1) smooth status") ;

    for (i = 0; i < 5 ; i++) 
	gsl_test_rel(w->alist[i],a[i],1e-15,"qagse(f1) smooth alist") ;

    for (i = 0; i < 5 ; i++) 
	gsl_test_rel(w->blist[i],b[i],1e-15,"qagse(f1) smooth blist") ;

    for (i = 0; i < 5 ; i++) 
	gsl_test_rel(w->rlist[i],r[i],1e-15,"qagse(f1) smooth rlist") ;

    for (i = 0; i < 5 ; i++) 
	gsl_test_rel(w->elist[i],e[i],1e-6,"qagse(f1) smooth elist") ;

    for (i = 0; i < 5 ; i++) 
	gsl_test_int((int)w->iord[i],iord[i]-1,"qagse(f1) smooth iord") ;

    gsl_integration_workspace_free (w) ;

  }

  /* Test f11 using an absolute error bound */

  {
    int status = 0, i; size_t last = 0,  neval = 0;
    double result = 0, abserr=0;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000) ;

    /* All results are for GSL_IEEE_MODE=double-precision */

    double exp_result = -5.908755278982136588E+03 ;
    double exp_abserr = 1.299646281053874554E-10 ;
    int exp_neval  =     357;
    int exp_ier    =       0;
    int exp_last   =       9;

    double a[9] = { 1.000000000000000000E+00,
		    5.005000000000000000E+02,
		    2.507500000000000000E+02,
		    1.258750000000000000E+02,
		    6.343750000000000000E+01,
		    3.221875000000000000E+01,
		    1.660937500000000000E+01,
		    8.804687500000000000E+00,
		    4.902343750000000000E+00 } ;
    double b[9] = { 4.902343750000000000E+00,
		    1.000000000000000000E+03,
		    5.005000000000000000E+02,
		    2.507500000000000000E+02,
		    1.258750000000000000E+02,
		    6.343750000000000000E+01,
		    3.221875000000000000E+01,
		    1.660937500000000000E+01,
		    8.804687500000000000E+00 } ;
    double r[9] = { -3.890977835520834649E+00,
		    -3.297343675805121620E+03,
		    -1.475904154146372775E+03,
		    -6.517404019686431411E+02,
		    -2.829354222635842007E+02,
		    -1.201692001973227519E+02,
		    -4.959999906099650246E+01,
		    -1.971441499411640308E+01,
		    -7.457032710459004399E+00 } ;
    double e[9] = { 6.448276035006137169E-11,
		    3.660786868980994028E-11,
		    1.638582774073219226E-11,
		    7.235772003440423011E-12,
		    3.141214202790722909E-12,
		    1.334146129098576244E-12,
		    5.506706097890446534E-13,
		    2.188739744348345039E-13,
		    8.278969410534525339E-14 } ;
    int iord[9] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 } ;

    double alpha = 2.0 ;
    gsl_function f = { &f11, &alpha } ;

    status = gsl_integration_qagse (&f, 1.0, 1000.0, 1e-7, 0.0, 
				    w, &last, 
				    &result, &abserr, &neval) ;

    gsl_test_rel(result,exp_result,1e-15,"qagse(f11) smooth result") ;
    gsl_test_rel(abserr,exp_abserr,1e-3,"qagse(f11) smooth abserr") ;
    gsl_test_int((int)neval,exp_neval,"qagse(f11) smooth neval") ;  
    gsl_test_int((int)last,exp_last,"qagse(f11) smooth last") ;  
    gsl_test_int(status,exp_ier,"qagse(f11) smooth status") ;

    for (i = 0; i < 9 ; i++) 
	gsl_test_rel(w->alist[i],a[i],1e-15,"qagse(f11) smooth alist") ;

    for (i = 0; i < 9 ; i++) 
	gsl_test_rel(w->blist[i],b[i],1e-15,"qagse(f11) smooth blist") ;

    for (i = 0; i < 9 ; i++) 
	gsl_test_rel(w->rlist[i],r[i],1e-15,"qagse(f11) smooth rlist") ;

    for (i = 0; i < 9 ; i++) 
	gsl_test_rel(w->elist[i],e[i],1e-5,"qagse(f11) smooth elist") ;

    for (i = 0; i < 9 ; i++) 
	gsl_test_int((int)w->iord[i],iord[i]-1,"qagse(f11) smooth iord");

    gsl_integration_workspace_free (w) ;

  }


  return gsl_test_summary() ;
} 

void
my_error_handler (const char *reason, const char *file, int line, int err)
{
  if (0) printf ("(caught [%s:%d: %s (%d)])\n", file, line, reason, err) ;
}
