#include <math.h>

int
gsl_integration_qk15 (double (*f) (double x),
		      double a, double b,
		      double * result, double * abserr,
		      double * resabs, double * resasc)
{
  /* Computes I=Integral of F over (a,b), with error estimate and
     J=integral of fabs(f) over (a,b), using 15-point gauss-kronrod rules

     Reimplementation of public domain SLATEC/QUADPACK (original
     authors, Robert Piessens and Elise de Doncker)

     Input parameters, f - function defining the integrand, a - lower
     limit of integration, b - upper limit of integration

     Returned parameters,

     result - approximation to the integral I by applying the 15-point
     kronrod rule (resk) obtained by optimal addition
     of abscissae to the 7-point gauss rule(resg).
  
     abserr - Estimate of the modulus of the absolute error,
     which should not exceed fabs(I-result)
     
     resabs - Approximation to the integral J
    
     resasc - Approximation to the integral of abs(f-I/(b-a)) over (a,b) 

     The abscissae and weights are given for the interval (-1,1).
     Because of symmetry only the positive abscissae and their
     corresponding weights are given. 

  /* Gauss quadrature weights and kronrod quadrature abscissae and
     weights as evaluated with 80 decimal digit arithmetic by
     L. W. Fullerton, Bell Labs, Nov. 1981. */

  const double wg[4] = /* weights of the 7-point gauss rule */
  {
    0.129484966168869693270611432679082,
    0.279705391489276667901467771423780,
    0.381830050505118944950369775488975,
    0.417959183673469387755102040816327
  };

  const double xgk[8] = /* weights of the 15-point kronrod rule */
  {
    0.991455371120812639206854697526329,
    0.949107912342758524526189684047851,
    0.864864423359769072789712788640926,
    0.741531185599394439863864773280788,
    0.586087235467691130294144838258730,
    0.405845151377397166906606412076961,
    0.207784955007898467600689403773245,
    0.000000000000000000000000000000000
  };

  const double wgk[8] =   /* abscissae of the 15-point kronrod rule */
  {
    0.022935322010529224963732008058970,   
    0.063092092629978553290700663189204, 
    0.104790010322250183839876322541518,
    0.140653259715525918745189590510238,
    0.169004726639267902826583426598550,
    0.190350578064785409913256402421014,
    0.204432940075298892414161999234649,
    0.209482141084727828012999174891714
  };

  /* xgk[1], xgk[3], ...  abscissae of the 7-point gauss rule. 
     
     xgk[0], xgk[2], ...  abscissae which are optimally added to the
     7-point gauss rule */
  
  double fv1[8], fv2[8];

  const double center = 0.5 * (a + b);
  const double half_length = 0.5 * (b - a);
  const double abs_half_length = fabs (half_length);
  const double fc = (*f) (center);

  /* result of 7-point gauss formula and of 15-point kronrod formula */

  double result_gauss7 = fc * wg[3];	
  double result_kronrod15 = fc * wgk[7];

  double tmp_resabs = fabs (result_kronrod15);
  double tmp_resasc = 0 ;

  int j ;

  for (j = 0; j < 3; j++)
    {
      const int jtw = j * 2 + 1;   /* j=1,2,3 jtw=2, 4, 6 */
      const double abscissa = half_length * xgk[jtw];
      const double fval1 = (*f) (center - abscissa);
      const double fval2 = (*f) (center + abscissa);
      const double fsum = fval1 + fval2;
      fv1[jtw] = fval1;
      fv2[jtw] = fval2;
      result_gauss7 +=  wg[j] * fsum;
      result_kronrod15 +=  wgk[jtw] * fsum;
      tmp_resabs += wgk[jtw] * (abs (fval1) + abs (fval2));
    }

  for (j = 0; j < 7; j += 2)
    {
      const double abscissa = half_length * xgk[j];
      const double fval1 = (*f)(center - abscissa);
      const double fval2 = (*f)(center + abscissa);
      fv1[j] = fval1;
      fv2[j] = fval2;
      result_kronrod15 += wgk[j] * (fval1 + fval2);
      tmp_resabs +=  wgk[j] * (fabs (fval1) + fabs (fval2));
    } ;

  {
    const double res_kronrod_mean = result_kronrod15 * 0.5;
    tmp_resasc = wgk[7] * fabs (fc - res_kronrod_mean);
    for (j = 0; j < 7; j++)
      {
	tmp_resasc += wgk[j] * (fabs (fv1[j] - res_kronrod_mean) 
				+ fabs (fv2[j] - res_kronrod_mean));
      }
  }

  *result = result_kronrod15 * half_length;

  tmp_resabs *= abs_half_length;
  tmp_resasc *= abs_half_length;
  
  *resabs = tmp_resabs ;
  *resasc = tmp_resasc ;

  {
    double tmp_abserr = fabs ((result_kronrod15 - result_gauss7) * half_length);
    if (tmp_resasc != 0 && tmp_abserr != 0)
      {
	double scale = pow((200 * tmp_abserr / tmp_resasc), 1.5) ;
	if (scale > 1)
	  {
	    *abserr = tmp_resasc * scale ;
	  }
	else 
	  {
	    *abserr = tmp_resasc ;
	  }
      }
    if (tmp_resabs > DBL_MIN / (50 * DBL_EPSILON))
      {
	double err = (DBL_EPSILON * 50) * tmp_resabs ;
	if (err > *abserr) 
	  {
	    *abserr = err ;
	  }
      }
  }
}
