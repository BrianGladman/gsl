/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_gamma.h"
#include "gsl_sf_coupling.h"

#define locMax(a,b)  ((a)>(b) ? (a) : (b))
#define locMin(a,b)  ((a)<(b) ? (a) : (b))


/* See: [Thompson, Atlas for Computing Mathematical Functions] */


static
int
twice(double x)
{
  if(x >= 0.0)
    return (int) floor(2.0*x + 0.1);
  else 
    return (int)  ceil(2.0*x - 0.1);
}


static
int
delta(int ta, int tb, int tc, double * d)
{
  double f1, f2, f3, f4;
  int status = 0;
  status += gsl_sf_fact_impl((ta + tb - tc)/2, &f1);
  status += gsl_sf_fact_impl((ta + tc - tb)/2, &f2);
  status += gsl_sf_fact_impl((tb + tc - ta)/2, &f3);
  status += gsl_sf_fact_impl((tb + tc + ta)/2, &f4);
  if(status != 0) {
    *d = 0.0;
    return GSL_EOVRFLW;
  }
  *d = f1 * f2 * f3 / f4;
  return GSL_SUCCESS;
}


static
int
triangle_selection_fails(int two_ja, int two_jb, int two_jc)
{
  return ((two_jb < abs(two_ja - two_jc)) || (two_jb > two_ja + two_jc));
}


static
int
m_selection_fails(int two_ja, int two_jb, int two_jc,
                  int two_ma, int two_mb, int two_mc)
{
  return (   abs(two_ma) > two_ja 
          || abs(two_mb) > two_jb
	  || abs(two_mc) > two_jc
	  || GSL_IS_ODD(two_ja + two_ma)
	  || GSL_IS_ODD(two_jb + two_mb)
	  || GSL_IS_ODD(two_jc + two_mc)
	  );
}

int
gsl_sf_coupling_3j_impl(double ja, double jb, double jc,
                        double ma, double mb, double mc,
			double * result)
{
  int two_ja = twice(ja);
  int two_jb = twice(jb);
  int two_jc = twice(jc);

  int two_ma = twice(ma);
  int two_mb = twice(mb);
  int two_mc = twice(mc);

  if(   triangle_selection_fails(two_ja, two_jb, two_jc)
     || m_selection_fails(two_ja, two_jb, two_jc, two_ma, two_mb, two_mc)
     ) {
    *result = 0.0;
    return GSL_SUCCESS;
  }
  else {
    double n1_a, n1_b, n3_a, n3_b;
    double d1_a, d1_b, d2_a, d2_b, d3_a, d3_b;
    double n1, n2, n3;
    double d1, d2, d3;
    double norm;
    double sign = (GSL_IS_ODD((two_ja - two_jb - two_jc)/2) ? -1.0 : 1.0);
    int tk, tkmin, tkmax;
    double sum = 0.0;
    double phase;
    int status = 0;
    status += gsl_sf_fact_impl((two_jc + two_ja - two_jb)/2, &n1_a);
    status += gsl_sf_fact_impl((two_jc - two_ja + two_jb)/2, &n1_b);
    status += gsl_sf_fact_impl((two_ja + two_jb - two_jc)/2, &n2);
    status += gsl_sf_fact_impl((two_jc - two_mc)/2, &n3_a);
    status += gsl_sf_fact_impl((two_jc + two_mc)/2, &n3_b);
    status += gsl_sf_fact_impl((two_ja + two_jb + two_jc)/2, &d1);
    status += gsl_sf_fact_impl((two_ja - two_ma)/2, &d2_a);
    status += gsl_sf_fact_impl((two_ja + two_ma)/2, &d2_b);
    status += gsl_sf_fact_impl((two_jb - two_mb)/2, &d3_a);
    status += gsl_sf_fact_impl((two_jb + two_mb)/2, &d3_b);
    if(status != GSL_SUCCESS) {
      *result = 0.0;
      return GSL_EOVRFLW;
    }
    n1 = n1_a * n1_b;
    n3 = n3_a * n3_b;
    d2 = d2_a * d2_b;
    d3 = d3_a * d3_b;
    
    norm = sign * sqrt(n1*n2*n3)/sqrt(d1*d2*d3);
    
    tkmin = locMax(0, two_jb - two_ja - two_jc);
    tkmax = locMin(two_jc - two_ja + two_jb, two_jc - two_mc);
    
    phase = GSL_IS_ODD((tkmin + two_jb + two_mb)/2) ? -1.0 : 1.0;

    for(tk=tkmin; tk<=tkmax; tk += 2) {
      double term;

      status = 0;
      status += gsl_sf_fact_impl((two_ja + two_jb + two_jc - tk)/2, &n1);
      status += gsl_sf_fact_impl((two_ja - two_ma + tk)/2, &n2);
      status += gsl_sf_fact_impl(tk/2, &d1_a);
      status += gsl_sf_fact_impl((two_jc - two_ja + two_jb - tk)/2, &d1_b);
      status += gsl_sf_fact_impl((two_jc - two_mc - tk)/2, &d2);
      status += gsl_sf_fact_impl((two_ja - two_jb + two_mc + tk)/2, &d3);

      if(status != 0) {
        *result = 0.0;
	return GSL_EOVRFLW;
      }

      term = phase * n1 * n2 / (d1 * d2 * d3);
      phase = -phase;
      sum += norm * term;
    }
    
    *result = sum;
    return GSL_SUCCESS;
  }
}


int
gsl_sf_coupling_6j_impl(double ja, double jb, double jc,
                        double jd, double je, double jf,
			double * result)
{
  int two_ja = twice(ja);
  int two_jb = twice(jb);
  int two_jc = twice(jc);
  int two_jd = twice(jd);
  int two_je = twice(je);
  int two_jf = twice(jf);
  
  if(   triangle_selection_fails(two_ja, two_jb, two_je)
     || triangle_selection_fails(two_ja, two_jc, two_jf)
     || triangle_selection_fails(two_jb, two_jd, two_jf)
     || triangle_selection_fails(two_jc, two_jd, two_je)
     ) {
    *result = 0.0;
    return GSL_SUCCESS;
  }
  else {
    double n1;
    double d1, d2, d3, d4, d5, d6;
    double norm;
    int tk, tkmin, tkmax;
    double phase;
    double sum = 0.0;
    int status = 0;
    status += delta(two_ja, two_jb, two_je, &d1);
    status += delta(two_ja, two_jc, two_jf, &d2);
    status += delta(two_jb, two_jd, two_jf, &d3);
    status += delta(two_jc, two_jd, two_je, &d4);
    if(status != GSL_SUCCESS) {
      *result = 0.0;
      return GSL_EOVRFLW;
    }
    norm = sqrt(d1) * sqrt(d2) * sqrt(d3) * sqrt(d4);
    
    tkmin = locMax3(0,
                   two_ja + two_jd - two_je - two_jf,
                   two_jb + two_jc - two_je - two_jf);

    tkmax = locMin4(two_ja + two_jb + two_jc + two_jd + 2,
                    two_ja + two_jb - two_je,
		    two_jc + two_jd - two_je,
		    two_ja + two_jc - two_jf,
		    two_jb + two_jd - two_jf);

    phase = GSL_IS_ODD((two_ja + two_jb + two_jc + two_jd + tkmin)/2)
            ? -1.0
	    :  1.0;

    for(tk=tkmin; tk<=tkmax; tk += 2) {
      double term;
      double d1_a, d1_b;
      status = 0;
      
      status += gsl_sf_fact_impl((two_ja + two_jb + two_jc + two_jd - tk)/2 + 1, &n1);
      status += gsl_sf_fact_impl(tk/2, &d1_a);
      status += gsl_sf_fact_impl((two_je + two_jf - two_ja - two_jd + tk)/2, &d1_b);
      status += gsl_sf_fact_impl((two_je + two_jf - two_jb - two_jc + tk)/2, &d2);
      status += gsl_sf_fact_impl((two_ja + two_jb - two_je - tk)/2, &d3);
      status += gsl_sf_fact_impl((two_jc + two_jd - two_je - tk)/2, &d4);
      status += gsl_sf_fact_impl((two_ja + two_jc - two_jf - tk)/2, &d5);
      status += gsl_sf_fact_impl((two_jb + two_jd - two_jf - tk)/2, &d6);
      
      if(status != GSL_SUCCESS) {
        *result = 0.0;
	return GSL_EOVRFLW;
      }

      d1 = d1_a * d1_b;
      
      term  = phase * n1 / (d1*d2*d3) / (d4*d5*d6);
      phase = -phase;
      sum  += norm * term;
    }
    
    *result = sum;
    return GSL_SUCCESS;
  }
}


int gsl_sf_coupling_3j_e(double ja, double jb, double jc,
                         double ma, double mb, double mc,
		         double * result)
{
  int status = gsl_sf_coupling_3j_impl(ja, jb, jc, ma, mb, mc, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_coupling_3j_e", status);
  }
  return status;
}


int gsl_sf_coupling_6j_e(double ja, double jb, double jc,
                         double jd, double je, double jf,
			 double * result)
{
  int status = gsl_sf_coupling_6j_impl(ja, jb, jc, jd, je, jf, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_coupling_6j_e", status);
  }
  return status;
}



double gsl_sf_coupling_3j(double ja, double jb, double jc,
                          double ma, double mb, double mc)
{
  double y;
  int status = gsl_sf_coupling_3j_impl(ja, jb, jc, ma, mb, mc, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_coupling_3j", status);
  }
  return y;
}


double gsl_sf_coupling_6j(double ja, double jb, double jc,
                          double jd, double je, double jf)
{
  double y;
  int status = gsl_sf_coupling_6j_impl(ja, jb, jc, jd, je, jf, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_coupling_6j", status);
  }
  return y;
}
