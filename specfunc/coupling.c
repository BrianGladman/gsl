/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_gamma.h"
#include "gsl_sf_coupling.h"


#ifdef HAVE_INLINE
inline
#endif
static
int locMax3(const int a, const int b, const int c)
{
  int d = GSL_MAX(a, b);
  return GSL_MAX(d, c);
}

#ifdef HAVE_INLINE
inline
#endif
static
int locMin3(const int a, const int b, const int c)
{
  int d = GSL_MIN(a, b);
  return GSL_MIN(d, c);
}

#ifdef HAVE_INLINE
inline
#endif
static
int locMin5(const int a, const int b, const int c, const int d, const int e)
{
  int f = GSL_MIN(a, b);
  int g = GSL_MIN(c, d);
  int h = GSL_MIN(f, g);
  return GSL_MIN(e, h);
}


/* See: [Thompson, Atlas for Computing Mathematical Functions] */

static
int
delta(int ta, int tb, int tc, double * d)
{
  double f1, f2, f3, f4;
  int status = 0;
  status += gsl_sf_fact_impl((ta + tb - tc)/2, &f1);
  status += gsl_sf_fact_impl((ta + tc - tb)/2, &f2);
  status += gsl_sf_fact_impl((tb + tc - ta)/2, &f3);
  status += gsl_sf_fact_impl((ta + tb + tc)/2 + 1, &f4);
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


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_coupling_3j_impl(int two_ja, int two_jb, int two_jc,
                        int two_ma, int two_mb, int two_mc,
			double * result)
{
  if(two_ja < 0 || two_jb < 0 || two_jc < 0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(   triangle_selection_fails(two_ja, two_jb, two_jc)
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
    double sign = (GSL_IS_ODD((two_ja - two_jb - two_mc)/2) ? -1.0 : 1.0);
    int tk, tkmin, tkmax;
    double sum = 0.0;
    double phase;
    int status = 0;
    status += gsl_sf_fact_impl((two_jc + two_ja - two_jb)/2, &n1_a);
    status += gsl_sf_fact_impl((two_jc - two_ja + two_jb)/2, &n1_b);
    status += gsl_sf_fact_impl((two_ja + two_jb - two_jc)/2, &n2);
    status += gsl_sf_fact_impl((two_jc - two_mc)/2, &n3_a);
    status += gsl_sf_fact_impl((two_jc + two_mc)/2, &n3_b);
    status += gsl_sf_fact_impl((two_ja + two_jb + two_jc)/2 + 1, &d1);
    status += gsl_sf_fact_impl((two_ja - two_ma)/2, &d2_a);
    status += gsl_sf_fact_impl((two_ja + two_ma)/2, &d2_b);
    status += gsl_sf_fact_impl((two_jb - two_mb)/2, &d3_a);
    status += gsl_sf_fact_impl((two_jb + two_mb)/2, &d3_b);

    if(status != 0) {
      *result = 0.0;
      return GSL_EOVRFLW;
    }

    n1 = n1_a * n1_b;
    n3 = n3_a * n3_b;
    d2 = d2_a * d2_b;
    d3 = d3_a * d3_b;

    norm = sign * sqrt(n1*n2*n3)/sqrt(d1*d2*d3);

    tkmin = GSL_MAX(0, two_jb - two_ja - two_mc);
    tkmax = GSL_MIN(two_jc - two_ja + two_jb, two_jc - two_mc);
    
    phase = GSL_IS_ODD((tkmin + two_jb + two_mb)/2) ? -1.0 : 1.0;

    for(tk=tkmin; tk<=tkmax; tk += 2) {
      double term;

      status = 0;
      status += gsl_sf_fact_impl((two_jb + two_jc + two_ma - tk)/2, &n1);
      status += gsl_sf_fact_impl((two_ja - two_ma + tk)/2, &n2);
      status += gsl_sf_fact_impl(tk/2, &d1_a);
      status += gsl_sf_fact_impl((two_jc - two_ja + two_jb - tk)/2, &d1_b);
      status += gsl_sf_fact_impl((two_jc - two_mc - tk)/2, &d2);
      status += gsl_sf_fact_impl((two_ja - two_jb + two_mc + tk)/2, &d3);

      if(status != 0) {
        *result = 0.0;
	return GSL_EOVRFLW;
      }

      d1 = d1_a * d1_b;

      term = phase * n1 * n2 / (d1 * d2 * d3);
      phase = -phase;
      sum += norm * term;
    }
    
    *result = sum;
    return GSL_SUCCESS;
  }
}


int
gsl_sf_coupling_6j_impl(int two_ja, int two_jb, int two_jc,
                        int two_jd, int two_je, int two_jf,
			double * result)
{
  if(   two_ja < 0 || two_jb < 0 || two_jc < 0
     || two_jd < 0 || two_je < 0 || two_je < 0
     ) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(   triangle_selection_fails(two_ja, two_jb, two_je)
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

    tkmax = locMin5(two_ja + two_jb + two_jc + two_jd + 2,
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


int
gsl_sf_coupling_9j_impl(int two_ja, int two_jb, int two_jc,
                        int two_jd, int two_je, int two_jf,
			int two_jg, int two_jh, int two_ji,
			double * result)
{
  if(   two_ja < 0 || two_jb < 0 || two_jc < 0
     || two_jd < 0 || two_je < 0 || two_jf < 0
     || two_jg < 0 || two_jh < 0 || two_ji < 0
     ) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(   triangle_selection_fails(two_ja, two_jb, two_jc)
          || triangle_selection_fails(two_jd, two_je, two_jf)
          || triangle_selection_fails(two_jg, two_jh, two_ji)
          || triangle_selection_fails(two_ja, two_jd, two_jg)
          || triangle_selection_fails(two_jb, two_je, two_jh)
          || triangle_selection_fails(two_jc, two_jf, two_ji)
     ) {
    *result = 0.0;
    return GSL_SUCCESS;
  }
  else {
    int tk;
    int tkmin = locMax3(abs(two_ja-two_ji), abs(two_jh-two_jd), abs(two_jb-two_jf));
    int tkmax = locMin3(two_ja + two_ji, two_jh + two_jd, two_jb + two_jf);
    double sum = 0.0;
    double phase;
    for(tk=tkmin; tk<=tkmax; tk += 2) {
      double s1, s2, s3;
      double term;
      int status = 0;
      status += gsl_sf_coupling_6j_impl(two_ja, two_ji, two_jd,  two_jh, tk, two_jg,  &s1);
      status += gsl_sf_coupling_6j_impl(two_jb, two_jf, two_jh,  two_jd, tk, two_je,  &s2);
      status += gsl_sf_coupling_6j_impl(two_ja, two_ji, two_jb,  two_jf, tk, two_jc,  &s3);
      if(status != GSL_SUCCESS) {
        *result = 0.0;
	return GSL_EOVRFLW;
      }
      term = s1 * s2 * s3;
      sum += (tk + 1) * term;
    }
    
    phase = GSL_IS_ODD(tkmin) ? -1.0 : 1.0;
    *result = phase * sum;
    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_coupling_3j_e(int two_ja, int two_jb, int two_jc,
                         int two_ma, int two_mb, int two_mc,
		         double * result)
{
  int status = gsl_sf_coupling_3j_impl(two_ja, two_jb, two_jc,
                                       two_ma, two_mb, two_mc,
				       result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_coupling_3j_e", status);
  }
  return status;
}


int gsl_sf_coupling_6j_e(int two_ja, int two_jb, int two_jc,
                         int two_jd, int two_je, int two_jf,
			 double * result)
{
  int status = gsl_sf_coupling_6j_impl(two_ja, two_jb, two_jc,
                                       two_jd, two_je, two_jf,
				       result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_coupling_6j_e", status);
  }
  return status;
}


int gsl_sf_coupling_9j_e(int two_ja, int two_jb, int two_jc,
                         int two_jd, int two_je, int two_jf,
			 int two_jg, int two_jh, int two_ji,
			 double * result)
{
  int status = gsl_sf_coupling_9j_impl(two_ja, two_jb, two_jc,
                                       two_jd, two_je, two_jf,
				       two_jg, two_jh, two_ji,
				       result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_coupling_9j_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*/

double gsl_sf_coupling_3j(int two_ja, int two_jb, int two_jc,
                          int two_ma, int two_mb, int two_mc)
{
  double y;
  int status = gsl_sf_coupling_3j_impl(two_ja, two_jb, two_jc,
                                       two_ma, two_mb, two_mc,
				       &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_coupling_3j", status);
  }
  return y;
}


double gsl_sf_coupling_6j(int two_ja, int two_jb, int two_jc,
                          int two_jd, int two_je, int two_jf)
{
  double y;
  int status = gsl_sf_coupling_6j_impl(two_ja, two_jb, two_jc,
                                       two_jd, two_je, two_jf,
				       &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_coupling_6j", status);
  }
  return y;
}

double gsl_sf_coupling_9j(int two_ja, int two_jb, int two_jc,
                          int two_jd, int two_je, int two_jf,
			  int two_jg, int two_jh, int two_ji)
{
  double y;
  int status = gsl_sf_coupling_9j_impl(two_ja, two_jb, two_jc,
                                       two_jd, two_je, two_jf,
				       two_jg, two_jh, two_ji,
				       &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_coupling_9j", status);
  }
  return y;
}
