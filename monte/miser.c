/* Based on the algorithm described in Numerical recipes, 
   but we make a few changes.  

   
   Seperate vectors of the upper and lower limits are used, 
   as in our version of vegas.  Also, arrays are zero based.

   */
/* Author: MJB */
/* RCS: $Id$ */

#include <config.h>
#include <math.h>
#include <stdlib.h>

/* gsl headers */
#include <gsl_math.h>
#include <gsl_vector.h>
#include <gsl_errno.h>
#include <gsl_rng.h>

#include <gsl_monte.h>
#include <gsl_miser.h>
#include <gsl_monte_plain.h>

#define PFAC 0.1
#define TINY 1.0e-30
#define BIG 1.0e30

#define myMAX(a,b) ((a) >= (b) ? (a) : (b))
#define myMIN(a,b) ((a) <= (b) ? (a) : (b))

#define COPYSIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

size_t min_calls = 15;
size_t min_calls_per_bisection = 60;
  
double dither;

static long iran;
static long idum;

void ranpt(const gsl_rng * r, gsl_vector* pt, double* xl, double* xh, size_t n);

inline void ranpt(const gsl_rng * r, gsl_vector* pt, double* xl, double* xh, size_t n)
{
  size_t j;

  for (j=0; j < n; j++)
    gsl_vector_set(pt, j, xl[j]+(xh[j]-xl[j])*gsl_rng_uniform(r));
}


int gsl_monte_miser(const gsl_rng * r, 
		    double (*func)(double []), double xl[], double xh[], 
		    size_t num_dim, size_t calls, double *avg, double *var)
{
  int status = 0;

  gsl_vector *xl_temp, *xh_temp;

  size_t n, npre, calls_l, calls_r;
  size_t j;
  int j_bisect;
  double avg_l, var_l;
  double frac_l, fval;
  double rgn_l, rgn_m, rgn_r, s, sig_l, sig_l_bisect, sig_r, sig_r_bisect;
  double sum, sum_bisect;

  gsl_vector *fmax_l, *fmax_r, *fmin_l, *fmin_r;
  gsl_vector *pt, *rmid;

  pt = gsl_vector_alloc(num_dim);

  if (calls < min_calls_per_bisection) {
    status = gsl_monte_plain(r, func, xl, xh, num_dim, calls, avg, var);
  }
  else {
    rmid = gsl_vector_alloc(num_dim);

    npre = myMAX((size_t)(calls*PFAC), min_calls);

    fmax_l = gsl_vector_alloc(num_dim);
    fmax_r = gsl_vector_alloc(num_dim);
    fmin_l = gsl_vector_alloc(num_dim);
    fmin_r = gsl_vector_alloc(num_dim);

    /* bisect the integration region, with some fuzz */
    for (j = 0; j < num_dim; j++) {
      iran = (iran*2661+36979) % 175000;
      s = COPYSIGN(dither, (double)(iran-87500));
      gsl_vector_set(rmid, j, (0.5+s)*xl[j]+(0.5-s)*xh[j]);
      fmin_l->data[j] = fmin_r->data[j] = BIG;
      fmax_l->data[j] = fmax_r->data[j] = -BIG;
    }
    
    /* Compute some function values to find the half-spaces with the 
       largest variance.  First we find the largest and smallest value
       in each half..
    */
    for (n = 1; n <= npre; n++) {
      ranpt(r, pt, xl, xh, num_dim);
      fval=(*func)(pt->data);
      for (j = 0; j < num_dim; j++) {
	if (pt->data[j] <= rmid->data[j]) {
	  fmin_l->data[j] = myMIN(fmin_l->data[j], fval);
	  fmax_l->data[j] = myMAX(fmax_l->data[j], fval);
	}
	else {
	  fmin_r->data[j] = myMIN(fmin_r->data[j], fval);
	  fmax_r->data[j] = myMAX(fmax_r->data[j], fval);
	}
      }
    }
    /* Now find direction with the largest variance */
    sum_bisect = BIG;
    j_bisect = -1;
    sig_l_bisect = sig_r_bisect = 1.0;
    for (j = 0; j < num_dim;j++) {
      if (fmax_l->data[j] > fmin_l->data[j] && 
	  fmax_r->data[j] > fmin_r->data[j]) {
	sig_l = myMAX(TINY, pow(fmax_l->data[j]-fmin_l->data[j], 2.0/3.0));
	sig_r = myMAX(TINY, pow(fmax_r->data[j]-fmin_r->data[j], 2.0/3.0));
	sum = sig_l+sig_r;
	if (sum <= sum_bisect) {
	  sum_bisect = sum;
	  j_bisect = j;
	  sig_l_bisect = sig_l;
	  sig_r_bisect = sig_r;
	}
      }
    }

    gsl_vector_free(fmin_r);
    gsl_vector_free(fmin_l);
    gsl_vector_free(fmax_r);
    gsl_vector_free(fmax_l);

    if (j_bisect < 0) 
      j_bisect = (num_dim*iran)/175000; /* All were same, so chose
					   direction at random */
    rgn_l = xl[j_bisect];
    rgn_m = rmid->data[j_bisect];
    rgn_r = xh[j_bisect];
    
    frac_l = fabs((rgn_m-rgn_l)/(rgn_r-rgn_l));
    calls_l = (unsigned long)
      (min_calls + (calls-npre-2*min_calls)*frac_l*sig_l_bisect
       /(frac_l*sig_l_bisect+(1.0-frac_l)*sig_r_bisect));
    calls_r = calls - npre - calls_l;

    xl_temp = gsl_vector_alloc(num_dim);
    xh_temp = gsl_vector_alloc(num_dim);

    for (j = 0; j < num_dim;j++) {
      xl_temp->data[j] = xl[j];
      xh_temp->data[j] = xh[j];
    }

    xh_temp->data[j_bisect] = rmid->data[j_bisect];
    status = gsl_monte_miser(r, func, xl_temp->data, xh_temp->data, 
		   num_dim, calls_l, &avg_l, &var_l);
    xl_temp->data[j_bisect] = rmid->data[j_bisect];
    xh_temp->data[j_bisect] = xh[j_bisect];
    status = gsl_monte_miser(r, func, xl_temp->data, xh_temp->data, 
		   num_dim, calls_r, avg, var);

    *avg = frac_l*avg_l+(1-frac_l)*(*avg);
    *var = frac_l*frac_l*var_l+(1-frac_l)*(1-frac_l)*(*var);

    gsl_vector_free(xl_temp);
    gsl_vector_free(xh_temp);
    gsl_vector_free(rmid);
  }

  gsl_vector_free(pt);

  return status;
}
