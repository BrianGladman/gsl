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

#define SQR(a) ((a)*(a))
#define myMAX(a,b) ((a) >= (b) ? (a) : (b))
#define myMIN(a,b) ((a) <= (b) ? (a) : (b))

/* these should be in the state structure */
size_t min_calls = 15;
size_t min_calls_per_bisection = 60;
double precond_frac =  0.1;
double ALPHA = 2.0/3.0; 
double dither;


int gsl_monte_miser(const gsl_rng * r, 
		    double (*func)(double []), double xl[], double xu[], 
		    size_t num_dim, size_t calls, double *res, double *err)
{
  int status = 0;

  gsl_vector *xl_tmp, *xu_tmp;

  size_t n, pre_calls, calls_l, calls_r;
  size_t i;
  int i_bisect;
  double res_l, err_l;
  double fraction_l;
  double f;
  double xbi_l, xbi_m, xbi_r, s;
  double sigma_l, sigma_r;
  double weight_l, weight_r;
  double sum, sum_bisect;

  gsl_vector *fmax_l, *fmax_r, *fmin_l, *fmin_r;
  gsl_vector *x, *x_mid;

  x = gsl_vector_alloc(num_dim);

  if (calls < min_calls_per_bisection) {
    status = gsl_monte_plain(r, func, xl, xu, num_dim, calls, res, err);
  }
  else {
    pre_calls = myMAX((size_t)(calls*precond_frac), min_calls);

    x_mid = gsl_vector_alloc(num_dim);
    fmax_l = gsl_vector_alloc(num_dim);
    fmax_r = gsl_vector_alloc(num_dim);
    fmin_l = gsl_vector_alloc(num_dim);
    fmin_r = gsl_vector_alloc(num_dim);

    for (i = 0; i < num_dim; i++) {

      /* flip a coin to bisect the integration region with some fuzz */
      s = (0.5 - gsl_rng_uniform(r) >= 0.0 ? dither : -dither ); 
      x_mid->data[i] = (0.5 + s)*xl[i] + (0.5 - s)*xu[i];

      fmin_l->data[i] = fmin_r->data[i] = DBL_MAX;
      fmax_l->data[i] = fmax_r->data[i] = -DBL_MAX;
    }
    
    /* 
       The idea is to chose the direction to bisect based on which will
       give the smallest total variance.  We could (and may do so later)
       use MC to compute these variances.  But the NR guys simply estimate
       the variances by finding the min and max function values 
       for each half-region for each bisection.
    */
    for (n = 1; n <= pre_calls; n++) {
      for (i = 0; i < num_dim; i++)
	x->data[i] = xl[i] + (xu[i] - xl[i])*gsl_rng_uniform(r);
      f = (*func)(x->data);
      for (i = 0; i < num_dim; i++) {
	if (x->data[i] <= x_mid->data[i]) {
	  fmin_l->data[i] = myMIN(fmin_l->data[i], f);
	  fmax_l->data[i] = myMAX(fmax_l->data[i], f);
	}
	else {
	  fmin_r->data[i] = myMIN(fmin_r->data[i], f);
	  fmax_r->data[i] = myMAX(fmax_r->data[i], f);
	}
      }
    }
    /* Now find direction with the largest variance */
    sum = 0;
    sum_bisect = DBL_MAX;
    i_bisect = -1;
    weight_l = weight_r = 1.0;

    for (i = 0; i < num_dim;i++) {
      if (fmax_l->data[i] > fmin_l->data[i] && 
	  fmax_r->data[i] > fmin_r->data[i]) {
	sigma_l = myMAX(GSL_MACH_EPS, fmax_l->data[i] - fmin_l->data[i]);
	sigma_r = myMAX(GSL_MACH_EPS, fmax_r->data[i] - fmin_r->data[i]);
	sum = pow(sigma_l, ALPHA) + pow(sigma_r, ALPHA);

	if (sum <= sum_bisect) {
	  sum_bisect = sum;
	  i_bisect = i;
	  weight_l = pow(sigma_l, ALPHA);
	  weight_r = pow(sigma_r, ALPHA);
	}
      }
    }

    if (i_bisect < 0) {
      /* All were same, so chose direction at random */
      
      i_bisect = num_dim*gsl_rng_uniform(r); 
    }

    xbi_l = xl[i_bisect];
    xbi_m = x_mid->data[i_bisect];
    xbi_r = xu[i_bisect];

    /* get the actual fractional sizes of the two "halves" */
    fraction_l = fabs((xbi_m - xbi_l)/(xbi_r - xbi_l));
    calls_l = min_calls;
    calls_l += (size_t)
      (calls - pre_calls - 2*min_calls)*fraction_l*weight_l
       /(fraction_l*weight_l + (1.0 - fraction_l)*weight_r);
    calls_r = calls - pre_calls - calls_l;

    gsl_vector_free(fmin_r);
    gsl_vector_free(fmin_l);
    xl_tmp = fmax_r; /* reuse vectors */
    xu_tmp = fmax_l;
    for (i = 0; i < num_dim; i++) {
      xl_tmp->data[i] = xl[i];
      xu_tmp->data[i] = xu[i];
    }

    xu_tmp->data[i_bisect] = x_mid->data[i_bisect];
    status = gsl_monte_miser(r, func, xl_tmp->data, xu_tmp->data, 
		   num_dim, calls_l, &res_l, &err_l);
    xl_tmp->data[i_bisect] = x_mid->data[i_bisect];
    xu_tmp->data[i_bisect] = xu[i_bisect];
    status = gsl_monte_miser(r, func, xl_tmp->data, xu_tmp->data, 
		   num_dim, calls_r, res, err);

    *res += res_l;
    *err = sqrt( SQR(err_l) + SQR(*err) );   

    gsl_vector_free(xl_tmp);
    gsl_vector_free(xu_tmp);
    gsl_vector_free(x_mid);
  }

  gsl_vector_free(x);

  return status;
}
