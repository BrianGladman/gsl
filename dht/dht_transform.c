/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_math.h>
#include <gsl_sf_bessel.h>
#include "gsl_dht.h"



/* Handle internal allocation for a given length. */
static int
dht_transform_allocator(gsl_dht_transform * t, size_t N)
{
  t->j = (double *)malloc((N+1)*sizeof(double));
  if(t->j == 0) {
    return GSL_ENOMEM;
  }

  t->Jjj = (double *)malloc(N*(N-1)/2 * sizeof(double));
  if(t->Jjj == 0) {
    free(t->j);
    return GSL_ENOMEM;
  }

  t->J2 = (double *)malloc((N+1)*sizeof(double));
  if(t->J2 == 0) {
    free(t->Jjj);
    free(t->j);
    return GSL_ENOMEM;
  }

  return GSL_SUCCESS;
}


/* Handle internal calculation of Bessel zeros. */
static int
dht_bessel_zeros(gsl_dht_transform * t)
{
  size_t s;
  gsl_sf_result z;
  int stat_z = 0;
  t->j[0] = 0.0;
  for(s=1; s <= t->size; s++) {
    stat_z += gsl_sf_bessel_zero_Jnu_impl(t->nu, s, &z);
    t->j[s] = z.val;
  }
  if(stat_z != 0) {
    return GSL_EFAILED;
  }
  else {
    return GSL_SUCCESS;
  }
}



gsl_dht_transform *
gsl_dht_transform_new(size_t size, double nu, double xmax)
{
  if(size == 0) {
    GSL_ERROR_RETURN("gsl_dht_transform_new: size == 0", GSL_EDOM, 0);
  }
  else if(nu < 0.0) {
    GSL_ERROR_RETURN("gsl_dht_transform_new: nu < 0.0", GSL_EDOM, 0);
  }
  else if(xmax <= 0.0) {
    GSL_ERROR_RETURN("gsl_dht_transform_new: xmax <= 0.0", GSL_EDOM, 0);
  }
  else {
    gsl_dht_transform * t = (gsl_dht_transform *)malloc(sizeof(gsl_dht_transform));
    if(t == 0) {
      GSL_ERROR_RETURN("gsl_dht_transform_new: out of memory", GSL_ENOMEM, 0);
    }
    else {
      const size_t N = size + 1;
      int stat_al = dht_transform_allocator(t, N);
      if(stat_al != GSL_SUCCESS) {
        free(t);
	GSL_ERROR_RETURN("gsl_dht_transform_new: alloc failed", GSL_ENOMEM, 0);
      }
      else {
        int stat_re = gsl_dht_transform_recalc_impl(t, nu, xmax);
        if(stat_re != GSL_SUCCESS) {
          free(t);
          GSL_ERROR_RETURN("gsl_dht_transform_new: recalc failed", GSL_EFAILED, 0);
        }
      }
      return t;
    }
  }
}


int
gsl_dht_transform_recalc_impl(gsl_dht_transform * t, double nu, double xmax)
{
  if(xmax <= 0.0 || nu < 0.0) {
    return GSL_EDOM;
  }
  else {
    size_t n, m;
    int stat_bz = GSL_SUCCESS;
    int stat_J  = 0;

    if(nu != t->nu) {
      /* Recalculate Bessel zeros if necessary. */
      stat_bz = dht_bessel_zeros(t);
    }

    t->nu   = nu;
    t->xmax = xmax;
    t->kmax = t->j[t->size+1] / t->xmax;

    t->J2[0] = 0.0;
    for(n=1; n<t->size+1; n++) {
      gsl_sf_result J;
      stat_J += gsl_sf_bessel_Jnu_impl(t->nu + 1.0, t->j[n], &J);
      t->J2[n] = J.val * J.val;
    }

    /* J_nu(j[n] j[m] / j[N]) = Jjj[n(n-1)/2 + m - 1]
     */
    for(n=1; n<t->size+1; n++) {
      for(m=1; m<=n; m++) {
        double arg = t->j[n] * t->j[m] / t->j[t->size + 1];
        gsl_sf_result J;
        stat_J += gsl_sf_bessel_Jnu_impl(t->nu, arg, &J);
        t->Jjj[n*(n-1)/2 + m - 1] = J.val;
      }
    }

    if(stat_J != 0) {
      return stat_bz;
    }
    else {
      return GSL_EFAILED;
    }
  }
}


double gsl_dht_transform_x_sample(const gsl_dht_transform * t, int n)
{
  return t->j[n-1]/t->j[t->size+1] * t->xmax;
}


double gsl_dht_transform_k_sample(const gsl_dht_transform * t, int n)
{
  return t->j[n-1] / t->xmax;
}


void gsl_dht_transform_free(gsl_dht_transform * t)
{
  if(t != 0) {
    free(t->J2);
    free(t->Jjj);
    free(t->j);
    free(t);
  }
}


int
gsl_dht_transform_apply(const gsl_dht_transform * t, double * f_in, double * f_out)
{
  double r = t->xmax / t->j[t->size + 1];
  size_t m;
  size_t i;

  for(m=0; m<t->size; m++) {
    double sum = 0.0;
    double Y;
    for(i=0; i<t->size; i++) {
      size_t m_local; 
      size_t n_local;
      if(i < m) {
        m_local = i + 1;
	n_local = m + 1;
      }
      else {
        m_local = m + 1;
	n_local = i + 1;
      }
      Y = t->Jjj[n_local*(n_local-1)/2 + m_local - 1] / t->J2[m_local];
      sum += Y * f_in[i];
    }
    f_out[m] = sum * 2.0 * r*r;
  }

  return GSL_SUCCESS;
}


int
gsl_dht_transform_recalc_e(gsl_dht_transform * t, double nu, double xmax)
{
  int status = gsl_dht_transform_recalc_impl(t, nu, xmax);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_dht_transform_recalc_e", status);
  }
  return status;
}
