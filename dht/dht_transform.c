/* dht/dht_transform.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include "gsl_dht.h"



/* Handle internal allocation for a given length. */
static int
dht_transform_allocator(gsl_dht_transform * t, size_t size)
{
  t->j = (double *)malloc((size+2)*sizeof(double));
  if(t->j == 0) {
    return GSL_ENOMEM;
  }

  t->Jjj = (double *)malloc(size*(size+1)/2 * sizeof(double));
  if(t->Jjj == 0) {
    free(t->j);
    return GSL_ENOMEM;
  }

  t->J2 = (double *)malloc((size+1)*sizeof(double));
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
  for(s=1; s < t->size + 2; s++) {
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
    GSL_ERROR_VAL("gsl_dht_transform_new: size == 0", GSL_EDOM, 0);
  }
  else if(nu < 0.0) {
    GSL_ERROR_VAL("gsl_dht_transform_new: nu < 0.0", GSL_EDOM, 0);
  }
  else if(xmax <= 0.0) {
    GSL_ERROR_VAL("gsl_dht_transform_new: xmax <= 0.0", GSL_EDOM, 0);
  }
  else {
    gsl_dht_transform * t = (gsl_dht_transform *)malloc(sizeof(gsl_dht_transform));
    if(t == 0) {
      GSL_ERROR_VAL("gsl_dht_transform_new: out of memory", GSL_ENOMEM, 0);
    }
    else {
      t->nu   = -1.0; /* Make it clear that this needs to be calculated. */
      t->size = size;
      if(dht_transform_allocator(t, size) != GSL_SUCCESS) {
        free(t);
	GSL_ERROR_VAL("gsl_dht_transform_new: alloc failed", GSL_ENOMEM, 0);
      }
      else {
        if(gsl_dht_transform_recalc_impl(t, nu, xmax) != GSL_SUCCESS) {
	  free(t->J2);
	  free(t->Jjj);
	  free(t->j);
          free(t);
          GSL_ERROR_VAL("gsl_dht_transform_new: recalc failed", GSL_EFAILED, 0);
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
    double jN;

    if(nu != t->nu) {
      /* Recalculate Bessel zeros if necessary. */
      t->nu = nu;
      stat_bz = dht_bessel_zeros(t);
    }

    jN = t->j[t->size+1];

    t->xmax = xmax;
    t->kmax = jN / xmax;

    t->J2[0] = 0.0;
    for(m=1; m<t->size+1; m++) {
      gsl_sf_result J;
      stat_J += gsl_sf_bessel_Jnu_impl(nu + 1.0, t->j[m], &J);
      t->J2[m] = J.val * J.val;
    }

    /* J_nu(j[n] j[m] / j[N]) = Jjj[n(n-1)/2 + m - 1], 1 <= n,m <= size
     */
    for(n=1; n<t->size+1; n++) {
      for(m=1; m<=n; m++) {
        double arg = t->j[n] * t->j[m] / jN;
        gsl_sf_result J;
        stat_J += gsl_sf_bessel_Jnu_impl(nu, arg, &J);
        t->Jjj[n*(n-1)/2 + m - 1] = J.val;
      }
    }

    if(stat_J != 0) {
      return GSL_EFAILED;
    }
    else {
      return stat_bz;
    }
  }
}


double gsl_dht_transform_x_sample(const gsl_dht_transform * t, int n)
{
  return t->j[n+1]/t->j[t->size+1] * t->xmax;
}


double gsl_dht_transform_k_sample(const gsl_dht_transform * t, int n)
{
  return t->j[n+1] / t->xmax;
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
  const double jN = t->j[t->size + 1];
  const double r  = t->xmax / jN;
  size_t m;
  size_t i;

  for(m=0; m<t->size; m++) {
    double sum = 0.0;
    double Y;
    for(i=0; i<t->size; i++) {
      /* Need to find max and min so that we
       * address the symmetric Jjj matrix properly.
       * FIXME: we can presumably optimize this
       * by just running over the elements of Jjj
       * in a deterministic manner.
       */
      size_t m_local; 
      size_t n_local;
      if(i < m) {
        m_local = i;
	n_local = m;
      }
      else {
        m_local = m;
	n_local = i;
      }
      Y = t->Jjj[n_local*(n_local+1)/2 + m_local] / t->J2[i+1];
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
