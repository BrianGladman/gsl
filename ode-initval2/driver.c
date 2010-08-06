/* ode-initval/driver.c
 * 
 * Copyright (C) 2009 Tuomo Keskitalo
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Driver routine for gsl_odeiv2. This is a wrapper for low level GSL
   functions that allows a simple interface to step, control and
   evolve layers.
 */

#include <config.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_machine.h>

static gsl_odeiv2_driver *
driver_alloc (gsl_odeiv2_system * sys, const double hstart,
              const gsl_odeiv2_step_type * T)
{
  /* Allocates and initializes an ODE driver system. Step and evolve
     objects are allocated here, but control object is allocated in
     another function.
   */

  gsl_odeiv2_driver *state =
    (gsl_odeiv2_driver *) malloc (sizeof (gsl_odeiv2_driver));

  if (state == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for driver state",
                      GSL_ENOMEM);
    }

  if (sys == NULL)
    {
      GSL_ERROR_NULL ("gsl_odeiv2_system must be defined", GSL_EINVAL);
    }

  const size_t dim = sys->dimension;

  if (dim == 0)
    {
      GSL_ERROR_NULL ("gsl_odeiv2_system dimension must be a positive integer",
                      GSL_EINVAL);
    }

  state->sys = sys;

  state->s = gsl_odeiv2_step_alloc (T, dim);

  if (state->s == NULL)
    {
      free (state);
      GSL_ERROR_NULL ("failed to allocate step object", GSL_ENOMEM);
    }

  state->e = gsl_odeiv2_evolve_alloc (dim);

  if (state->e == NULL)
    {
      gsl_odeiv2_step_free (state->s);
      free (state);
      GSL_ERROR_NULL ("failed to allocate evolve object", GSL_ENOMEM);
    }

  if (hstart < 0.0)
    {
      state->h = hstart;
      state->hmax = -1.0 * GSL_DBL_MAX;
    }
  else if (hstart > 0.0)
    {
      state->h = hstart;
      state->hmax = GSL_DBL_MAX;
    }
  else
    {
      GSL_ERROR_NULL ("invalid hstart", GSL_EINVAL);
    }

  state->hmin = 0.0;
  state->nmax = 0;
  state->n = 0;
  state->c = NULL;

  return state;
}

int
gsl_odeiv2_driver_set_hmin (gsl_odeiv2_driver * D, const double hmin)
{
  /* Sets minimum allowed step size (hmin) for driver */

  if ((D->h < 0.0 && hmin > 0.0) || (D->h > 0.0 && hmin < 0.0))
    {
      GSL_ERROR_NULL ("min step direction must match step direction",
                      GSL_EINVAL);
    }

  if (hmin >= 0.0 || hmin < 0.0)
    {
      D->hmin = hmin;
    }
  else
    {
      GSL_ERROR_NULL ("invalid hmin", GSL_EINVAL);
    }

  return GSL_SUCCESS;
}

int
gsl_odeiv2_driver_set_hmax (gsl_odeiv2_driver * D, const double hmax)
{
  /* Sets maximum allowed step size (hmax) for driver */

  if ((D->h < 0.0 && hmax > 0.0) || (D->h > 0.0 && hmax < 0.0))
    {
      GSL_ERROR_NULL ("max step direction must match step direction",
                      GSL_EINVAL);
    }

  if (hmax > 0.0 || hmax < 0.0)
    {
      D->hmax = hmax;
    }
  else
    {
      GSL_ERROR_NULL ("invalid hmax", GSL_EINVAL);
    }

  return GSL_SUCCESS;
}

int
gsl_odeiv2_driver_set_nmax (gsl_odeiv2_driver * D, const size_t nmax)
{
  /* Sets maximum number of allowed steps (nmax) for driver */

  D->nmax = nmax;

  return GSL_SUCCESS;
}


gsl_odeiv2_driver *
gsl_odeiv2_driver_alloc_y_new (gsl_odeiv2_system * sys,
                              const gsl_odeiv2_step_type * T,
                              const double hstart,
                              const double epsabs, const double epsrel)
{
  /* Initializes an ODE driver system with control object of type y_new. */

  gsl_odeiv2_driver *state = driver_alloc (sys, hstart, T);

  if (state == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate driver object", GSL_ENOMEM);
    }

  if (epsabs >= 0.0 && epsrel >= 0.0)
    {
      state->c = gsl_odeiv2_control_y_new (epsabs, epsrel);

      if (state->c == NULL)
        {
          gsl_odeiv2_driver_free (state);
          GSL_ERROR_NULL ("failed to allocate control object", GSL_ENOMEM);
        }
    }
  else
    {
      gsl_odeiv2_driver_free (state);
      GSL_ERROR_NULL ("epsabs and epsrel must be positive", GSL_EINVAL);
    }

  gsl_odeiv2_step_set_control (state->s, state->c);

  return state;
}

gsl_odeiv2_driver *
gsl_odeiv2_driver_alloc_scaled_new (gsl_odeiv2_system * sys,
                                   const gsl_odeiv2_step_type * T,
                                   const double hstart,
                                   const double epsabs, const double epsrel,
                                   const double a_y, const double a_dydt,
                                   const double scale_abs[])
{
  /* Initializes an ODE driver system with control object of type
     scaled_new. 
   */

  gsl_odeiv2_driver *state = driver_alloc (sys, hstart, T);

  if (state == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate driver object", GSL_ENOMEM);
    }

  if (epsabs >= 0.0 && epsrel >= 0.0)
    {
      state->c = gsl_odeiv2_control_scaled_new (epsabs, epsrel, a_y, a_dydt,
                                               scale_abs, sys->dimension);

      if (state->c == NULL)
        {
          gsl_odeiv2_driver_free (state);
          GSL_ERROR_NULL ("failed to allocate control object", GSL_ENOMEM);
        }
    }
  else
    {
      gsl_odeiv2_driver_free (state);
      GSL_ERROR_NULL ("epsabs and epsrel must be positive", GSL_EINVAL);
    }

  gsl_odeiv2_step_set_control (state->s, state->c);

  return state;
}

int
gsl_odeiv2_driver_apply (gsl_odeiv2_driver * D, double *t,
                        const double t1, double y[])
{
  /* Main driver function that evolves the system from t to t1. In
     beginning vector y contains the values of dependent variables at
     t. This function returns values at t=t1 in y. In case of
     unrecoverable error, y and t contains the values after the last
     successful step.
   */

  int s;

  D->n = 0;

  while (*t < t1)
    {
      s =
        gsl_odeiv2_evolve_apply (D->e, D->c, D->s, D->sys, t, t1, &(D->h), y);

      if (s != GSL_SUCCESS)
        {
          return s;
        }

      /* Check for maximum allowed steps */

      if ((D->nmax > 0) && (D->n > D->nmax))
        {
          return GSL_EMAXITER;
        }

      /* Set step size if maximum size is exceeded */

      if (fabs (D->h) > fabs (D->hmax))
        {
          D->h = D->hmax;
        }

      /* Check for too small step size */

      if (fabs (D->h) < fabs (D->hmin))
        {
          return GSL_ENOPROG;
        }

      D->n++;
    }

  return GSL_SUCCESS;
}

int
gsl_odeiv2_driver_reset (gsl_odeiv2_driver * D)
{
  /* Reset the driver. Resets evolve and step objects. */

  {
    int s = gsl_odeiv2_evolve_reset (D->e);

    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  {
    int s = gsl_odeiv2_step_reset (D->s);

    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  return GSL_SUCCESS;
}


void
gsl_odeiv2_driver_free (gsl_odeiv2_driver * state)
{
  if (state->c != NULL)
    {
      gsl_odeiv2_control_free (state->c);
    }

  gsl_odeiv2_evolve_free (state->e);
  gsl_odeiv2_step_free (state->s);
  free (state);
}
