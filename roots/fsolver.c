/* roots/fsolver.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Reid Priedhorsky, Brian Gough
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

#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

gsl_root_fsolver *
gsl_root_fsolver_alloc (const gsl_root_fsolver_type * T, 
			 gsl_function * f, gsl_interval x)
{
  int status;

  gsl_root_fsolver * s = (gsl_root_fsolver *) malloc (sizeof (gsl_root_fsolver));

  if (s == 0)
    {
      GSL_ERROR_RETURN ("failed to allocate space for root solver struct",
			GSL_ENOMEM, 0);
    };

  s->state = malloc (T->size);

  if (s->state == 0)
    {
      free (s);		/* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for root solver state",
			GSL_ENOMEM, 0);
    };

  s->type = T ;

  status = gsl_root_fsolver_set (s, f, x); /* seed the generator */

  if (status != GSL_SUCCESS)
    {
      free (s->state);
      free (s);		/* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to set solver", status, 0);
    };

  return s;
}

int
gsl_root_fsolver_set (gsl_root_fsolver * s, gsl_function * f, gsl_interval x)
{
  s->function = f;
  s->root = 0;
  s->interval = x;

  if (x.lower > x.upper)
    {
      GSL_ERROR ("invalid interval (lower > upper)", GSL_EINVAL);
    }

  return (s->type->set) (s->state, s->function, &(s->root), &(s->interval));
}

int
gsl_root_fsolver_iterate (gsl_root_fsolver * s)
{
  return (s->type->iterate) (s->state, 
			     s->function, &(s->root), &(s->interval));
}

void
gsl_root_fsolver_free (gsl_root_fsolver * s)
{
  free (s->state);
  free (s);
}

const char *
gsl_root_fsolver_name (const gsl_root_fsolver * s)
{
  return s->type->name;
}

double
gsl_root_fsolver_root (const gsl_root_fsolver * s)
{
  return s->root;
}

gsl_interval
gsl_root_fsolver_interval (const gsl_root_fsolver * s)
{
  return s->interval;
}

