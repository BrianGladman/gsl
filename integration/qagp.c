/* integration/qagp.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include "integration.h"

int
gsl_integration_qagp (const gsl_function *f,
		      double * pts, size_t npts,
		      double epsabs, double epsrel, size_t limit,
		      gsl_integration_workspace * workspace,
		      double * result, double * abserr)
{
  int status = gsl_integration_qagp_impl (f, pts, npts, 
					  epsabs, epsrel, limit,
					  workspace,
					  result, abserr,
					  &gsl_integration_qk21) ;
  
  return status ;
}

