/* integration/qag.c
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
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include "integration.h"

int
gsl_integration_qag (const gsl_function *f,
		     double a, double b,
		     double epsabs, double epsrel, size_t limit,
		     int key,
		     gsl_integration_workspace * workspace,
		     double * result, double * abserr)
{
  int status ;
  gsl_integration_rule * integration_rule = gsl_integration_qk15 ;

  if (key < GSL_INTEG_GAUSS15)
    {
      key = GSL_INTEG_GAUSS15 ;
    } 
  else if (key > GSL_INTEG_GAUSS61) 
    {
      key = GSL_INTEG_GAUSS61 ;
    }

  switch (key) 
    {
    case GSL_INTEG_GAUSS15:
      integration_rule = gsl_integration_qk15 ;
      break ;
    case GSL_INTEG_GAUSS21:
      integration_rule = gsl_integration_qk21 ;
      break ;
    case GSL_INTEG_GAUSS31:
      integration_rule = gsl_integration_qk31 ; 
      break ;
    case GSL_INTEG_GAUSS41:
      integration_rule = gsl_integration_qk41 ;
      break ;      
    case GSL_INTEG_GAUSS51:
      integration_rule = gsl_integration_qk51 ;
      break ;      
    case GSL_INTEG_GAUSS61:
      integration_rule = gsl_integration_qk61 ;
      break ;      
    default:
      GSL_ERROR("value of key does specify a known integration rule", 
		GSL_EINVAL) ;
    }

  status = gsl_integration_qag_impl (f, a, b, epsabs, epsrel, limit,
				     workspace, 
				     result, abserr, 
				     integration_rule) ;

  return status ;
}

