/* integration/integration.h
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

void gsl_integration_qk (const int n, const double xgk[], 
                         const double wg[], const double wgk[],
                         double fv1[], double fv2[],
                         const gsl_function *f, double a, double b,
                         double * result, double * abserr, 
                         double * resabs, double * resasc);

int
gsl_integration_qag_impl (const gsl_function *f,
			  const double a, const double b,
			  const double epsabs, const double epsrel,
			  const size_t limit,
			  gsl_integration_workspace * workspace,
			  double * result, double * abserr,
			  gsl_integration_rule * q) ;

void
gsl_integration_qelg (size_t * n, double epstab[], 
		      double * result, double * abserr,
		      double res3la[], size_t * nres) ;

void
gsl_integration_qcheb (gsl_function * f, double a, double b, 
		       double *cheb12, double *cheb24);

int
gsl_integration_qags_impl (const gsl_function *f, 
			   const double a, const double b, 
			   const double epsabs, const double epsrel, const size_t limit,
			   gsl_integration_workspace * workspace,
			   double * result, double * abserr,
			   gsl_integration_rule * q) ;

int
gsl_integration_qagp_impl (const gsl_function *f,
			   const double *pts, const size_t npts,
			   const double epsabs, const double epsrel, const size_t limit,
			   gsl_integration_workspace * workspace,
			   double *result, double *abserr,
			   gsl_integration_rule * q);

void gsl_integration_qcheb (gsl_function * f, double a, double b, 
			    double *cheb12, double *cheb24);

