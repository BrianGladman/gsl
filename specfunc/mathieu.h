/* specfunc/mathieu.h
 * 
 * Copyright (C) 2002 Lowell Johnson
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

/* Author:  L. Johnson */

#ifndef _MATHIEU_H_
#define _MATHIEU_H_

#include <gsl/gsl_sf_mathieu.h>

#define NUM_MATHIEU_COEFF 100


/* Compute an array of characteristic (eigen) values from the recurrence
   matrices for the Mathieu equations. */
int gsl_sf_mathieu_c_charv_array(double qq, gsl_sf_mathieu_workspace *work);
int gsl_sf_mathieu_s_charv_array(double qq, gsl_sf_mathieu_workspace *work);

/* Compute the characteristic value for a Mathieu function of order n and
   type ntype. */
int gsl_sf_mathieu_c_charv(int order, double qq, double *aa);
int gsl_sf_mathieu_s_charv(int order, double qq, double *aa);

/* Compute the Fourier coefficients for a Mathieu function. */
int gsl_sf_mathieu_c_coeff(int order, double qq, double aa, double coeff[]);
int gsl_sf_mathieu_s_coeff(int order, double qq, double aa, double coeff[]);

#endif /* !_MATHIEU_H_ */
