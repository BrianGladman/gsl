/* fit/linear.c
 * 
 * Copyright (C) 2000 Brian Gough
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
#include <gsl_fit.h>

/* Fit the weighted data (x_i, w_i, y_i) to the linear relationship 

   Y = c0 + c1 x

   returning, 

   c0, c1  --  coefficients
   s0, s1  --  the standard deviations of c0 and c1,
   r       --  the correlation coefficient between c0 and c1,
   chisq   --  weighted sum of squares of residuals */

int 
gsl_fit_linear (const double * x, 
                const double * w,  const double * y, 
                size_t n,
                double * c0, double * c1, 
                double * s0, double * s1, double * r, 
                double * chisq)
{
  
}


