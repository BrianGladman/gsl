/* sys/infnan.c
 * 
 * Copyright (C) 2001 Brian Gough
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
#include <math.h>

int gsl_isnan (const double x);
int gsl_isinf (const double x);
int gsl_isreal (const double x);

int
gsl_isnan (const double x)
{
  int status = (x != x);
  return status;
}

int
gsl_isinf (const double x)
{
  if (x > DBL_MAX)
    return +1;
  else if (x < -DBL_MAX)
    return -1;
  else
    return 0;
}

int
gsl_isreal (const double x)
{
  const double y = x - x;
  int status = (y == y);
  return status;
}
