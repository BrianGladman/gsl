/* monte/init_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Michael Booth
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

BASE * 
FUNCTION(gsl_monte_vector,alloc) (const size_t n)
{
  BASE * v ;

  if (n == 0)
    {
      GSL_ERROR_RETURN ("vector length n must be positive integer", 
			GSL_EDOM, 0) ;
    }

  v = (BASE *) malloc(n * sizeof(BASE)) ;

  if (v == 0) 
    {
      GSL_ERROR_RETURN ("failed to allocate space for vector data", 
			GSL_ENOMEM, 0);
    }
  return v ;
}

BASE *
FUNCTION(gsl_monte_vector,calloc) (const size_t n)
{
  size_t i ;

  BASE * v = (BASE *) FUNCTION(gsl_monte_vector,alloc) (n) ;
  
  if (v == 0) 
    return 0 ;

  for (i = 0 ; i < n; i++)  /* initialize vector to zero */
    {
      v[i] = ZERO ;
    }

  return v ;
}


int
FUNCTION(gsl_monte_vector,free) (BASE * v)
{
  if ( v == (BASE *) NULL) {
    GSL_ERROR("Attempt to free null pointer", GSL_EFAULT);
  }
  free(v) ;
  return 0;
}

