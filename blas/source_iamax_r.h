/* blas/source_iamax_r.h
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

/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */

  BASE_TYPE max = 0.0;
  size_t i = 0;
  CBLAS_INDEX n;
  CBLAS_INDEX result = 0;
  for(n=0; n<N; n++) {
    if(fabs(X[i]) > max) {
      max = fabs(X[i]);
      result = n;
    }
    i += incX;
  }
  return result;
