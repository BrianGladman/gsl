/* blas/source_axpy_r.h
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

  size_t i;

  if(incX == 1 && incY == 1) {
    const size_t m = N % 4;
    for(i=0; i<m; i++) {
      Y[i] += alpha * X[i];
    }
    for(i=m; i+3<N; i += 4) {
      Y[i]   += alpha*X[i];
      Y[i+1] += alpha*X[i+1];
      Y[i+2] += alpha*X[i+2];
      Y[i+3] += alpha*X[i+3];
    }
  }
  else {
    size_t ix = 0;
    size_t iy = 0;
    for(i=0; i<N; i++) {
      Y[iy] += alpha*X[ix];
      ix += incX;
      iy += incY;
    }
  }
