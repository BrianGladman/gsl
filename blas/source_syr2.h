/* blas/source_syr2.h
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

  size_t i, j;
  size_t ix, iy;
  size_t jx = 0;
  size_t jy = 0;

  if(Uplo == CblasUpper) {
    for(j=0; j<N; j++) {
      const BASE_TYPE tmp1 = alpha * Y[jy];
      const BASE_TYPE tmp2 = alpha * X[jx];
      ix = jx;
      iy = jy;
      for(i=j; i<N; i++) {
        A[lda*j + i] += X[ix]*tmp1 + Y[iy]*tmp2;
	ix += incX;
	iy += incY;
      }
      jx += incX;
      jy += incY;
    }
  }
  else {
    for(j=0; j<N; j++) {
      const BASE_TYPE tmp1 = alpha * Y[jy];
      const BASE_TYPE tmp2 = alpha * X[jx];
      ix = 0;
      iy = 0;
      for(i=0; i<=j; i++) {
        A[lda*j + i] += X[ix]*tmp1 + Y[iy]*tmp2;
        ix += incX;
        iy += incY;
      }
      jx += incX;
      jy += incY;
    }
  }
