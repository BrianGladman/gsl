/* blas/source_symv.h
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
  size_t ix, iy, jx, jy;

  iy = 0;
  for(i=0; i<N; i++) {
    Y[iy] *= beta;
    iy += incY;
  }

  if(Uplo == CblasUpper) {
    jx = 0;
    jy = 0;
    for(j=0; j<N; j++) {
      BASE_TYPE tmp1 = alpha * X[jx];
      BASE_TYPE tmp2 = 0.0;
      Y[jy] += tmp1*A[lda*j + j];
      ix = jx;
      iy = jy;
      for(i=j+1; i<N; i++) {
        ix += incX;
        iy += incY;
        Y[iy] += tmp1 * A[lda*j + i];
	tmp2  += A[lda*j + i] * X[ix];
      }
      Y[jy] += alpha * tmp2;
      jx += incX;
      jy += incY;
    }
  }
  else {
    jx = 0;
    jy = 0;
    for(j=0; j<N; j++) {
      BASE_TYPE tmp1 = alpha * X[jx];
      BASE_TYPE tmp2 = 0.0;
      ix = 0;
      iy = 0;
      for(i=0; i<j; i++) {
        Y[iy] += tmp1 * A[lda*j + i];
	tmp2  += A[lda*j + i] * X[ix];
	ix += incX;
        iy += incY;
      }
      Y[jy] += tmp1*A[lda*j + j] + alpha * tmp2;
      jx += incX;
      jy += incY;
    }
  }
