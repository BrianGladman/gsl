/* blas/source_syr.h
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
  size_t ix, jx;

  if(Uplo == CblasUpper) {
    jx = 0;
    for(j=0; j<N; j++) {
      BASE_TYPE tmp = alpha * X[jx];
      ix = jx;
      for(i=j; i<N; i++) {
        A[lda*j + i] += X[ix] * tmp;
        ix += incX;
      }
      jx += incX;
    }
  }
  else {
    jx = 0;
    for(j=0; j<N; j++) {
      BASE_TYPE tmp = alpha * X[jx];
      ix = 0;
      for(i=0; i<=j; i++) {
        A[lda*j + i] += X[ix]*tmp;
        ix += incX;
      }
      jx += incX;
    }
  }
