/* blas/source_tbmv_r.h
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

  int nounit = ( Diag == CblasNonUnit );
  size_t i;
  size_t j;

  if(TransA == CblasNoTrans) {
    /* form  x := A*x */

    if(Uplo == CblasUpper) {
      size_t ix = 0;
      for(i=0; i<N; i++) {
        BASE_TYPE temp = 0.0;
        size_t jx = (i+1)*incX;
        for(j=i+1; j<GSL_MIN(N,i+K+1); j++) {
	  temp += X[jx] * A[lda * i + j];
	  jx += incX;
	}
	if(nounit) {
	  X[ix] = temp + X[ix] * A[lda * i + i];
	}
	else {
	  X[ix] += temp;
	}
	ix += incX;
      }
    }
    else {
      size_t ix = (N-1)*incX;
      for(i=0; i<N; i++) {
        BASE_TYPE temp = 0.0;
	const size_t j_min = ( K>N-1-i ? 0 : N-1-i-K );
	size_t jx = j_min * incX;
        for(j=j_min; j<N-1-i; j++) {
	  temp += X[jx] * A[lda * (N-1-i) + j];
	  jx += incX;
	}
	if(nounit) {
	  X[ix] = temp + X[ix] * A[lda * (N-1-i) + N-1-i];
	}
	else {
	  X[ix] += temp;
	}
	ix -= incX;
      }
    }
  }
  else {
    /* form  x := A'*x */

    if(Uplo == CblasUpper) {
      size_t ix = (N-1)*incX;
      for(i=0; i<N; i++) {
        BASE_TYPE temp = 0.0;
	const size_t j_min = ( K>N-1-i ? 0 : N-1-i-K );
	size_t jx = j_min * incX;
        for(j=j_min; j<N-1-i; j++) {
          temp += X[jx] * A[lda * j + N-1-i];
          jx += incX;
        }
	if(nounit) {
          X[ix] = temp + X[ix] * A[lda * (N-1-i) + (N-1-i)];
        }
	else {
          X[ix] += temp;
        }
        ix -= incX;
      }
    }
    else {
      size_t ix = 0;
      for(i=0; i<N; i++) {
        BASE_TYPE temp = 0.0;
        size_t jx = (i+1)*incX;
        for(j=i+1; j<GSL_MIN(N,i+K+1); j++) {
	  temp += X[jx] * A[lda * j + i];
	  jx += incX;
	}
	if(nounit) {
	  X[ix] = temp + X[ix] * A[lda * i + i];
	}
	else {
	  X[ix] += temp;
	}
	ix += incX;
      }
    }
  }
