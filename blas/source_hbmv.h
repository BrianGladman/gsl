/* blas/source_hbmv.h
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
  const BASE_TYPE conj = -1.0;

  for(i=0; i<N; i++) {
    const BASE_TYPE tmpR = REAL(Y, incY, i) * REAL0(beta) - IMAG(Y, incY, i) * IMAG0(beta);
    const BASE_TYPE tmpI = REAL(Y, incY, i) * IMAG0(beta) + IMAG(Y, incY, i) * REAL0(beta);
    REAL(Y, incY, i) = tmpR;
    IMAG(Y, incY, i) = tmpI;
  }

  if(Uplo == CblasUpper) {
    for(j=0; j<N; j++) {
      const BASE_TYPE tmp1R = REAL0(alpha) * REAL(X, incX, j) - IMAG0(alpha) * IMAG(X, incX, j);
      const BASE_TYPE tmp1I = REAL0(alpha) * IMAG(X, incX, j) + IMAG0(alpha) * REAL(X, incX, j);
      BASE_TYPE tmp2R = 0.0;
      BASE_TYPE tmp2I = 0.0;
      REAL(Y, incY, j) += tmp1R * REAL(A, 1, lda*j + j) /* - tmp1I * IMAG(A, 1, lda*j + j) */;
      IMAG(Y, incY, j) += /* tmp1R * IMAG(A, 1, lda*j + j) + */ tmp1I * REAL(A, 1, lda*j + j);
      for(i=j+1; i<GSL_MIN(N,j+K+1); i++) {
	REAL(Y, incY, i) +=        tmp1R * REAL(A, 1, lda*j + i) - conj * tmp1I * IMAG(A, 1, lda*j + i);
	IMAG(Y, incY, i) += conj * tmp1R * IMAG(A, 1, lda*j + i) + tmp1I * REAL(A, 1, lda*j + i);
	tmp2R += REAL(A, 1, lda*j + i) * REAL(X, incX, i) - IMAG(A, 1, lda*j + i) * IMAG(X, incX, i);
	tmp2I += REAL(A, 1, lda*j + i) * IMAG(X, incX, i) + IMAG(A, 1, lda*j + i) * REAL(X, incX, i);
      }
      REAL(Y, incY, j) += REAL0(alpha) * tmp2R - IMAG0(alpha) * tmp2I;
      IMAG(Y, incY, j) += REAL0(alpha) * tmp2I + IMAG0(alpha) * tmp2R;
    }
  }
  else {
    for(j=0; j<N; j++) {
      BASE_TYPE tmp1R = REAL0(alpha) * REAL(X, incX, j) - IMAG0(alpha) * IMAG(X, incX, j);
      BASE_TYPE tmp1I = REAL0(alpha) * IMAG(X, incX, j) + IMAG0(alpha) * REAL(X, incX, j);
      BASE_TYPE tmp2R = 0.0;
      BASE_TYPE tmp2I = 0.0;
      const size_t i_min = ( j>K ? j-K : 0 );
      for(i=i_min; i<j; i++) {
        REAL(Y, incY, i) +=        tmp1R * REAL(A, 1, lda*j + i) - conj * tmp1I * IMAG(A, 1, lda*j + i);
	IMAG(Y, incY, i) += conj * tmp1R * IMAG(A, 1, lda*j + i) + tmp1I * REAL(A, 1, lda*j + i);
	tmp2R += REAL(A, 1, lda*j + i) * REAL(X, incX, i) - IMAG(A, 1, lda*j + i) * IMAG(X, incX, i);
	tmp2I += REAL(A, 1, lda*j + i) * IMAG(X, incX, i) + IMAG(A, 1, lda*j + i) * REAL(X, incX, i);
      }
      REAL(Y, incY, j) += tmp1R * REAL(A, 1, lda*j + j) /* - tmp1I * IMAG(A, 1, lda*j + j) */;
      IMAG(Y, incY, j) += /* tmp1R  * IMAG(A, 1, lda*j + j) */ + tmp1I * REAL(A, 1, lda*j + j);
      REAL(Y, incY, j) += REAL0(alpha) * tmp2R - IMAG0(alpha) * tmp2I;
      IMAG(Y, incY, j) += REAL0(alpha) * tmp2I + IMAG0(alpha) * tmp2R;
    }
  }
