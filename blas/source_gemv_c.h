/* blas/source_gemv_c.h
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
  size_t lenX, lenY;

  if(TransA == CblasNoTrans) {
    lenX = N;
    lenY = M;
  }
  else {
    lenX = M;
    lenY = N;
  }

  for(i=0; i<lenY; i++) {
    const BASE_TYPE tmpR = REAL(Y, incY, i) * REAL0(beta) - IMAG(Y, incY, i) * IMAG0(beta);
    const BASE_TYPE tmpI = REAL(Y, incY, i) * IMAG0(beta) - IMAG(Y, incY, i) * REAL0(beta);
    REAL(Y, incY, i) = tmpR;
    IMAG(Y, incY, i) = tmpI;
  }

  if(TransA == CblasNoTrans) {
    for(i=0; i<lenY; i++) {
      BASE_TYPE dotR = 0.0;
      BASE_TYPE dotI = 0.0;
      for(j=0; j<lenX; j++) {
        dotR += REAL(X, incX, j) * REAL(A, 1, lda*i + j) - IMAG(X, incX, j) * IMAG(A, 1, lda*i + j);
	dotI += REAL(X, incX, j) * IMAG(A, 1, lda*i + j) + IMAG(X, incX, j) * REAL(A, 1, lda*i + j);
      }
      REAL(Y, incY, i) += REAL0(alpha) * dotR - IMAG0(alpha) * dotI;
      IMAG(Y, incY, i) += REAL0(alpha) * dotI + IMAG0(alpha) * dotR;
    }
  }
  else {
    for(j=0; j<lenX; j++) {
      BASE_TYPE tmpR = REAL0(alpha) * REAL(X, incX, j) - IMAG0(alpha) * IMAG(X, incX, j);
      BASE_TYPE tmpI = REAL0(alpha) * IMAG(X, incX, j) + IMAG0(alpha) * REAL(X, incX, j);
      for(i=0; i<lenY; i++) {
        REAL(Y, incY, i) += REAL(A, 1, lda*j + i) * tmpR - IMAG(A, 1, lda*j + i) * tmpI;
	IMAG(Y, incY, i) += REAL(A, 1, lda*j + i) * tmpI + IMAG(A, 1, lda*j + i) * tmpR;
      }
    }
  }
