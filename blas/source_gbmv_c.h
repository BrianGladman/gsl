/* blas/source_gbmv_c.h
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

  if(REAL0(alpha) == 0.0 && IMAG0(alpha) == 0.0 && REAL0(beta) == 1.0 && IMAG0(beta) == 0.0) return;  

  if(TransA == CblasNoTrans) {
    lenX = N;
    lenY = M;
  }
  else {
    lenX = M;
    lenY = N;
  }

  /* form  y := beta*y */
  if(! (REAL0(beta) == 1.0 && IMAG0(beta) == 0.0) ) {
    for(i=0; i<lenY; i++) {
      BASE_TYPE tmpR = REAL(Y, incY, i) * REAL0(beta) - IMAG(Y, incY, i) * IMAG0(beta);
      BASE_TYPE tmpI = REAL(Y, incY, i) * IMAG0(beta) + IMAG(Y, incY, i) * REAL0(beta);
      REAL(Y, incY, i) = tmpR;
      IMAG(Y, incY, i) = tmpI;
    }
  }

  if(REAL0(alpha) == 0.0 && IMAG0(alpha) == 0.0) return;

  if(TransA == CblasNoTrans) {
    /* form  y := alpha*A*x + y */
    for(i=0; i<lenY; i++) {
      BASE_TYPE tmpR = 0.0;
      BASE_TYPE tmpI = 0.0;
      const size_t j0 = ( i > KL ? i-KL : 0 );
      for(j=j0; j<GSL_MIN(lenX, i+KU+1); j++) {
        tmpR += REAL(X, incX, j) * REAL(A, 1, lda*i + j) - IMAG(X, incX, j) * IMAG(A, 1, lda*i + j);
	tmpI += REAL(X, incX, j) * IMAG(A, 1, lda*i + j) + IMAG(X, incX, j) * REAL(A, 1, lda*i + j);
      }
      REAL(Y, incY, i) += REAL0(alpha) * tmpR - IMAG0(alpha) * tmpI;
      IMAG(Y, incY, i) += REAL0(alpha) * tmpI + IMAG0(alpha) * tmpR;
    }
  }
  else {
    /* form  y := alpha*A'*x + y */
    for(j=0; j<lenX; j++) {
      BASE_TYPE tmpR = REAL0(alpha) * REAL(X, incX, j) - IMAG0(alpha) * IMAG(X, incX, j);
      BASE_TYPE tmpI = REAL0(alpha) * IMAG(X, incX, j) + IMAG0(alpha) * REAL(X, incX, j);
      const size_t i0 = ( j > KU ? j-KU : 0 );
      for(i=i0; i<GSL_MIN(lenY, j+KL+1); i++) {
        REAL(Y, incY, i) += tmpR * REAL(A, 1, lda*i+j) - tmpI * IMAG(A, 1, lda*i+j);
	IMAG(Y, incY, i) += tmpR * IMAG(A, 1, lda*i+j) + tmpI * REAL(A, 1, lda*i+j);
      }
    }
  }
