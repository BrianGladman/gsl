/* blas/source_tXsv_c.h
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

  const int nonunit = (Diag == CblasNonUnit);
  size_t i, j;
  size_t id;

  if(N == 0) return;

  if(TransA == CblasNoTrans) {
    /* form  x := inv( A )*x */

    if(Uplo == CblasUpper) {
      /* backsubstitution */

      if(nonunit) {
        const BASE_TYPE aR = ACCESS_UP_CR(MATRIX_VAR_NAME,N,LDA,N-1,N-1);
	const BASE_TYPE aI = ACCESS_UP_CI(MATRIX_VAR_NAME,N,LDA,N-1,N-1);
	const BASE_TYPE xR = REAL(X,incX,N-1);
	const BASE_TYPE xI = IMAG(X,incX,N-1);
	const BASE_TYPE den = aR*aR + aI*aI;
	REAL(X,incX,N-1) = (xR*aR + xI*aI)/den;
	IMAG(X,incX,N-1) = (xI*aR - aI*xR)/den;
      }

      for(id=0; id<N-1; id++) {
        BASE_TYPE tmpR = REAL(X,incX,N-2-id);
	BASE_TYPE tmpI = IMAG(X,incX,N-2-id);
        i = N-2-id;
	for(j=i+1; j<GSL_MIN(N,i+KBAND+1); j++) {
	  const BASE_TYPE AijR = ACCESS_UP_CR(MATRIX_VAR_NAME,N,LDA,i,j);
	  const BASE_TYPE AijI = ACCESS_UP_CI(MATRIX_VAR_NAME,N,LDA,i,j);
	  const BASE_TYPE xR = REAL(X,incX,j);
	  const BASE_TYPE xI = IMAG(X,incX,j);
	  tmpR -= AijR*xR - AijI*xI;
	  tmpI -= AijR*xI + AijI*xR;
	}
	if(nonunit) {
          const BASE_TYPE aR = ACCESS_UP_CR(MATRIX_VAR_NAME,N,LDA,i,i);
	  const BASE_TYPE aI = ACCESS_UP_CI(MATRIX_VAR_NAME,N,LDA,i,i);
	  const BASE_TYPE den = aR*aR + aI*aI;
	  REAL(X,incX,i) = (tmpR*aR + tmpI*aI)/den;
	  IMAG(X,incX,i) = (tmpI*aR - tmpR*aI)/den;
	}
	else {
	  REAL(X,incX,i) = tmpR;
	  IMAG(X,incX,i) = tmpI;
	}
      }
    }
    else {
      /* forward substitution */

      if(nonunit) {
        const BASE_TYPE aR = ACCESS_LO_CR(MATRIX_VAR_NAME,N,LDA,0,0);
	const BASE_TYPE aI = ACCESS_LO_CI(MATRIX_VAR_NAME,N,LDA,0,0);
	const BASE_TYPE xR = REAL(X,incX,0);
	const BASE_TYPE xI = IMAG(X,incX,0);
	const BASE_TYPE den = aR*aR + aI*aI;
	REAL(X,incX,0) = (xR*aR + xI*aI)/den;
	IMAG(X,incX,0) = (xI*aR - aI*xR)/den;
      }

      for(i=1; i<N; i++) {
        BASE_TYPE tmpR = REAL(X, incX, i);
	BASE_TYPE tmpI = IMAG(X, incX, i);
	const size_t j0 = (i > KBAND ? i-KBAND : 0 );
	for(j=j0; j<i; j++) {
	  const BASE_TYPE AijR = ACCESS_LO_CR(MATRIX_VAR_NAME,N,LDA,i,j);
	  const BASE_TYPE AijI = ACCESS_LO_CI(MATRIX_VAR_NAME,N,LDA,i,j);
	  const BASE_TYPE xR = REAL(X,incX,j);
	  const BASE_TYPE xI = IMAG(X,incX,j);
	  tmpR -= AijR*xR - AijI*xI;
	  tmpI -= AijR*xI + AijI*xR;
	}
	if(nonunit) {
          const BASE_TYPE aR = ACCESS_LO_CR(MATRIX_VAR_NAME,N,LDA,i,i);
	  const BASE_TYPE aI = ACCESS_LO_CI(MATRIX_VAR_NAME,N,LDA,i,i);
	  const BASE_TYPE den = aR*aR + aI*aI;
	  REAL(X,incX,i) = (tmpR*aR + tmpI*aI)/den;
	  IMAG(X,incX,i) = (tmpI*aR - tmpR*aI)/den;
	}
	else {
	  REAL(X,incX,i) = tmpR;
	  IMAG(X,incX,i) = tmpI;
	}
      }
    }
  }
  else {
    /* form  x := inv( A' )*x */

    if(Uplo == CblasUpper) {
      /* forward substitution */

      if(nonunit) {
        const BASE_TYPE aR = ACCESS_UP_CR(MATRIX_VAR_NAME,N,LDA,0,0);
	const BASE_TYPE aI = ACCESS_UP_CI(MATRIX_VAR_NAME,N,LDA,0,0);
	const BASE_TYPE xR = REAL(X,incX,0);
	const BASE_TYPE xI = IMAG(X,incX,0);
	const BASE_TYPE den = aR*aR + aI*aI;
	REAL(X,incX,0) = (xR*aR + xI*aI)/den;
	IMAG(X,incX,0) = (xI*aR - aI*xR)/den;
      }

      for(i=1; i<N; i++) {
        BASE_TYPE tmpR = REAL(X,incX,i);
	BASE_TYPE tmpI = IMAG(X,incX,i);
	const size_t j0 = ( i > KBAND ? i-KBAND : 0 );
	for(j=j0; j<i; j++) {
	  const BASE_TYPE AijR = ACCESS_UP_CR(MATRIX_VAR_NAME,N,LDA,j,i);
	  const BASE_TYPE AijI = ACCESS_UP_CI(MATRIX_VAR_NAME,N,LDA,j,i);
	  const BASE_TYPE xR = REAL(X,incX,j);
	  const BASE_TYPE xI = IMAG(X,incX,j);
	  tmpR -= AijR*xR - AijI*xI;
	  tmpI -= AijR*xI + AijI*xR;
	}
	if(nonunit) {
          const BASE_TYPE aR = ACCESS_UP_CR(MATRIX_VAR_NAME,N,LDA,i,i);
	  const BASE_TYPE aI = ACCESS_UP_CI(MATRIX_VAR_NAME,N,LDA,i,i);
	  const BASE_TYPE den = aR*aR + aI*aI;
	  REAL(X,incX,i) = (tmpR*aR + tmpI*aI)/den;
	  IMAG(X,incX,i) = (tmpI*aR - tmpR*aI)/den;
	}
	else {
	  REAL(X,incX,i) = tmpR;
	  IMAG(X,incX,i) = tmpI;
	}
      }
    }
    else {
      /* backsubstitution */

      if(nonunit) {
        const BASE_TYPE aR = ACCESS_LO_CR(MATRIX_VAR_NAME,N,LDA,N-1,N-1);
	const BASE_TYPE aI = ACCESS_LO_CI(MATRIX_VAR_NAME,N,LDA,N-1,N-1);
	const BASE_TYPE xR = REAL(X,incX,N-1);
	const BASE_TYPE xI = IMAG(X,incX,N-1);
	const BASE_TYPE den = aR*aR + aI*aI;
	REAL(X,incX,N-1) = (xR*aR + xI*aI)/den;
	IMAG(X,incX,N-1) = (xI*aR - aI*xR)/den;
      }

      for(id=0; id<N-1; id++) {
        BASE_TYPE tmpR = REAL(X,incX,N-2-id);
	BASE_TYPE tmpI = IMAG(X,incX,N-2-id);
        i = N-2-id;
	for(j=i+1; j<GSL_MIN(N,i+KBAND+1); j++) {
	  const BASE_TYPE AijR = ACCESS_LO_CR(MATRIX_VAR_NAME,N,LDA,j,i);
	  const BASE_TYPE AijI = ACCESS_LO_CI(MATRIX_VAR_NAME,N,LDA,j,i);
	  const BASE_TYPE xR = REAL(X,incX,j);
	  const BASE_TYPE xI = IMAG(X,incX,j);
	  tmpR -= AijR*xR - AijI*xI;
	  tmpI -= AijR*xI + AijI*xR;
	}
	if(nonunit) {
          const BASE_TYPE aR = ACCESS_LO_CR(MATRIX_VAR_NAME,N,LDA,i,i);
	  const BASE_TYPE aI = ACCESS_LO_CI(MATRIX_VAR_NAME,N,LDA,i,i);
	  const BASE_TYPE den = aR*aR + aI*aI;
	  REAL(X,incX,i) = (tmpR*aR + tmpI*aI)/den;
	  IMAG(X,incX,i) = (tmpI*aR - tmpR*aI)/den;
	}
	else {
	  REAL(X,incX,i) = tmpR;
	  IMAG(X,incX,i) = tmpI;
	}
      }
    }
  }
