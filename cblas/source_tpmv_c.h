/* blas/source_tpmv_c.h
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


{
    size_t i, j;
    size_t ix, jx;
    const int conj = (TransA == CblasConjTrans) ? -1 : 1;
    const int Trans = (TransA != CblasConjTrans) ? TransA : CblasTrans;
    const int nounit = (Diag == CblasNonUnit);

    if ((order == CblasRowMajor && Trans == CblasNoTrans && Uplo == CblasUpper)
	|| (order == CblasColMajor && Trans == CblasTrans && Uplo == CblasLower)) {
	/* form  x:= A*x */
      
      ix = OFFSET(N, incX);
      for (i = 0; i < N; i++) {
        BASE atmp_r = REAL(Ap, TPUP(N, i, i));
        BASE atmp_i = conj * IMAG(Ap, TPUP(N, i, i));
        BASE temp_r;
        BASE temp_i;
        if (nounit) {
          BASE x_real = REAL(X, ix);
          BASE x_imag = IMAG(X, ix);
          temp_r = atmp_r * x_real - atmp_i * x_imag;
          temp_i = atmp_r * x_imag + atmp_i * x_real;
        } else {
          temp_r = REAL(X, ix);
          temp_i = IMAG(X, ix);
        }
        
        jx = OFFSET(N, incX) + (i + 1) * incX;
        for (j = i + 1; j < N; j++) {
          BASE atmp_r = REAL(Ap, TPUP(N, i, j));
          BASE atmp_i = conj * IMAG(Ap, TPUP(N, i, j));
          BASE x_real = REAL(X, jx);
          BASE x_imag = IMAG(X, jx);
          temp_r += atmp_r * x_real - atmp_i * x_imag;
          temp_i += atmp_r * x_imag + atmp_i * x_real;
          jx += incX;
        }
        
        REAL(X, ix) = temp_r;
        IMAG(X, ix) = temp_i;
        ix += incX;
      }
    } else if ((order == CblasRowMajor && Trans == CblasNoTrans && Uplo == CblasLower)
               || (order == CblasColMajor && Trans == CblasTrans && Uplo == CblasUpper)) {

      ix = OFFSET(N, incX) + incX * (N - 1);
      for (i = N;  i > 0 && i--;) {
        BASE atmp_r = REAL(Ap, TPLO(N, i, i));
        BASE atmp_i = conj * IMAG(Ap, TPLO(N, i, i));
        BASE temp_r;
        BASE temp_i;
        if (nounit) {
          BASE x_real = REAL(X, ix);
          BASE x_imag = IMAG(X, ix);
          temp_r = atmp_r * x_real - atmp_i * x_imag;
          temp_i = atmp_r * x_imag + atmp_i * x_real;
        } else {
          temp_r = REAL(X, ix);
          temp_i = IMAG(X, ix);
        }
        
        jx = OFFSET(N, incX);
        for (j = 0; j < i; j++) {
          BASE atmp_r = REAL(Ap, TPLO(N, i, j));
          BASE atmp_i = conj * IMAG(Ap, TPLO(N, i, j));
          BASE x_real = REAL(X, jx);
          BASE x_imag = IMAG(X, jx);
          temp_r += atmp_r * x_real - atmp_i * x_imag;
          temp_i += atmp_r * x_imag + atmp_i * x_real;
          jx += incX;
        }
        
        REAL(X, ix) = temp_r;
        IMAG(X, ix) = temp_i;
        ix -= incX;
      }
    }
    else if ((order == CblasRowMajor && Trans == CblasTrans && Uplo == CblasUpper)
	       || (order == CblasColMajor && Trans == CblasNoTrans && Uplo == CblasLower)) {
	/* form  x := A'*x */
      
      ix = OFFSET(N, incX) + incX * (N - 1);
      for (i = N; i > 0 && i--;) {
        BASE atmp_r = REAL(Ap, TPUP(N, i, i));
        BASE atmp_i = conj * IMAG(Ap, TPUP(N, i, i));
        BASE temp_r;
        BASE temp_i;
        if (nounit) {
          BASE x_real = REAL(X, ix);
          BASE x_imag = IMAG(X, ix);
          temp_r = atmp_r * x_real - atmp_i * x_imag;
          temp_i = atmp_r * x_imag + atmp_i * x_real;
        } else {
          temp_r = REAL(X, ix);
          temp_i = IMAG(X, ix);
        }
        jx = OFFSET(N, incX);
        for (j = 0; j < i; j++) {
          BASE x_real = REAL(X, jx);
          BASE x_imag = IMAG(X, jx);
          BASE atmp_r = REAL(Ap, TPUP(N, j, i));
          BASE atmp_i = conj * IMAG(Ap, TPUP(N, j, i));
          temp_r += atmp_r * x_real - atmp_i * x_imag;
          temp_i += atmp_r * x_imag + atmp_i * x_real;
          jx += incX;
        }
        
        REAL(X, ix) = temp_r;
        IMAG(X, ix) = temp_i;
        ix -= incX;
      }
    } else if ((order == CblasRowMajor && Trans == CblasTrans && Uplo == CblasLower)
	       || (order == CblasColMajor && Trans == CblasNoTrans && Uplo == CblasUpper)) {

      ix = OFFSET(N, incX);
      for (i = 0; i < N; i++) {
        BASE atmp_r = REAL(Ap, TPLO(N, i, i));
        BASE atmp_i = conj * IMAG(Ap, TPLO(N, i, i));
        BASE temp_r;
        BASE temp_i;
        if (nounit) {
          BASE x_real = REAL(X, ix);
          BASE x_imag = IMAG(X, ix);
          temp_r = atmp_r * x_real - atmp_i * x_imag;
          temp_i = atmp_r * x_imag + atmp_i * x_real;
        } else {
          temp_r = REAL(X, ix);
          temp_i = IMAG(X, ix);
        }
        jx = OFFSET(N, incX) + (i + 1) * incX;
        for (j = i + 1; j < N; j++) {
          BASE x_real = REAL(X, jx);
          BASE x_imag = IMAG(X, jx);
          BASE atmp_r = REAL(Ap, TPLO(N, j, i));
          BASE atmp_i = conj * IMAG(Ap, TPLO(N, j, i));
          temp_r += atmp_r * x_real - atmp_i * x_imag;
          temp_i += atmp_r * x_imag + atmp_i * x_real;
          jx += incX;
        }
        REAL(X, ix) = temp_r;
        IMAG(X, ix) = temp_i;
        ix += incX;
      }
    } else {
	BLAS_ERROR("unrecognized operation");
    }
}
