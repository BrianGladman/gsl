/* blas/source_her2.h
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

{
  const int conj = (order == CblasColMajor) ? -1 : 1;
    size_t i, j;
    size_t ix, jx;

    const BASE alpha_real = REAL0(alpha), alpha_imag = IMAG0(alpha);

    if (alpha_real == 0.0 && alpha_imag == 0.0)
      return;

    if ((order == CblasRowMajor && Uplo == CblasUpper)
        || (order == CblasColMajor && Uplo == CblasLower)) {
      size_t ix = OFFSET(N, incX);
      size_t iy = OFFSET(N, incY);
      for (i = 0; i < N; i++) {
        const BASE Xi_real = REAL(X, ix);
        const BASE Xi_imag = IMAG(X, ix);
        /* tmp1 = alpha Xi */
	const BASE tmp1_real = alpha_real * Xi_real - alpha_imag * Xi_imag;
	const BASE tmp1_imag = alpha_imag * Xi_real + alpha_real * Xi_imag;

        const BASE Yi_real = REAL(Y, iy);
        const BASE Yi_imag = IMAG(Y, iy);
        /* tmp2 = conj(alpha) Yi */
	const BASE tmp2_real = alpha_real * Yi_real + alpha_imag * Yi_imag;
	const BASE tmp2_imag = -alpha_imag * Yi_real + alpha_real * Yi_imag;

	size_t jx = ix + incX;
	size_t jy = iy + incY;

        /* Aij = alpha*Xi*conj(Yj) + conj(alpha)*Yi*conj(Xj) */

        REAL(A, lda * i + i) += 2*(tmp1_real * Yi_real + tmp1_imag * Yi_imag);
        IMAG(A, lda * i + i) = 0;

	for (j = i+1 ; j < N; j++) {
          const BASE Xj_real = REAL(X, jx);
          const BASE Xj_imag = IMAG(X, jx);
          const BASE Yj_real = REAL(Y, jy);
          const BASE Yj_imag = IMAG(Y, jy);
          REAL(A, lda * i + j) += ((tmp1_real * Yj_real + tmp1_imag * Yj_imag)
                                   + (tmp2_real * Xj_real + tmp2_imag * Xj_imag));
          IMAG(A, lda * i + j) += conj * ((tmp1_imag * Yj_real - tmp1_real * Yj_imag) 
                                          + (tmp2_imag * Xj_real - tmp2_real * Xj_imag));
          jx += incX;
          jy += incY;
	}
	ix += incX;
        iy += incY;
      }
    } else if ((order == CblasRowMajor && Uplo == CblasLower)
               || (order == CblasColMajor && Uplo == CblasUpper)) {

      size_t ix = OFFSET(N, incX);
      size_t iy = OFFSET(N, incY);
      for (i = 0; i < N; i++) {
        const BASE Xi_real = REAL(X, ix);
        const BASE Xi_imag = IMAG(X, ix);
	const BASE tmp1_real = alpha_real * Xi_real - alpha_imag * Xi_imag;
	const BASE tmp1_imag = alpha_imag * Xi_real + alpha_real * Xi_imag;

        const BASE Yi_real = REAL(Y, iy);
        const BASE Yi_imag = IMAG(Y, iy);
	const BASE tmp2_real = alpha_real * Yi_real + alpha_imag * Yi_imag;
	const BASE tmp2_imag = -alpha_imag * Yi_real + alpha_real * Yi_imag;

	size_t jx = OFFSET(N, incX);
	size_t jy = OFFSET(N, incY);

        /* Aij = alpha*Xi*conj(Yj) + conj(alpha)*Yi*conj(Xj) */

	for (j = 0 ; j < i; j++) {
          const BASE Xj_real = REAL(X, jx);
          const BASE Xj_imag = IMAG(X, jx);
          const BASE Yj_real = REAL(Y, jy);
          const BASE Yj_imag = IMAG(Y, jy);
          REAL(A, lda * i + j) += ((tmp1_real * Yj_real + tmp1_imag * Yj_imag)
                                   + (tmp2_real * Xj_real + tmp2_imag * Xj_imag));
          IMAG(A, lda * i + j) += conj * ((tmp1_imag * Yj_real - tmp1_real * Yj_imag) 
                                          + (tmp2_imag * Xj_real - tmp2_real * Xj_imag));
          jx += incX;
          jy += incY;
	}
        
        REAL(A, lda * i + i) += 2 * (tmp1_real * Yi_real + tmp1_imag * Yi_imag);
        IMAG(A, lda * i + i) = 0;

        ix += incX;
	iy += incY;
      }
    } else {
      BLAS_ERROR("unrecognized operation");
    }
}
