/* blas/source_tbsv_r.h
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
    const int nonunit = (Diag == CblasNonUnit);
    size_t ix, jx;
    size_t i, j;
    size_t id;

    if (N == 0)
	return;

    /* form  x := inv( A )*x */

    if ((order == CblasRowMajor && TransA == CblasNoTrans && Uplo == CblasUpper)
        || (order == CblasColMajor && TransA == CblasTrans && Uplo == CblasLower)) {
      /* backsubstitution */
      ix = OFFSET(N, incX) + incX * (N - 1);
      if (nonunit) {
        X[ix] = X[ix] / A[lda * K + (N - 1)];
      }
      ix -= incX;
      for (i = N - 1; i > 0 && i--; ) {
        BASE tmp = X[ix];
        const size_t j_min = i + 1;
        const size_t j_max = GSL_MIN (N, i + K + 1);
        size_t jx = OFFSET (N, incX) + j_min * incX;
        for (j = j_min; j < j_max; j++) {
          const BASE Aij = A[lda * (K+i-j) + j];
          tmp -= Aij * X[jx];
          jx += incX;
        }
        if (nonunit) {
          X[ix] = tmp / A[lda * K + i];
        } else {
          X[ix] = tmp;
        }
        ix -= incX;
      }
    } else if ((order == CblasRowMajor && TransA == CblasNoTrans && Uplo == CblasLower)
               || (order == CblasColMajor && TransA == CblasTrans && Uplo == CblasUpper)) {

      /* forward substitution */
      ix = OFFSET(N, incX);
      if (nonunit) {
        X[ix] = X[ix] / A[lda * 0 + 0];
      }
      ix += incX;
      for (i = 1; i < N; i++) {
        BASE tmp = X[ix];
        const size_t j_min = (K > i ? 0 : i - K);
        const size_t j_max = i;
        size_t jx = OFFSET (N, incX) + j_min * incX;
        for (j = j_min; j < i; j++) {
          const BASE Aij = A[lda * (i-j) + j];
          tmp -= Aij * X[jx];
          jx += incX;
        }
        if (nonunit) {
          X[ix] = tmp / A[lda * 0 + i];
        } else {
          X[ix] = tmp;
        }
        ix += incX;
      }
    } else if ((order == CblasRowMajor && TransA == CblasTrans && Uplo == CblasUpper)
               || (order == CblasColMajor && TransA == CblasNoTrans && Uplo == CblasLower)) {
      
      /* form  x := inv( A' )*x */
      
      /* forward substitution */
      ix = OFFSET(N, incX);
      if (nonunit) {
        X[ix] = X[ix] / A[0 + lda * 0];
      }
      ix += incX;
      for (i = 1; i < N; i++) {
        BASE tmp = X[ix];
        const size_t j_min = (K > i ? 0 : i - K);
        const size_t j_max = i;
        size_t jx = OFFSET (N, incX) + j_min * incX;
        for (j = j_min; j < j_max; j++) {
          const BASE Aji = A[(i-j) + lda * j];
          tmp -= Aji * X[jx];
          jx += incX;
        }
        if (nonunit) {
          X[ix] = tmp / A[0 + lda * i];
        } else {
          X[ix] = tmp;
        }
        ix += incX;
      }
    } else if ((order == CblasRowMajor && TransA == CblasTrans && Uplo == CblasLower)
               || (order == CblasColMajor && TransA == CblasNoTrans && Uplo == CblasUpper)) {

      /* backsubstitution */
      ix = OFFSET(N, incX) + (N-1) * incX;
      if (nonunit) {
        X[ix] = X[ix] / A[K + (N - 1) * lda];
      }
      ix -= incX;
      for (i = N-1; i > 0 && i--;) {
        BASE tmp = X[ix];
        const size_t j_min = i + 1;
        const size_t j_max = GSL_MIN (N, i + K + 1);
        size_t jx = OFFSET (N, incX) + j_min * incX;
        for (j = j_min; j < j_max; j++) {
          const BASE Aji = A[(K+i-j) + lda * j];
          tmp -= Aji * X[jx];
          jx += incX;
        }
        if (nonunit) {
          X[ix] = tmp / A[K + lda * i];
        } else {
          X[ix] = tmp;
        }
        ix -= incX;
      }
    } else {
      BLAS_ERROR ("unrecognized operation");
    }
}
