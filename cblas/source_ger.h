/* blas/source_ger.h
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
    size_t i, j;

    if (order == CblasRowMajor) {
      size_t ix = OFFSET(M, incX);
      for (i = 0; i < M; i++) {
	const BASE tmp = alpha * X[ix];
	size_t jy = OFFSET(N, incY);
	for (j = 0; j < N; j++) {
          A[lda * i + j] += Y[jy] * tmp;
          jy += incY;
	}
	ix += incX;
      }
    } else if (order == CblasColMajor) {
      size_t jy = OFFSET(N, incY);
      for (j = 0; j < N; j++) {
	const BASE tmp = alpha * Y[jy];
	size_t ix = OFFSET(M, incX);
	for (i = 0; i < M; i++) {
          A[i + lda * j] += X[ix] * tmp;
          ix += incX;
	}
	jy += incY;
      }
    } else {
      BLAS_ERROR ("unrecognized operation");
    }
}
