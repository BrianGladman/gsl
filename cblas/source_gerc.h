/* blas/source_gerc.h
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
  const BASE alpha_real = REAL0(alpha), alpha_imag = IMAG0(alpha);

  if (order == CblasRowMajor) {
    size_t ix = OFFSET(M, incX);
    for (i = 0; i < M; i++) {
      const BASE X_real = REAL(X, ix);
      const BASE X_imag = IMAG(X, ix);
      const BASE tmp_real = alpha_real * X_real - alpha_imag * X_imag;
      const BASE tmp_imag = alpha_imag * X_real + alpha_real * X_imag;
      size_t jy = OFFSET(N, incY);
      for (j = 0; j < N; j++) {
	const BASE Y_real = REAL(Y, jy);
	const BASE Y_imag = (-1.0) * IMAG(Y, jy);
	REAL(A, lda * i + j) += Y_real * tmp_real - Y_imag * tmp_imag;
	IMAG(A, lda * i + j) += Y_imag * tmp_real + Y_real * tmp_imag;
	jy += incY;
      }
      ix += incX;
    }
  } else if (order == CblasColMajor) {
    size_t jy = OFFSET(N, incY);
    for (j = 0; j < N; j++) {
      const BASE Y_real = REAL(Y, jy);
      const BASE Y_imag = (-1.0) * IMAG(Y, jy);
      const BASE tmp_real = alpha_real * Y_real - alpha_imag * Y_imag;
      const BASE tmp_imag = alpha_imag * Y_real + alpha_real * Y_imag;
      size_t ix = OFFSET(M, incX);
      for (i = 0; i < M; i++) {
	const BASE X_real = REAL(X, ix);
	const BASE X_imag = IMAG(X, ix);
	REAL(A, i + lda * j) += X_real * tmp_real - X_imag * tmp_imag;
	IMAG(A, i + lda * j) += X_imag * tmp_real + X_real * tmp_imag;
	ix += incX;
      }
      jy += incY;
    }
  } else {
    BLAS_ERROR("unrecognized operation");
  }
}
