/* blas/source_hpr.h
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

    if (alpha == 0.0)
      return;

    if ((order == CblasRowMajor && Uplo == CblasUpper)
        || (order == CblasColMajor && Uplo == CblasLower)) {
      size_t ix = OFFSET(N, incX);
      for (i = 0; i < N; i++) {
	const BASE tmp_real = alpha * REAL(X, ix);
	const BASE tmp_imag = alpha * conj * IMAG(X, ix);
	size_t jx = ix;

        {
          const BASE X_real = REAL(X, jx);
          const BASE X_imag = -conj * IMAG(X, jx);
          REAL(Ap, TPUP(N,i,i)) += X_real * tmp_real - X_imag * tmp_imag;
          IMAG(Ap, TPUP(N,i,i)) = 0;
          jx += incX;
        }

	for (j = i+1 ; j < N; j++) {
          const BASE X_real = REAL(X, jx);
          const BASE X_imag = conj * IMAG(X, jx);
          REAL(Ap, TPUP(N,i,j)) += X_real * tmp_real - X_imag * tmp_imag;
          IMAG(Ap, TPUP(N,i,j)) += X_imag * tmp_real + X_real * tmp_imag;
          jx += incX;
	}
	ix += incX;
      }
    } else if ((order == CblasRowMajor && Uplo == CblasLower)
               || (order == CblasColMajor && Uplo == CblasUpper)) {
      size_t ix = OFFSET(N, incX);
      for (i = 0; i < N; i++) {
	const BASE tmp_real = alpha * REAL(X,ix);
	const BASE tmp_imag = alpha * conj * IMAG(X,ix);
	size_t jx = OFFSET(N, incX);
	for (j = 0 ; j < i; j++) {
          const BASE X_real = REAL(X, jx);
          const BASE X_imag = conj * IMAG(X, jx);
          REAL(Ap, TPLO(N,i,j)) += X_real * tmp_real - X_imag * tmp_imag;
          IMAG(Ap, TPLO(N,i,j)) += X_imag * tmp_real + X_real * tmp_imag;
          jx += incX;
	}

        {
          const BASE X_real = REAL(X, jx);
          const BASE X_imag = -conj * IMAG(X, jx);
          REAL(Ap, TPLO(N,i,i)) += X_real * tmp_real - X_imag * tmp_imag;
          IMAG(Ap, TPLO(N,i,i)) = 0;
          jx += incX;
        }

	ix += incX;
      }
    } else {
      BLAS_ERROR("unrecognized operation");
    }
}
