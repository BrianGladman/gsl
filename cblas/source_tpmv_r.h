/* blas/source_tpmv_r.h
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
    const int nounit = (Diag == CblasNonUnit);

    if ((order == CblasRowMajor && TransA == CblasNoTrans && Uplo == CblasUpper)
	|| (order == CblasColMajor && TransA == CblasTrans && Uplo == CblasLower)) {
	/* form  x:= A*x */

      ix = OFFSET(N, incX);
      for (i = 0; i < N; i++) {
        BASE atmp = Ap[TPUP(N, i, i)];
        BASE temp = (nounit ? X[ix] * atmp : X[ix]);
        jx = OFFSET(N, incX) + (i + 1) * incX;
        for (j = i + 1; j < N; j++) {
          atmp = Ap[TPUP(N, i, j)];
          temp += atmp * X[jx];
          jx += incX;
        }
        X[ix] = temp;
        ix += incX;
      }
    } else if ((order == CblasRowMajor && TransA == CblasNoTrans && Uplo == CblasLower)
               || (order == CblasColMajor && TransA == CblasTrans && Uplo == CblasUpper)) {

      ix = OFFSET(N, incX) + (N - 1) * incX;
      for (i = N; i > 0 && i--;) {
        BASE atmp = Ap[TPLO(N, i, i)];
        BASE temp = (nounit ? X[ix] * atmp : X[ix]);
        size_t jx = OFFSET(N, incX);
        for (j = 0; j < i; j++) {
          atmp = Ap[TPLO(N, i, j)];
          temp += atmp * X[jx];
          jx += incX;
        }
        X[ix] = temp;
        ix -= incX;
      }
    } else if ((order == CblasRowMajor && TransA == CblasTrans && Uplo == CblasUpper)
	       || (order == CblasColMajor && TransA == CblasNoTrans && Uplo == CblasLower)) {
	/* form  x := A'*x */
      
      ix = OFFSET(N, incX) + (N - 1) * incX;
      for (i = N; i > 0 && i--;) {
        BASE atmp = Ap[TPUP(N, i, i)];
        BASE temp = (nounit ? X[ix] * atmp : X[ix]);
        size_t jx = OFFSET(N, incX);
        for (j = 0; j < i; j++) {
          atmp = Ap[TPUP(N, j, i)];
          temp += atmp * X[jx];
          jx += incX;
        }
        X[ix] = temp;
        ix -= incX;
      }
    } else if ((order == CblasRowMajor && TransA == CblasTrans && Uplo == CblasLower)
	       || (order == CblasColMajor && TransA == CblasNoTrans && Uplo == CblasUpper)) {
      
      ix = OFFSET(N, incX);
      for (i = 0; i < N; i++) {
        BASE atmp = Ap[TPLO(N, i, i)];
        BASE temp = (nounit ? X[ix] * atmp : X[ix]);
        jx = OFFSET(N, incX) + (i + 1) * incX;
        for (j = i + 1; j < N; j++) {
          atmp = Ap[TPLO(N, j, i)];
          temp += atmp * X[jx];
          jx += incX;
        }
        X[ix] = temp;
        ix += incX;
      }
    } else {
      BLAS_ERROR ("unrecognized operation");
    }
}

