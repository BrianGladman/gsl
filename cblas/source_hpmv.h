/* blas/source_hpmv.h
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
    size_t i, j, k;
    size_t kk = 0;
    const BASE conj = -1.0;
    const BASE aR = REAL0(alpha);
    const BASE aI = IMAG0(alpha);
    const BASE bR = REAL0(beta);
    const BASE bI = IMAG0(beta);

    for (i = 0; i < N; i++) {
	const BASE yR = REAL(Y, incY, i);
	const BASE yI = IMAG(Y, incY, i);
	REAL(Y, incY, i) = yR * bR - yI * bI;
	IMAG(Y, incY, i) = yR * bI + yI * bR;
    }

    if (Uplo == CblasUpper) {
	for (j = 0; j < N; j++) {
	    BASE tmp1R = aR * REAL(X, incX, j) - aI * IMAG(X, incX, j);
	    BASE tmp1I = aR * IMAG(X, incX, j) + aI * REAL(X, incX, j);
	    BASE tmp2R = 0.0;
	    BASE tmp2I = 0.0;
	    REAL(Y, incY, j) +=
		tmp1R * REAL(Ap, 1, kk) - tmp1I * IMAG(Ap, 1, kk);
	    IMAG(Y, incY, j) +=
		tmp1R * IMAG(Ap, 1, kk) + tmp1I * REAL(Ap, 1, kk);
	    i = j;
	    for (k = kk + 1; k < kk + N - j; k++) {
		const BASE apkR = REAL(Ap, 1, k);
		const BASE apkI = IMAG(Ap, 1, k);
		i++;
		REAL(Y, incY, i) += tmp1R * apkR - tmp1I * conj * apkI;
		IMAG(Y, incY, i) += tmp1I * apkR + tmp1R * conj * apkI;
		tmp2R += apkR * REAL(X, incX, i) - apkI * IMAG(X, incX, i);
		tmp2I += apkR * IMAG(X, incX, i) + apkI * REAL(X, incX, i);
	    }
	    REAL(Y, incY, j) += aR * tmp2R - aI * tmp2I;
	    IMAG(Y, incY, j) += aR * tmp2I + aI * tmp2R;
	    kk += N - j;
	}
    } else {
	for (j = 0; j < N; j++) {
	    BASE tmp1R = aR * REAL(X, incX, j) - aI * IMAG(X, incX, j);
	    BASE tmp1I = aR * IMAG(X, incX, j) + aI * REAL(X, incX, j);
	    BASE tmp2R = 0.0;
	    BASE tmp2I = 0.0;
	    i = 0;
	    for (k = kk; k < kk + j; k++) {
		const BASE apkR = REAL(Ap, 1, k);
		const BASE apkI = IMAG(Ap, 1, k);
		REAL(Y, incY, i) += tmp1R * apkR - tmp1I * conj * apkI;
		IMAG(Y, incY, i) += tmp1I * apkR + tmp1R * conj * apkI;
		tmp2R += apkR * REAL(X, incX, i) - apkI * IMAG(X, incX, i);
		tmp2I += apkR * IMAG(X, incX, i) + apkI * REAL(X, incX, i);
		i++;
	    }
	    REAL(Y, incY, j) +=
		tmp1R * REAL(Ap, 1, kk + j) - tmp1I * IMAG(Ap, 1, kk + j);
	    IMAG(Y, incY, j) +=
		tmp1R * IMAG(Ap, 1, kk + j) + tmp1I * REAL(Ap, 1, kk + j);
	    REAL(Y, incY, j) += aR * tmp2R - aI * tmp2I;
	    IMAG(Y, incY, j) += aI * tmp2R + aR * tmp2I;
	    kk += j + 1;
	}
    }
}
