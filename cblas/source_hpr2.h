/* blas/source_hpr2.h
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
    size_t k = 0;
    const BASE aR = REAL0(alpha);
    const BASE aI = IMAG0(alpha);
    const BASE conj = -1.0;

    if (Uplo == CblasUpper) {
	for (i = 0; i < N; i++) {
	    const BASE XiR = REAL(X, incX, i);
	    const BASE XiI = IMAG(X, incX, i);
	    const BASE YiR = REAL(Y, incY, i);
	    const BASE YiI = IMAG(Y, incY, i);
	    for (j = i; j < N; j++) {
		const BASE XjR = REAL(X, incX, j);
		const BASE XjI = IMAG(X, incX, j);
		const BASE YjR = REAL(Y, incY, j);
		const BASE YjI = IMAG(Y, incY, j);
		const BASE tmpijR = XiR * YjR - conj * XiI * YjI;
		const BASE tmpijI = XiI * YjR + conj * XiR * YjI;
		const BASE tmpjiR = XjR * YiR - conj * XjI * YiI;
		const BASE tmpjiI = XjR * YiI + conj * XjI * YiR;
		REAL(Ap, 1, k) +=
		    aR * tmpijR - aI * tmpijI + aR * tmpjiR -
		    conj * aI * tmpjiI;
		IMAG(Ap, 1, k) +=
		    aR * tmpijI + aI * tmpijR + aR * tmpjiI +
		    conj * aI * tmpjiR;
		k++;
	    }
	}
    } else {
	for (i = 0; i < N; i++) {
	    const BASE XiR = REAL(X, incX, i);
	    const BASE XiI = IMAG(X, incX, i);
	    const BASE YiR = REAL(Y, incY, i);
	    const BASE YiI = IMAG(Y, incY, i);
	    for (j = 0; j <= i; j++) {
		const BASE XjR = REAL(X, incX, j);
		const BASE XjI = IMAG(X, incX, j);
		const BASE YjR = REAL(Y, incY, j);
		const BASE YjI = IMAG(Y, incY, j);
		const BASE tmpijR = XiR * YjR - conj * XiI * YjI;
		const BASE tmpijI = XiI * YjR + conj * XiR * YjI;
		const BASE tmpjiR = XjR * YiR - conj * XjI * YiI;
		const BASE tmpjiI = XjR * YiI + conj * XjI * YiR;
		REAL(Ap, 1, k) +=
		    aR * tmpijR - aI * tmpijI + aR * tmpjiR -
		    conj * aI * tmpjiI;
		IMAG(Ap, 1, k) +=
		    aR * tmpijI + aI * tmpijR + aR * tmpjiI +
		    conj * aI * tmpjiR;
		k++;
	    }
	}
    }
}
