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

/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */

{
    size_t i, j;
    const BASE conj = -1.0;

    for (j = 0; j < M; j++) {
	const BASE tmpR =
	    REAL0(alpha) * REAL(X, incX, j) - IMAG0(alpha) * IMAG(X, incX,
								  j);
	const BASE tmpI =
	    REAL0(alpha) * IMAG(X, incX, j) + IMAG0(alpha) * REAL(X, incX,
								  j);
	for (i = 0; i < N; i++) {
	    REAL(A, 1, lda * j + i) +=
		REAL(Y, incY, i) * tmpR - conj * IMAG(Y, incY, i) * tmpI;
	    IMAG(A, 1, lda * j + i) +=
		REAL(Y, incY, i) * tmpI + conj * IMAG(Y, incY, i) * tmpR;
	}
    }
}
