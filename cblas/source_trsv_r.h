/* blas/source_tXsv_r.h
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
    const int nonunit = (Diag == CblasNonUnit);
    size_t ix, jx;
    size_t i, j;
    size_t id;

    if (N == 0)
	return;

    if (TransA == CblasNoTrans) {
	/* form  x := inv( A )*x */

	if (Uplo == CblasUpper) {
	    /* backsubstitution */

	    if (nonunit) {
		const size_t max_ix = incX * (N - 1);
		X[max_ix] =
		    X[max_ix] / ACCESS_UP(MATRIX_VAR_NAME, N, LDA, N - 1,
					  N - 1);
	    }

	    ix = incX * (N - 2);
	    for (id = 0; id < N - 1; id++) {
		BASE tmp = X[ix];
		i = N - 2 - id;
		jx = ix + incX;
		for (j = i + 1; j < GSL_MIN(N, i + KBAND + 1); j++) {
		    const BASE Aij =
			ACCESS_UP(MATRIX_VAR_NAME, N, LDA, i, j);
		    tmp -= Aij * X[jx];
		    jx += incX;
		}
		if (nonunit) {
		    X[ix] = tmp / ACCESS_UP(MATRIX_VAR_NAME, N, LDA, i, i);
		} else {
		    X[ix] = tmp;
		}
		ix -= incX;
	    }
	} else {
	    /* forward substitution */

	    if (nonunit) {
		X[0] = X[0] / ACCESS_LO(MATRIX_VAR_NAME, N, LDA, 0, 0);
	    }

	    ix = incX;
	    for (i = 1; i < N; i++) {
		BASE tmp = X[ix];
		const size_t j0 = (i > KBAND ? i - KBAND : 0);
		jx = j0 * incX;
		for (j = j0; j < i; j++) {
		    const BASE Aij =
			ACCESS_LO(MATRIX_VAR_NAME, N, LDA, i, j);
		    tmp -= Aij * X[jx];
		    jx += incX;
		}
		if (nonunit) {
		    X[ix] = tmp / ACCESS_LO(MATRIX_VAR_NAME, N, LDA, i, i);
		} else {
		    X[ix] = tmp;
		}
		ix += incX;
	    }
	}
    } else {
	/* form  x := inv( A' )*x */

	if (Uplo == CblasUpper) {
	    /* forward substitution */

	    if (nonunit) {
		X[0] = X[0] / ACCESS_UP(MATRIX_VAR_NAME, N, LDA, 0, 0);
	    }

	    ix = incX;
	    for (i = 1; i < N; i++) {
		BASE tmp = X[ix];
		const size_t j0 = (i > KBAND ? i - KBAND : 0);
		jx = j0 * incX;
		for (j = j0; j < i; j++) {
		    const BASE Aji =
			ACCESS_UP(MATRIX_VAR_NAME, N, LDA, j, i);
		    tmp -= Aji * X[jx];
		    jx += incX;
		}
		if (nonunit) {
		    X[ix] = tmp / ACCESS_UP(MATRIX_VAR_NAME, N, LDA, i, i);
		} else {
		    X[ix] = tmp;
		}
		ix += incX;
	    }
	} else {
	    /* backsubstitution */

	    if (nonunit) {
		const size_t max_ix = incX * (N - 1);
		X[max_ix] =
		    X[max_ix] / ACCESS_LO(MATRIX_VAR_NAME, N, LDA, N - 1,
					  N - 1);
	    }

	    ix = incX * (N - 2);
	    for (id = 0; id < N - 1; id++) {
		BASE tmp = X[ix];
		i = N - 2 - id;
		jx = ix + incX;
		for (j = i + 1; j < GSL_MIN(N, i + KBAND + 1); j++) {
		    const BASE Aji =
			ACCESS_LO(MATRIX_VAR_NAME, N, LDA, j, i);
		    tmp -= Aji * X[jx];
		    jx += incX;
		}
		if (nonunit) {
		    X[ix] = tmp / ACCESS_LO(MATRIX_VAR_NAME, N, LDA, i, i);
		} else {
		    X[ix] = tmp;
		}
		ix -= incX;
	    }
	}
    }
}
