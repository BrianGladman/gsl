/* blas/source_nrm2_c.h
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
    BASE scale = 0.0;
    BASE ssq = 1.0;
    size_t n;
    size_t i = 0;

    if (N == 0 || incX < 1) {
	return 0;
    }

    for (n = 0; n < N; n++) {
	const BASE xi = REAL(X, i);
	const BASE yi = IMAG(X, i);

	if (xi != 0) {
	    const BASE axi = fabs(xi);

	    if (scale < axi) {
		ssq = 1.0 + ssq * (scale / axi) * (scale / axi);
		scale = axi;
	    } else {
		ssq += (axi / scale) * (axi / scale);
	    }
	}

	if (yi != 0) {
	    const BASE ayi = fabs(yi);

	    if (scale < ayi) {
		ssq = 1.0 + ssq * (scale / ayi) * (scale / ayi);
		scale = ayi;
	    } else {
		ssq += (ayi / scale) * (ayi / scale);
	    }
	}

	i += incX;
    }

    return scale * sqrt(ssq);
}
