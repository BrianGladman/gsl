/* blas/source_axpy_c.h
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

  if (fabs(REAL0(alpha)) + fabs(IMAG0(alpha)) != 0.0 ) {
    size_t i;
    for(i=0; i<N; i++) {
      const BASE_TYPE tmpR = REAL0(alpha) * REAL(X, incX, i) - IMAG0(alpha) * IMAG(X, incX, i);
      const BASE_TYPE tmpI = REAL0(alpha) * IMAG(X, incX, i) + IMAG0(alpha) * REAL(X, incX, i);
      REAL(Y, incY, i) += tmpR;
      IMAG(Y, incY, i) += tmpI;
    }
  }
