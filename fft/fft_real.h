/* fft/fft_real.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
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

#include "complex_internal.h"

int gsl_fft_real_pass_2 (const double in[],
			 size_t istride,
			 double out[],
			 size_t ostride,
			 const size_t product,
			 const size_t n,
			 const gsl_complex twiddle[]);

int gsl_fft_real_pass_3 (const double in[], 
			 size_t istride,
			 double out[],
			 size_t ostride,
			 const size_t product,
			 const size_t n,
			 const gsl_complex twiddle1[],
			 const gsl_complex twiddle2[]);

int gsl_fft_real_pass_4 (const double in[],
			 size_t istride,
			 double out[],
			 size_t ostride,
			 const size_t product,
			 const size_t n,
			 const gsl_complex twiddle1[],
			 const gsl_complex twiddle2[],
			 const gsl_complex twiddle3[]);

int gsl_fft_real_pass_5 (const double in[],
			 size_t istride,
			 double out[],
			 size_t ostride,
			 const size_t product,
			 const size_t n,
			 const gsl_complex twiddle1[],
			 const gsl_complex twiddle2[],
			 const gsl_complex twiddle3[],
			 const gsl_complex twiddle4[]);

int gsl_fft_real_pass_6 (const double in[], 
			 size_t istride,
			 double out[],
			 size_t ostride,
			 size_t product, size_t n,
			 gsl_complex * twiddle1, gsl_complex * twiddle2,
			 gsl_complex * twiddle3, gsl_complex * twiddle4,
			 gsl_complex * twiddle5);

int gsl_fft_real_pass_n (const double in[], 
			 size_t istride,
			 double out[],
			 size_t ostride,
			 const size_t factor,
			 const size_t product,
			 const size_t n,
			 const gsl_complex twiddle[]);
