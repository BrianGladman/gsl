/* fft/signals.h
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

int FUNCTION(fft_signal,complex_pulse) (size_t k, 
					size_t n,
					size_t stride,
					double z_real, double z_imag,
					BASE data[],
					BASE fft[]);

int FUNCTION(fft_signal,complex_constant) (size_t n,
					   size_t stride,
					   double z_real,
					   double z_imag,
					   BASE data[],
					   BASE fft[]);

int FUNCTION(fft_signal,complex_exp) (int k,
				      size_t n,
				      size_t stride,
				      double z_real,
				      double z_imag,
				      BASE data[],
				      BASE fft[]);


int FUNCTION(fft_signal,complex_exppair) (int k1,
					  int k2,
					  size_t n,
					  size_t stride,
					  double z1_real,
					  double z1_imag,
					  double z2_real,
					  double z2_imag,
					  BASE data[],
					  BASE fft[]);

int FUNCTION(fft_signal,complex_noise) (size_t n,
					size_t stride,
					BASE data[],
					BASE fft[]);

int FUNCTION(fft_signal,real_noise) (size_t n,
				     size_t stride,
				     BASE data[],
				     BASE fft[]);

