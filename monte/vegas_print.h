/* monte/vegas_print.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Michael Booth
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

/* gsl_vegas_print.h */
/* $Id$ */

#ifndef __GSL_VEGAS_PRINT_H__
#define __GSL_VEGAS_PRINT_H__

void print_lim(gsl_monte_vegas_state* state, 
	       double xl[], double xu[], unsigned long dim);
void print_head(gsl_monte_vegas_state* state, 
		unsigned long num_dim, unsigned long calls, 
		int it_num, int bins, int boxes);
void print_res(gsl_monte_vegas_state* state, 
	       int itr, double res, double err, double cum_res, double cum_err, 
	       double chi_sq);
void print_grid(gsl_monte_vegas_state* state, unsigned long dim);

int vegas_open_log(gsl_monte_vegas_state* state);
int vegas_close_log(gsl_monte_vegas_state* state);


#endif
