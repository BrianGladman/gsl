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
