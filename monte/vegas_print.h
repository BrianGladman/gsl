/* gsl_vegas_print.h */
/* $Id$ */

#ifndef GSL_VEGAS_PRINT_H
#define GSL_VEGAS_PRINT_H

void prn_lim(double xl[], double xu[], unsigned long dim);
void prn_head(gsl_monte_vegas_state* state, 
	      unsigned long num_dim, unsigned long calls, 
	      int it_num, int bins, int boxes);
void prn_res(int itr, double res, double err, double cum_res, double cum_err, 
	     double chi_sq);
void prn_grid(gsl_monte_vegas_state* state, unsigned long dim);
void vegas_open_log(void);
void vegas_close_log(void);


#endif
