/* gsl_vegas_print.h */
/* $Id$ */

#ifndef GSL_VEGAS_PRINT_H
#define GSL_VEGAS_PRINT_H

void prn_lim(double a[], double b[], int m);
void prn_head(gsl_monte_vegas_state* state, 
	      int num_dim, int calls, int it_num, int bins, int boxes);
void prn_res(int a, double b, double c, double d, double e, double f);
void prn_grid(gsl_monte_vegas_state* state, int m);
void vegas_open_log(void);
void vegas_close_log(void);


#endif
