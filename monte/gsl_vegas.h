/* header for the gsl "vegas" routines.  Mike Booth, May 1998 */
/* RCS $Id$ */

#ifndef GSL_VEGAS_H
#define GSL_VEGAS_H

/* This will go away soon. */
#define GSL_V_BINS_MAX 50  /* even integer because will be divided by two. */
#define GSL_V_MAX_DIM 10

/* control variables */
extern double acc, alpha;
extern int    mode, verbose;
extern int    max_it_num;
extern int    calls;

/* function prototypes */
typedef double (*gsl_monte_f_T)(double *);


int gsl_monte_vegas(gsl_monte_f_T fxn, double xl[], double xu[], int num_dim,
		    double* tot_int, double* tot_sig, double* chi_sq_ptr);
int gsl_monte_vegas1(gsl_monte_f_T fxn, double xl[], double xu[], int num_dim,
		     double* tot_int, double* tot_sig, double* chi_sq_ptr);
int gsl_monte_vegas2(gsl_monte_f_T fxn, double xl[], double xu[], int num_dim,
		     double* tot_int, double* tot_sig, double* chi_sq_ptr);
int gsl_monte_vegas3(gsl_monte_f_T fxn, double xl[], double xu[], int num_dim,
		     double* tot_int, double* tot_sig, double* chi_sq_ptr);

#endif /* !GSL_VEGAS_H */

