#ifndef GSL_RANDIST_H
#define GSL_RANDIST_H
#include <gsl_rng.h>

int *gsl_ran_shuffle (int N, int *x);
int *gsl_ran_choose (int K, int N, int *x);

double gsl_ran_exponential (const gsl_rng * r, double mu);
double gsl_ran_cauchy (const gsl_rng * r, double mu);
double gsl_ran_chisq (const gsl_rng * r, unsigned int n);
double gsl_ran_erlang (const gsl_rng * r, double mu, unsigned int n);
double gsl_ran_laplace (const gsl_rng * r, double mu);
double gsl_ran_gaussian (const gsl_rng * r);
double gsl_ran_lognormal (const gsl_rng * r);
double gsl_ran_logistic (const gsl_rng * r);
double gsl_ran_pareto (const gsl_rng * r, double a);
double gsl_ran_gamma (const gsl_rng * r, double a);
double gsl_ran_weibull (const gsl_rng * r, double a);
double gsl_ran_flat (const gsl_rng * r, double a, double b);

unsigned int gsl_ran_poisson (const gsl_rng * r, double mu);
unsigned int gsl_ran_geometric (const gsl_rng * r, double p);
void gsl_ran_poisson_array (const gsl_rng * r, size_t n, unsigned int array[],
			    double mu);

#endif /* GSL_RANDIST_H */
