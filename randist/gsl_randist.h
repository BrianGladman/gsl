#ifndef GSL_RANDIST_H
#define GSL_RANDIST_H

int *gsl_ran_shuffle (int N, int *x);
int *gsl_ran_choose (int K, int N, int *x);

double gsl_ran_exponential (const gsl_rng * r, double mu);
double gsl_ran_cauchy (const gsl_rng * r, double mu);
double gsl_ran_laplace (const gsl_rng * r, double mu);
double gsl_ran_gaussian (const gsl_rng * r);
double gsl_ran_gamma (const gsl_rng * r, double a);
double gsl_ran_flat (const gsl_rng * r, double a, double b);

unsigned int gsl_ran_poisson (const gsl_rng * r, double mu);
void gsl_ran_poisson_array (const gsl_rng * r, size_t n, unsigned int array[],
			    double mu);

#endif /* GSL_RANDIST_H */
