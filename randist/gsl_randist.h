#ifndef _GSL_RANDIST_H_
#define _GSL_RANDIST_H_

int *gsl_ran_shuffle(int N, int *x);
int *gsl_ran_choose(int K, int N, int *x);

double gsl_ran_exponential(double mu);
int gsl_ran_poisson(double mu);

void gsl_ran_poisson_array(double mu, int N, int *p);

double gsl_ran_gaussian();

double gsl_ran_gammaint(int a);
double gsl_ran_gamma(double a);

#endif /* _GSL_RANDIST_H_ */
