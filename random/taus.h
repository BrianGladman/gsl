#ifndef _gsl_ran_taus_RANDOM_H_
#define _gsl_ran_taus_RANDOM_H_

double gsl_ran_taus_uniform(void);
void gsl_ran_taus_seed(int);

unsigned long gsl_ran_taus_random(void);
double gsl_ran_taus_max();

void *gsl_ran_taus_getRandomState(void);
void gsl_ran_taus_setRandomState(void *);

unsigned long gsl_ran_taus_random_wstate(void *);
double gsl_ran_taus_uniform_wstate(void *);
void gsl_ran_taus_seed_wstate(void *, int);

void gsl_ran_taus_printState(void *);

#endif /* _gsl_ran_taus_RANDOM_H_ */

