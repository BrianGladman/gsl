#ifndef _gsl_ran_RANDOM_H_
#define _gsl_ran_RANDOM_H_

double gsl_ran_uniform(void);
void gsl_ran_seed(int);

unsigned long gsl_ran_random(void);
double gsl_ran_max();

void *gsl_ran_getRandomState(void);
void gsl_ran_setRandomState(void *);

unsigned long gsl_ran_random_wstate(void *);
double gsl_ran_uniform_wstate(void *);
void gsl_ran_seed_wstate(void *, int);

void gsl_ran_printState(void *);

#endif /* _gsl_ran_RANDOM_H_ */

