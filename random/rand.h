/* $Id$ */
#ifndef _gsl_ran_rand_RANDOM_H_
#define _gsl_ran_rand_RANDOM_H_

double gsl_ran_rand_uniform(void);
void gsl_ran_rand_seed(int);

unsigned long gsl_ran_rand_random(void);
double gsl_ran_rand_max();

void *gsl_ran_rand_getRandomState(void);
void gsl_ran_rand_setRandomState(void *);

unsigned long gsl_ran_rand_random_wstate(void *);
double gsl_ran_rand_uniform_wstate(void *);
void gsl_ran_rand_seed_wstate(void *, int);

void gsl_ran_rand_printState(void *);

#endif /* _gsl_ran_rand_RANDOM_H_ */

