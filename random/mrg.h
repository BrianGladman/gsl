#ifndef _gsl_ran_mrg_RANDOM_H_
#define _gsl_ran_mrg_RANDOM_H_

double gsl_ran_mrg_uniform(void);
void gsl_ran_mrg_seed(int);

unsigned long gsl_ran_mrg_random(void);
double gsl_ran_mrg_max();

void *gsl_ran_mrg_getRandomState(void);
void gsl_ran_mrg_setRandomState(void *);

unsigned long gsl_ran_mrg_random_wstate(void *);
double gsl_ran_mrg_uniform_wstate(void *);
void gsl_ran_mrg_seed_wstate(void *, int);

void gsl_ran_mrg_printState(void *);

#endif /* _gsl_ran_mrg_RANDOM_H_ */

