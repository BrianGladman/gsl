/* $Id$ */
#ifndef _gsl_ran_cmrg_RANDOM_H_
#define _gsl_ran_cmrg_RANDOM_H_

double gsl_ran_cmrg_uniform(void);
void gsl_ran_cmrg_seed(int);

unsigned long gsl_ran_cmrg_random(void);
double gsl_ran_cmrg_max();

void *gsl_ran_cmrg_getRandomState(void);
void gsl_ran_cmrg_setRandomState(void *);

unsigned long gsl_ran_cmrg_random_wstate(void *);
double gsl_ran_cmrg_uniform_wstate(void *);
void gsl_ran_cmrg_seed_wstate(void *, int);

void gsl_ran_cmrg_printState(void *);

#endif /* _gsl_ran_cmrg_RANDOM_H_ */

