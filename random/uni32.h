#ifndef _gsl_ran_uni32_RANDOM_H_
#define _gsl_ran_uni32_RANDOM_H_

double gsl_ran_uni32_uniform(void);
void gsl_ran_uni32_seed(int);

unsigned long gsl_ran_uni32_random(void);
double gsl_ran_uni32_max();

void *gsl_ran_uni32_getRandomState(void);
void gsl_ran_uni32_setRandomState(void *);

unsigned long gsl_ran_uni32_random_wstate(void *);
double gsl_ran_uni32_uniform_wstate(void *);
void gsl_ran_uni32_seed_wstate(void *, int);

void gsl_ran_uni32_printState(void *);

#endif /* _gsl_ran_uni32_RANDOM_H_ */

