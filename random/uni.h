#ifndef _gsl_ran_uni_RANDOM_H_
#define _gsl_ran_uni_RANDOM_H_

double gsl_ran_uni_uniform(void);
void gsl_ran_uni_seed(int);

unsigned long gsl_ran_uni_random(void);
double gsl_ran_uni_max();

void *gsl_ran_uni_getRandomState(void);
void gsl_ran_uni_setRandomState(void *);

unsigned long gsl_ran_uni_random_wstate(void *);
double gsl_ran_uni_uniform_wstate(void *);
void gsl_ran_uni_seed_wstate(void *, int);

void gsl_ran_uni_printState(void *);

#endif /* _gsl_ran_uni_RANDOM_H_ */

