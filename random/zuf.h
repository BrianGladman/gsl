#ifndef _gsl_ran_zuf_RANDOM_H_
#define _gsl_ran_zuf_RANDOM_H_

double gsl_ran_zuf_uniform(void);
void gsl_ran_zuf_seed(int);

unsigned long gsl_ran_zuf_random(void);
double gsl_ran_zuf_max();

void *gsl_ran_zuf_getRandomState(void);
void gsl_ran_zuf_setRandomState(void *);

unsigned long gsl_ran_zuf_random_wstate(void *);
double gsl_ran_zuf_uniform_wstate(void *);
void gsl_ran_zuf_seed_wstate(void *, int);

void gsl_ran_zuf_printState(void *);

#endif /* _gsl_ran_zuf_RANDOM_H_ */

