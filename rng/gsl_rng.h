#ifndef GSL_RNG_H
#define GSL_RNG_H

#include <stdlib.h>

typedef struct
  {
    const char *name;
    unsigned long int max;
    unsigned long int min;
    size_t size;
    void *state;
    void (*set) (void *state, unsigned long int seed);
    unsigned long int (*get) (void *state);
  }
gsl_rng;

typedef struct
  {
    const char *name;
    unsigned long int max;
    unsigned long int min;
    size_t size;
    void (*set) (void *state, unsigned long int seed);
    unsigned long int (*get) (void *state);
  }
gsl_rng_type;

/* These structs also need to appear in default.c so you can select
   them via the environment variable GSL_RNG_TYPE */

extern const gsl_rng_type *gsl_rng_bsdrand;
extern const gsl_rng_type *gsl_rng_cmrg;
extern const gsl_rng_type *gsl_rng_minstd;
extern const gsl_rng_type *gsl_rng_mrg;
extern const gsl_rng_type *gsl_rng_mt19937;
extern const gsl_rng_type *gsl_rng_r250;
extern const gsl_rng_type *gsl_rng_ran0;
extern const gsl_rng_type *gsl_rng_ran1;
extern const gsl_rng_type *gsl_rng_ran2;
extern const gsl_rng_type *gsl_rng_ran3;
extern const gsl_rng_type *gsl_rng_rand;
extern const gsl_rng_type *gsl_rng_randu;
extern const gsl_rng_type *gsl_rng_ranlux389;
extern const gsl_rng_type *gsl_rng_ranlux;
extern const gsl_rng_type *gsl_rng_taus;
extern const gsl_rng_type *gsl_rng_tt800;
extern const gsl_rng_type *gsl_rng_uni32;
extern const gsl_rng_type *gsl_rng_uni;
extern const gsl_rng_type *gsl_rng_vax;
extern const gsl_rng_type *gsl_rng_zuf;

extern const gsl_rng_type *gsl_rng_default;
extern unsigned long int gsl_rng_default_seed;

unsigned long int gsl_rng_get (const gsl_rng * r);
double gsl_rng_uniform (const gsl_rng * r);
double gsl_rng_uniform_pos (const gsl_rng * r);
double gsl_rng_uniform_gt0_lt1 (const gsl_rng * r)

gsl_rng *gsl_rng_alloc (const gsl_rng_type * T);
gsl_rng *gsl_rng_cpy (gsl_rng * dest, const gsl_rng * src);
gsl_rng *gsl_rng_clone (const gsl_rng * r);

void gsl_rng_free (gsl_rng * r);

void gsl_rng_set (const gsl_rng * r, unsigned long int seed);
unsigned long int gsl_rng_max (const gsl_rng * r);
unsigned long int gsl_rng_min (const gsl_rng * r);
const char *gsl_rng_name (const gsl_rng * r);
void gsl_rng_print_state (const gsl_rng * r);

const gsl_rng_type * gsl_rng_env_setup (void);

#ifdef HAVE_INLINE
extern inline unsigned long int
gsl_rng_get (const gsl_rng * r)
{
  return (r->get) (r->state);
}

extern inline double
gsl_rng_uniform (const gsl_rng * r)
{
  unsigned long int k = (r->get) (r->state);
  unsigned long int max = r->max;

  return k / (1.0 + max);
}

extern inline double
gsl_rng_uniform_pos (const gsl_rng * r)
{
  unsigned long int max = r->max;
  unsigned long int k;
  
  do 
    {
      k = (r->get) (r->state);
    }
  while (k == 0) ;
    
  return k / ((double) max);
}

extern inline double
gsl_rng_uniform_gt0_lt1 (const gsl_rng * r)
{
  unsigned long int max = r->max;
  unsigned long int k;
  volatile double x;

  do 
    {
      k = (r->get) (r->state);
      x = k / ((double) max)
    }
  while (x == 0 || x == 1) ;

  return x;
}
#endif /* HAVE_INLINE */

#endif /* GSL_RNG_H */
