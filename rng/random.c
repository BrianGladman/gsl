#include <config.h>
#include <math.h>
#include <stdlib.h>
#include <gsl_rng.h>

/* This is the BSD random() generator. The sequence is

  
 */

inline unsigned long int random0_get (void *vstate);
inline unsigned long int random32_get (void *vstate);
inline unsigned long int random64_get (void *vstate);
inline unsigned long int random128_get (void *vstate);
inline unsigned long int random256_get (void *vstate);

double random0_get_double (void *vstate);
double random32_get_double (void *vstate);
double random64_get_double (void *vstate);
double random128_get_double (void *vstate);
double random256_get_double (void *vstate);

void glibc2random0_set (void *state, unsigned long int s);
void glibc2random32_set (void *state, unsigned long int s);
void glibc2random64_set (void *state, unsigned long int s);
void glibc2random128_set (void *state, unsigned long int s);
void glibc2random256_set (void *state, unsigned long int s);

void libc5random0_set (void *state, unsigned long int s);
void libc5random32_set (void *state, unsigned long int s);
void libc5random64_set (void *state, unsigned long int s);
void libc5random128_set (void *state, unsigned long int s);
void libc5random256_set (void *state, unsigned long int s);

void bsdrandom0_set (void *state, unsigned long int s);
void bsdrandom32_set (void *state, unsigned long int s);
void bsdrandom64_set (void *state, unsigned long int s);
void bsdrandom128_set (void *state, unsigned long int s);
void bsdrandom256_set (void *state, unsigned long int s);

void bsd_initialize (long int * x, int n, unsigned long int s);
void libc5_initialize (long int * x, int n, unsigned long int s);
void glibc2_initialize (long int * x, int n, unsigned long int s);

typedef struct
  {
    long int x;
  }
random0_state_t;

typedef struct
  {
    int i, j;
    long int x[7];
  }
random32_state_t;

typedef struct
  {
    int i, j;
    long int x[15];
  }
random64_state_t;

typedef struct
  {
    int i, j;
    long int x[31];
  }
random128_state_t;

typedef struct
  {
    int i, j;
    long int x[63];
  }
random256_state_t;

inline unsigned long int
random0_get (void *vstate)
{
  random0_state_t *state = (random0_state_t *) vstate;

  state->x = (1103515245 * state->x + 12345) & 0x7fffffffUL;
  return state->x;
}

inline long int
random_get_impl (int * i, int * j, int n, long int * x)
{
  long int k ;

  x[*i] += x[*j] ;
  k = (x[*i] >> 1) & 0x7FFFFFFF ;
  
  (*i)++ ;
  if (*i == n)
    *i = 0 ;
  
  (*j)++ ;
  if (*j == n)
    *j = 0 ;

  return k ;
}

inline unsigned long int
random32_get (void *vstate)
{
  random32_state_t *state = (random32_state_t *) vstate;
  unsigned long int k = random_get_impl (&state->i, &state->j, 7, state->x) ; 
  return k ;
}

inline unsigned long int
random64_get (void *vstate)
{
  random64_state_t *state = (random64_state_t *) vstate;
  long int k = random_get_impl (&state->i, &state->j, 15, state->x) ; 
  return k ;
}

inline unsigned long int
random128_get (void *vstate)
{
  random128_state_t *state = (random128_state_t *) vstate;
  unsigned long int k = random_get_impl (&state->i, &state->j, 31, state->x) ; 
  return k ;
}

inline unsigned long int
random256_get (void *vstate)
{
  random256_state_t *state = (random256_state_t *) vstate;
  long int k = random_get_impl (&state->i, &state->j, 63, state->x) ; 
  return k ;
}

double
random0_get_double (void *vstate)
{
  return random0_get (vstate) / 2147483648.0 ;
}

double
random32_get_double (void *vstate)
{
  return random32_get (vstate) / 2147483648.0 ;
}

double
random64_get_double (void *vstate)
{
  return random64_get (vstate) / 2147483648.0 ;
}

double
random128_get_double (void *vstate)
{
  return random128_get (vstate) / 2147483648.0 ;
}

double
random256_get_double (void *vstate)
{
  return random256_get (vstate) / 2147483648.0 ;
}

void
bsdrandom0_set (void *vstate, unsigned long int s)
{
  random0_state_t *state = (random0_state_t *) vstate;
  
  if (s == 0) 
    s = 1;

  state->x = s;
}

void
bsdrandom32_set (void *vstate, unsigned long int s)
{
  random32_state_t *state = (random32_state_t *) vstate;
  int i;

  bsd_initialize (state->x, 7, s) ;

  state->i = 3;
  state->j = 0;
  
  for (i = 0 ; i < 10 * 7 ; i++)
    random32_get (state) ; 
}

void
bsdrandom64_set (void *vstate, unsigned long int s)
{
  random64_state_t *state = (random64_state_t *) vstate;
  int i;

  bsd_initialize (state->x, 15, s) ;

  state->i = 1;
  state->j = 0;
  
  for (i = 0 ; i < 10 * 15 ; i++)
    random64_get (state) ; 
}

void
bsdrandom128_set (void *vstate, unsigned long int s)
{
  random128_state_t *state = (random128_state_t *) vstate;
  int i;

  bsd_initialize (state->x, 31, s) ;

  state->i = 3;
  state->j = 0;
  
  for (i = 0 ; i < 10 * 31 ; i++)
    random128_get (state) ; 
}

void
bsdrandom256_set (void *vstate, unsigned long int s)
{
  random256_state_t *state = (random256_state_t *) vstate;
  int i;

  bsd_initialize (state->x, 63, s) ;

  state->i = 1;
  state->j = 0;
  
  for (i = 0 ; i < 10 * 63 ; i++)
    random256_get (state) ; 
}

void 
bsd_initialize (long int * x, int n, unsigned long int s)
{
  int i; 

  if (s == 0)
    s = 1 ;

  x[0] = s;

  for (i = 1 ; i < n ; i++)
    x[i] = 1103515245 * x[i-1] + 12345 ;
}

void 
libc5_initialize (long int * x, int n, unsigned long int s)
{
  int i; 

  if (s == 0)
    s = 1 ;

  x[0] = s;

  for (i = 1 ; i < n ; i++)
    x[i] = 1103515145 * x[i-1] + 12345 ;
}

void 
glibc2_initialize (long int * x, int n, unsigned long int s)
{
  int i; 

  if (s == 0)
    s = 1 ;

  x[0] = s;

  for (i = 1 ; i < n ; i++)
    {
      const long int h = s / 127773;
      const long int t = 16807 * (s - h * 127773) - h * 2836;
      if (t < 0)
	{
	  s = t + 2147483647 ;
	}
      else
	{
	  s = t ;
	}

    x[i] = s ;
    }
}

void
glibc2random0_set (void *vstate, unsigned long int s)
{
  random0_state_t *state = (random0_state_t *) vstate;
  
  if (s == 0) 
    s = 1;

  state->x = s;
}

void
glibc2random32_set (void *vstate, unsigned long int s)
{
  random32_state_t *state = (random32_state_t *) vstate;
  int i;

  glibc2_initialize (state->x, 7, s) ;

  state->i = 3;
  state->j = 0;
  
  for (i = 0 ; i < 10 * 7 ; i++)
    random32_get (state) ; 
}

void
glibc2random64_set (void *vstate, unsigned long int s)
{
  random64_state_t *state = (random64_state_t *) vstate;
  int i;

  glibc2_initialize (state->x, 15, s) ;

  state->i = 1;
  state->j = 0;
  
  for (i = 0 ; i < 10 * 15 ; i++)
    random64_get (state) ; 
}

void
glibc2random128_set (void *vstate, unsigned long int s)
{
  random128_state_t *state = (random128_state_t *) vstate;
  int i;

  glibc2_initialize (state->x, 31, s) ;

  state->i = 3;
  state->j = 0;
  
  for (i = 0 ; i < 10 * 31 ; i++)
    random128_get (state) ; 
}

void
glibc2random256_set (void *vstate, unsigned long int s)
{
  random256_state_t *state = (random256_state_t *) vstate;
  int i;

  glibc2_initialize (state->x, 63, s) ;

  state->i = 1;
  state->j = 0;
  
  for (i = 0 ; i < 10 * 63 ; i++)
    random256_get (state) ; 
}


void
libc5random0_set (void *vstate, unsigned long int s)
{
  random0_state_t *state = (random0_state_t *) vstate;
  
  if (s == 0) 
    s = 1;

  state->x = s;
}

void
libc5random32_set (void *vstate, unsigned long int s)
{
  random32_state_t *state = (random32_state_t *) vstate;
  int i;

  libc5_initialize (state->x, 7, s) ;

  state->i = 3;
  state->j = 0;
  
  for (i = 0 ; i < 10 * 7 ; i++)
    random32_get (state) ; 
}

void
libc5random64_set (void *vstate, unsigned long int s)
{
  random64_state_t *state = (random64_state_t *) vstate;
  int i;

  libc5_initialize (state->x, 15, s) ;

  state->i = 1;
  state->j = 0;
  
  for (i = 0 ; i < 10 * 15 ; i++)
    random64_get (state) ; 
}

void
libc5random128_set (void *vstate, unsigned long int s)
{
  random128_state_t *state = (random128_state_t *) vstate;
  int i;

  libc5_initialize (state->x, 31, s) ;

  state->i = 3;
  state->j = 0;
  
  for (i = 0 ; i < 10 * 31 ; i++)
    random128_get (state) ; 
}

void
libc5random256_set (void *vstate, unsigned long int s)
{
  random256_state_t *state = (random256_state_t *) vstate;
  int i;

  libc5_initialize (state->x, 63, s) ;

  state->i = 1;
  state->j = 0;
  
  for (i = 0 ; i < 10 * 63 ; i++)
    random256_get (state) ; 
}



static const gsl_rng_type glibc2random_type =
{"glibc2random",			/* name */
 0x7fffffffUL,			/* RAND_MAX */
 0,				/* RAND_MIN */
 sizeof (random128_state_t),
 &glibc2random128_set,
 &random128_get,
 &random128_get_double};

static const gsl_rng_type glibc2random0_type =
{"glibc2random0",			/* name */
 0x7fffffffUL,			/* RAND_MAX */
 0,				/* RAND_MIN */
 sizeof (random0_state_t),
 &glibc2random0_set,
 &random0_get,
 &random0_get_double};

static const gsl_rng_type glibc2random32_type =
{"glibc2random32",			/* name */
 0x7fffffffUL,			/* RAND_MAX */
 0,				/* RAND_MIN */
 sizeof (random32_state_t),
 &glibc2random32_set,
 &random32_get,
 &random32_get_double};

static const gsl_rng_type glibc2random64_type =
{"glibc2random64",			/* name */
 0x7fffffffUL,			/* RAND_MAX */
 0,				/* RAND_MIN */
 sizeof (random64_state_t),
 &glibc2random64_set,
 &random64_get,
 &random64_get_double};

static const gsl_rng_type glibc2random128_type =
{"glibc2random128",			/* name */
 0x7fffffffUL,			/* RAND_MAX */
 0,				/* RAND_MIN */
 sizeof (random128_state_t),
 &glibc2random128_set,
 &random128_get,
 &random128_get_double};

static const gsl_rng_type glibc2random256_type =
{"glibc2random256",			/* name */
 0x7fffffffUL,			/* RAND_MAX */
 0,				/* RAND_MIN */
 sizeof (random256_state_t),
 &glibc2random256_set,
 &random256_get,
 &random256_get_double};

static const gsl_rng_type libc5random_type =
{"libc5random",			/* name */
 0x7fffffffUL,			/* RAND_MAX */
 0,				/* RAND_MIN */
 sizeof (random128_state_t),
 &libc5random128_set,
 &random128_get,
 &random128_get_double};

static const gsl_rng_type libc5random0_type =
{"libc5random0",			/* name */
 0x7fffffffUL,			/* RAND_MAX */
 0,				/* RAND_MIN */
 sizeof (random0_state_t),
 &libc5random0_set,
 &random0_get,
 &random0_get_double};

static const gsl_rng_type libc5random32_type =
{"libc5random32",			/* name */
 0x7fffffffUL,			/* RAND_MAX */
 0,				/* RAND_MIN */
 sizeof (random32_state_t),
 &libc5random32_set,
 &random32_get,
 &random32_get_double};

static const gsl_rng_type libc5random64_type =
{"libc5random64",			/* name */
 0x7fffffffUL,			/* RAND_MAX */
 0,				/* RAND_MIN */
 sizeof (random64_state_t),
 &libc5random64_set,
 &random64_get,
 &random64_get_double};

static const gsl_rng_type libc5random128_type =
{"libc5random128",			/* name */
 0x7fffffffUL,			/* RAND_MAX */
 0,				/* RAND_MIN */
 sizeof (random128_state_t),
 &libc5random128_set,
 &random128_get,
 &random128_get_double};

static const gsl_rng_type libc5random256_type =
{"libc5random256",			/* name */
 0x7fffffffUL,			/* RAND_MAX */
 0,				/* RAND_MIN */
 sizeof (random256_state_t),
 &libc5random256_set,
 &random256_get,
 &random256_get_double};

static const gsl_rng_type bsdrandom_type =
{"bsdrandom",			/* name */
 0x7fffffffUL,			/* RAND_MAX */
 0,				/* RAND_MIN */
 sizeof (random128_state_t),
 &bsdrandom128_set,
 &random128_get,
 &random128_get_double};

static const gsl_rng_type bsdrandom0_type =
{"bsdrandom0",			/* name */
 0x7fffffffUL,			/* RAND_MAX */
 0,				/* RAND_MIN */
 sizeof (random0_state_t),
 &bsdrandom0_set,
 &random0_get,
 &random0_get_double};

static const gsl_rng_type bsdrandom32_type =
{"bsdrandom32",			/* name */
 0x7fffffffUL,			/* RAND_MAX */
 0,				/* RAND_MIN */
 sizeof (random32_state_t),
 &bsdrandom32_set,
 &random32_get,
 &random32_get_double};

static const gsl_rng_type bsdrandom64_type =
{"bsdrandom64",			/* name */
 0x7fffffffUL,			/* RAND_MAX */
 0,				/* RAND_MIN */
 sizeof (random64_state_t),
 &bsdrandom64_set,
 &random64_get,
 &random64_get_double};

static const gsl_rng_type bsdrandom128_type =
{"bsdrandom128",		/* name */
 0x7fffffffUL,			/* RAND_MAX */
 0,				/* RAND_MIN */
 sizeof (random128_state_t),
 &bsdrandom128_set,
 &random128_get,
 &random128_get_double};

static const gsl_rng_type bsdrandom256_type =
{"bsdrandom256",		/* name */
 0x7fffffffUL,			/* RAND_MAX */
 0,				/* RAND_MIN */
 sizeof (random256_state_t),
 &bsdrandom256_set,
 &random256_get,
 &random256_get_double};

const gsl_rng_type *gsl_rng_libc5random    = &libc5random_type;
const gsl_rng_type *gsl_rng_libc5random0   = &libc5random0_type;
const gsl_rng_type *gsl_rng_libc5random32  = &libc5random32_type;
const gsl_rng_type *gsl_rng_libc5random64  = &libc5random64_type;
const gsl_rng_type *gsl_rng_libc5random128 = &libc5random128_type;
const gsl_rng_type *gsl_rng_libc5random256 = &libc5random256_type;

const gsl_rng_type *gsl_rng_glibc2random    = &glibc2random_type;
const gsl_rng_type *gsl_rng_glibc2random0   = &glibc2random0_type;
const gsl_rng_type *gsl_rng_glibc2random32  = &glibc2random32_type;
const gsl_rng_type *gsl_rng_glibc2random64  = &glibc2random64_type;
const gsl_rng_type *gsl_rng_glibc2random128 = &glibc2random128_type;
const gsl_rng_type *gsl_rng_glibc2random256 = &glibc2random256_type;

const gsl_rng_type *gsl_rng_bsdrandom    = &bsdrandom_type;
const gsl_rng_type *gsl_rng_bsdrandom0   = &bsdrandom0_type;
const gsl_rng_type *gsl_rng_bsdrandom32  = &bsdrandom32_type;
const gsl_rng_type *gsl_rng_bsdrandom64  = &bsdrandom64_type;
const gsl_rng_type *gsl_rng_bsdrandom128 = &bsdrandom128_type;
const gsl_rng_type *gsl_rng_bsdrandom256 = &bsdrandom256_type;




