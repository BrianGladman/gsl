/* Author: G. Jungman
 * RCS:    $Id$
 */
#include "gsl_qrng.h"


#define NIED2_MAX_DIMENSION 12


static size_t nied2_state_size(unsigned int dimension);
static int nied2_init(void * state, unsigned int dimension);
static int nied2_get(void * state, unsigned int dimension, double * v);


typedef struct
{
} nied2_state_t;


static size_t nied2_state_size(unsigned int dimension)
{
  return 0;
}


static int nied2_init(void * state, unsigned int dimension)
{
  return GSL_SUCCESS;
}


static int nied2_get(void * state, unsigned int dimension, double * v)
{
  return GSL_SUCCESS;
}


static const gsl_qrng_type nied2_type = 
{
  "niederreiter-base-2",
  NIED2_MAX_DIMENSION,
  nied2_state_size,
  nied2_init,
  nied2_get
};


const gsl_qrng_type * gsl_qrng_niederreiter_2 = &nied2_type;
