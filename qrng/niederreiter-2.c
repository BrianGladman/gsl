/* Author: G. Jungman
 * RCS:    $Id$
 */
#include "gsl_qrng.h"


#define NIED2_MAX_DIMENSION 12
#define NIED2_BIT_COUNT 30
#define NIED2_NBITS (NIED2_BIT_COUNT+1)


static size_t nied2_state_size(unsigned int dimension);
static int nied2_init(void * state, unsigned int dimension);
static int nied2_get(void * state, unsigned int dimension, double * v);


static const gsl_qrng_type nied2_type = 
{
  "niederreiter-base-2",
  NIED2_MAX_DIMENSION,
  nied2_state_size,
  nied2_init,
  nied2_get
};

const gsl_qrng_type * gsl_qrng_niederreiter_2 = &nied2_type;


typedef struct
{
  unsigned int sequence_count;
  int cj[NIED2_BIT_COUNT][NIED2_MAX_DIMENSION];
  int nextq[NIED2_MAX_DIMENSION];
} nied2_state_t;


static size_t nied2_state_size(unsigned int dimension)
{
  return sizeof(nied2_state_t);
}


static void calculate_cj(nied2_state_t * ns, unsigned int dimension)
{
}


static int nied2_init(void * state, unsigned int dimension)
{
  nied2_state_t * n_state = (nied2_state_t *) state;
  int i_dim;

  if(dimension < 1 || dimension > NIED2_MAX_DIMENSION) return GSL_EINVAL;

  calculate_cj(n_state, dimension);

  for(i_dim=0; i_dim<dimension; i_dim++) n_state->nextq[i_dim] = 0;
  n_state->sequence_count = 0;

  return GSL_SUCCESS;
}


static int nied2_get(void * state, unsigned int dimension, double * v)
{
  static const double recip = 1.0/(double)(1 << NIED2_NBITS); /* 2^(-nbits) */
  nied2_state_t * n_state = (nied2_state_t *) state;
  int r;
  int c;

  /* Load the result from the saved state. */
  for(i_dim=0; i_dim<dimension; i_dim++) {
    v[i_dim] = n_state->nextq[i_dim] * recip;
  }

  /* Find the position of the least-significant zero in sequence_count.
   * This is the bit that changes in the Gray-code representation as
   * the count is advanced.
   */
  r = 0;
  c = n_state->sequence_count;
  while(1) {
    if((c % 2) == 1) {
      ++r;
      c /= 2;
    }
    else break;
  }

  if(r >= NIED2_NBITS) return GSL_EFAILED; /* FIXME: better error code here */

  /* Calculate the next state. */
  for(i_dim=0; i_dim<dimension; i_dim++) {
    n_state->nextq[i_dim] ^= n_state->cj[r][i_dim];
  }

  n_state->sequence_count++;

  return GSL_SUCCESS;
}
