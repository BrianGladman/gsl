/* Authors: 
   G. Jungman
   J. Scott (james.scott@lexifi.com)
 */
/* Implementation for Sobol generator.
 * See
 *   [Bratley+Fox, TOMS 14, 88 (1988)]
 *   [Antonov+Saleev, USSR Comput. Maths. Math. Phys. 19, 252 (1980)]
 */
#include <config.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_rng.h>

#define SOBOL_BIT_COUNT 30

/* prototypes for generator type functions */
static size_t sobol_state_size(unsigned int dimension);
static int sobol_init(void * state, unsigned int dimension);
static int sobol_get(void * state, unsigned int dimension, double * v);

/* global Sobol generator type object */
static const gsl_qrng_type sobol_type = 
{
  "sobol",
  0,
  sobol_state_size,
  sobol_init,
  sobol_get
};
const gsl_qrng_type * gsl_qrng_sobol = &sobol_type;


/* Sobol generator state.
 *   sequence_count       = number of calls with this generator
 *   last_numerator_vec   = last generated numerator vector
 *   last_denominator_inv = 1/denominator for last numerator vector
 *   v_direction          = direction number table
 */
typedef struct
{
  unsigned int  sequence_count;
  double        last_denominator_inv;
  int           *last_numerator_vec;
  int           *v_direction[SOBOL_BIT_COUNT];
} sobol_state_t;


static size_t sobol_state_size(unsigned int dimension)
{
  return 
    sizeof(sobol_state_t) +      /* The struct */ 
    sizeof(int) * dimension +    /* for last_numerator_vec */ 
    sizeof(int) * dimension * SOBOL_BIT_COUNT; /* for the direction no.s */ 
}

/* l is a degree number, between 1 and the degree of the polynomial
   associated with the dimension being initialized */ 
static int sobol_init_direction(gsl_rng *r, int l)
{
  /* See Peter Jaeckel, Monte Carlo Methods in Finance, Wiley 2002, p86. 
   */ 
  int wkl = gsl_rng_uniform_int(r, 1<<l) | 1; 
  return wkl; 
}

void get_primitive_polynomials(int dimension, int *degree_table, int *primitive_polynomials); 

static int sobol_init(void * state, unsigned int dimension)
{
  sobol_state_t * s_state = (sobol_state_t *) state;
  unsigned int i_dim;
  int j, k;
  int ell;
  int *degree_table; 
  int *primitive_polynomials; 
  int *includ; 
  int max_degree; 
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default); 

  if(dimension < 1) {
    return GSL_EINVAL;
  }

  { /* initialize degree_table, primitive_polynomials, max_degree and includ */     
    degree_table = (int *) malloc(sizeof(int) * dimension); 
    if (degree_table==0) 
      GSL_ERROR_NULL ("allocation of degree table failed for sobol init", GSL_ENOMEM); 
    primitive_polynomials = (int *) malloc(sizeof(int) * dimension); 
    if (primitive_polynomials==0) {
      free(degree_table); 
      GSL_ERROR_NULL ("allocation of primitives failed for sobol init", GSL_ENOMEM); 
    }

    /* Generate the primitive polynomials */ 
    get_primitive_polynomials(dimension, degree_table, primitive_polynomials); 

    max_degree = degree_table[dimension-1]; 
    includ = (int *) malloc(sizeof(int) * max_degree); 
    if (includ==0) {
      free(degree_table); 
      free(primitive_polynomials); 
      GSL_ERROR_NULL ("allocation of 'includ' failed for sobol init", GSL_ENOMEM); 
    }
  }


  s_state->last_numerator_vec = (int *) ((char *) s_state + sizeof(sobol_state_t)); 

  /* Initialize direction table in dimension 0. */
  for(k=0; k<SOBOL_BIT_COUNT; k++) {
    s_state->v_direction[k] = 
      (int *) (((char *) s_state) + 
	       sizeof(sobol_state_t) + 
	       sizeof(int) * dimension + 
	       sizeof(int) * dimension * k); 
    s_state->v_direction[k][0] = 1;
  }

  /* Initialize in remaining dimensions. */
  for(i_dim=1; i_dim<dimension; i_dim++) {

    const int poly_index = i_dim;
    const int degree_i = degree_table[poly_index];

    /* Expand the polynomial bit pattern to separate
     * components of the logical array includ[].
     */
    int p_i = primitive_polynomials[poly_index];
    for(k = degree_i-1; k >= 0; k--) {
      includ[k] = ((p_i % 2) == 1);
      p_i /= 2;
    }

    /* Leading elements for dimension i are randomly initialized */
    for(j=0; j<degree_i; j++) s_state->v_direction[j][i_dim] = sobol_init_direction(r, j+1); 

    /* Calculate remaining elements for this dimension,
     * as explained in Bratley+Fox, section 2.
     */
    for(j=degree_i; j<SOBOL_BIT_COUNT; j++) {
      int newv = s_state->v_direction[j-degree_i][i_dim];
      ell = 1;
      for(k=0; k<degree_i; k++) {
        ell *= 2;
        if(includ[k]) newv ^= (ell * s_state->v_direction[j-k-1][i_dim]);
      }
      s_state->v_direction[j][i_dim] = newv;
    }
  }

  /* Multiply columns of v by appropriate power of 2. */
  ell = 1;
  for(j=SOBOL_BIT_COUNT-1-1; j>=0; j--) {
    ell *= 2;
    for(i_dim=0; i_dim<dimension; i_dim++) {
      s_state->v_direction[j][i_dim] *= ell;
    }
  }

  /* 1/(common denominator of the elements in v_direction) */
  s_state->last_denominator_inv = 1.0 /(2.0 * ell);

  /* final setup */
  s_state->sequence_count = 0;
  for(i_dim=0; i_dim<dimension; i_dim++) s_state->last_numerator_vec[i_dim] = 0;

  free(degree_table); 
  free(primitive_polynomials); 
  free(includ); 
  gsl_rng_free(r); 

  return GSL_SUCCESS;
}


static int sobol_get(void * state, unsigned int dimension, double * v)
{
  sobol_state_t * s_state = (sobol_state_t *) state;

  unsigned int i_dimension;

  /* Find the position of the least-significant zero in count. */
  int ell = 0;
  int c = s_state->sequence_count;
  while(1) {
    ++ell;
    if((c % 2) == 1) c /= 2;
    else break;
  }

  /* Check for exhaustion. */
  if(ell > SOBOL_BIT_COUNT) return GSL_EFAILED; /* FIXME: good return code here */

  for(i_dimension=0; i_dimension<dimension; i_dimension++) {
    const int direction_i     = s_state->v_direction[ell-1][i_dimension];
    const int old_numerator_i = s_state->last_numerator_vec[i_dimension];
    const int new_numerator_i = old_numerator_i ^ direction_i;
    s_state->last_numerator_vec[i_dimension] = new_numerator_i;
    v[i_dimension] = new_numerator_i * s_state->last_denominator_inv;
  }

  s_state->sequence_count++;

  return GSL_SUCCESS;
}
