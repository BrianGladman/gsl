/* 
   Authors: 
   G. Jungman
   J. Scott (james.scott@lexifi.com)
   
 */
/* Implement Niederreiter base 2 generator.
 * See:
 *   Bratley, Fox, Niederreiter, ACM Trans. Model. Comp. Sim. 2, 195 (1992)
 */
#include <config.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_rng.h>
#include <assert.h>

#define NIED2_BIT_COUNT 30
#define NIED2_NBITS (NIED2_BIT_COUNT+1)

/* Z_2 field operations */
#define NIED2_ADD(x,y) (((x)+(y))%2)
#define NIED2_MUL(x,y) (((x)*(y))%2)
#define NIED2_SUB(x,y) NIED2_ADD((x),(y))

void get_primitive_polynomials(int dimension, int *degree_table, int *primitive_polynomials);

static size_t nied2_state_size(unsigned int dimension);
static int nied2_init(void * state, unsigned int dimension);
static int nied2_get(void * state, unsigned int dimension, double * v);


static const gsl_qrng_type nied2_type = 
{
  "niederreiter-base-2",
  0,
  nied2_state_size,
  nied2_init,
  nied2_get
};

const gsl_qrng_type * gsl_qrng_niederreiter_2 = &nied2_type;


typedef struct
{
  unsigned int sequence_count;
  int *cj[NIED2_NBITS];
  int *nextq;
} nied2_state_t;


static size_t nied2_state_size(unsigned int dimension)
{
  return sizeof(nied2_state_t) +             /* the struct */ 
    sizeof(int) * dimension * NIED2_NBITS +  /* cj */ 
    sizeof(int) * dimension;                 /* nextq */ 
}


/* Multiply polynomials over Z_2.
 * Notice use of a temporary vector,
 * side-stepping aliasing issues when
 * one of inputs is the same as the output
 * [especially important in the original fortran version, I guess].
 */
static void poly_multiply(
  const int pa[], int pa_degree,
  const int pb[], int pb_degree,
  int pc[], int  * pc_degree, int max_degree
  )
{
  int j, k;
  int *pt; 
  int pt_degree; 

  pt = (int *) malloc(sizeof(int) * (max_degree+1)); 
  pt_degree = pa_degree + pb_degree;

  assert(pt_degree <= max_degree); 

  for(k=0; k<=pt_degree; k++) {
    int term = 0;
    for(j=0; j<=k; j++) {
      const int conv_term = NIED2_MUL(pa[k-j], pb[j]);
      term = NIED2_ADD(term, conv_term);
    }
    pt[k] = term;
  }

  for(k=0; k<=pt_degree; k++) {
    pc[k] = pt[k];
  }
  for(k=pt_degree+1; k<=max_degree; k++) {
    pc[k] = 0;
  }

  *pc_degree = pt_degree;
  free(pt); 
}

/* Calculate the values of the constants V(J,R) as
 * described in BFN section 3.3.
 *
 *   px = appropriate irreducible polynomial for current dimension
 *   pb = polynomial defined in section 2.3 of BFN.
 * pb is modified
 */
static void calculate_v(
  const  int px[],  int px_degree,
   int pb[],  int * pb_degree,
  int v[],  int maxv, int max_degree, gsl_rng *rng
  )
{
  const int nonzero_element = 1;    /* nonzero element of Z_2  */

  /* The polynomial ph is px**(J-1), which is the value of B on arrival.
   * In section 3.3, the values of Hi are defined with a minus sign:
   * don't forget this if you use them later !
   */
  int *ph;
  /* int ph_degree = *pb_degree; */
  int bigm = *pb_degree;      /* m from section 3.3 */
  int m;                      /* m from section 2.3 */
  int r, k, kj;

  ph = (int *) malloc(sizeof(int) * (max_degree+1)); 

  for(k=0; k<=max_degree; k++) {
    ph[k] = pb[k];
  }

  /* Now multiply B by PX so B becomes PX**J.
   * In section 2.3, the values of Bi are defined with a minus sign :
   * don't forget this if you use them later !
   */
   poly_multiply(px, px_degree, pb, *pb_degree, pb, pb_degree, max_degree);
   m = *pb_degree;

  /* Now choose a value of Kj as defined in section 3.3.
   * We must have 0 <= Kj < E*J = M.
   * The limit condition on Kj does not seem very relevant
   * in this program.
   */
  /* Quoting from BFN: "Our program currently sets each K_q
   * equal to eq. This has the effect of setting all unrestricted
   * values of v to 1."
   * Actually, it sets them to the arbitrary chosen value.
   * Whatever.
   */
   kj = gsl_rng_uniform_int(rng, bigm+1);    

  /* Now choose values of V in accordance with
   * the conditions in section 3.3.
   */
  for(r=0; r<kj; r++) {
    v[r] = 0;
  }
  v[kj] = 1;


  if(kj >= bigm) {
    for(r=kj+1; r<m; r++) {
      v[r] =  gsl_rng_uniform_int(rng, 2); /* arbitrary element */
    }
  }
  else {
    /* This block is never reached. */

    int term = NIED2_SUB(0, ph[kj]);

    for(r=kj+1; r<bigm; r++) {
      v[r] =  gsl_rng_uniform_int(rng, 2); /* arbitrary element */

      /* Check the condition of section 3.3,
       * remembering that the H's have the opposite sign.  [????????]
       */
      term = NIED2_SUB(term, NIED2_MUL(ph[r], v[r]));
    }

    /* Now v[bigm] != term. */
    v[bigm] = NIED2_ADD(nonzero_element, term);

    for(r=bigm+1; r<m; r++) {
      v[r] = gsl_rng_uniform_int(rng, 2); /* arbitrary element */ 
    }
  }

  /* Calculate the remaining V's using the recursion of section 2.3,
   * remembering that the B's have the opposite sign.
   */
  for(r=0; r<=maxv-m; r++) {
    int term = 0;
    for(k=0; k<m; k++) {
      term = NIED2_SUB(term, NIED2_MUL(pb[k], v[r+k]));
    }
    v[r+m] = term;
  }
  free(ph); 
}


static int calculate_cj(nied2_state_t * ns,  int dimension)
{
  int ci[NIED2_NBITS][NIED2_NBITS];
  int *v;
  int r;
  int i_dim;
  int max_degree, maxv; 
  int *pb, *px; 

  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default); 

  int *primitive_polynomials; 
  int *poly_degree; 

  poly_degree = (int *) malloc(sizeof(int) * (1+dimension)); 
  primitive_polynomials = (int *) malloc(sizeof(int) * (1+dimension)); 

  if (poly_degree==0 || primitive_polynomials==0) 
    GSL_ERROR_NULL ("allocation of degree table failed for niederreiter init", GSL_ENOMEM); 
  
  /* Generate the primitive polynomials.  */ 
  get_primitive_polynomials(dimension, poly_degree+1, primitive_polynomials+1); 

  poly_degree[0]=0;  poly_degree[1]=1;
  primitive_polynomials[0]=1; 
  primitive_polynomials[1]=2; 

  max_degree=NIED2_NBITS+poly_degree[dimension-1]; 
  maxv = max_degree; 

  v = (int *) malloc(sizeof(int) * (1+maxv));   
  pb = (int *) malloc(sizeof(int) * (1+max_degree)); 
  px = (int *) malloc(sizeof(int) * (1+max_degree)); 

  if (v==0 || pb==0 || px==0) 
    GSL_ERROR_NULL ("allocation of degree table failed for niederreiter init", GSL_ENOMEM); 

  for(i_dim=0; i_dim<dimension; i_dim++) {

    const int poly_index = i_dim + 1;
    int j, k;

    /* Niederreiter (page 56, after equation (7), defines two
     * variables Q and U.  We do not need Q explicitly, but we
     * do need U.
     */
    int u = 0;
     
    /* For each dimension, we need to calculate powers of an
     * appropriate irreducible polynomial, see Niederreiter
     * page 65, just below equation (19).
     * Copy the appropriate irreducible polynomial into PX,
     * and its degree into E.  Set polynomial B = PX ** 0 = 1.
     * M is the degree of B.  Subsequently B will hold higher
     * powers of PX.
     */
    int px_degree; 
    int pb_degree; 
    px_degree = poly_degree[poly_index];
    pb_degree = 0;
    
    for(k=0; k<=px_degree; k++) {
      px[k] = (primitive_polynomials[poly_index] >> k) & 1;
      pb[k] = 0;
    }

    for (;k<max_degree+1;k++) {
      px[k] = 0;
      pb[k] = 0;
    }

    pb[0] = 1;

    for(j=0; j<NIED2_NBITS; j++) {
      /* If U = 0, we need to set B to the next power of PX
       * and recalculate V.
       */
      if(u == 0) calculate_v(px, px_degree, pb, &pb_degree, v, maxv, max_degree, rng);

      /* Now C is obtained from V.  Niederreiter
       * obtains A from V (page 65, near the bottom), and then gets
       * C from A (page 56, equation (7)).  However this can be done
       * in one step.  Here CI(J,R) corresponds to
       * Niederreiter's C(I,J,R).
       */
      for(r=0; r<NIED2_NBITS; r++) {
        ci[r][j] = v[r+u];
      }

      /* Advance Niederreiter's state variables. */
      ++u;
      if(u == px_degree) u = 0;
    }

    /* The array CI now holds the values of C(I,J,R) for this value
     * of I.  We pack them into array CJ so that CJ(I,R) holds all
     * the values of C(I,J,R) for J from 1 to NBITS.
     */
    for(r=0; r<NIED2_NBITS; r++) {
      int term = 0;
      for(j=0; j<NIED2_NBITS; j++) {
        term = 2*term + ci[r][j];
      }
      ns->cj[r][i_dim] = term;
    }
  }
  free(primitive_polynomials);
  free(poly_degree); 
  free(v); 
  free(pb); 
  free(px); 
  gsl_rng_free(rng); 
  return GSL_SUCCESS; 
}


static int nied2_init(void * state, unsigned int dimension)
{
  nied2_state_t * n_state = (nied2_state_t *) state;
  unsigned int i_dim, i_bits;
  int ret; 

  if(dimension < 1) return GSL_EINVAL;

  for (i_bits=0; i_bits<NIED2_NBITS; i_bits++)
    n_state->cj[i_bits] = 
      (int *) ((char *) n_state + 
	       sizeof(nied2_state_t) + 
	       sizeof(int) * dimension * i_bits); 

  n_state->nextq = 
    (int *) ((char *) n_state + 
	     sizeof(nied2_state_t) + 
	     sizeof(int) * dimension * NIED2_NBITS); 

  if ((ret = calculate_cj(n_state, dimension)) != GSL_SUCCESS)
    return ret; 

  /* skip the '0th' state, which is 0 */ 
  for(i_dim=0; i_dim<dimension; i_dim++) 
    n_state->nextq[i_dim] = n_state->cj[0][i_dim]; 

  n_state->sequence_count = 1;

  return GSL_SUCCESS;
}


static int nied2_get(void * state, unsigned int dimension, double * v)
{
  static const double recip = 1.0/(double)(1U << NIED2_NBITS); /* 2^(-nbits) */
  nied2_state_t * n_state = (nied2_state_t *) state;
  int r;
  int c;
  unsigned int i_dim;

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
