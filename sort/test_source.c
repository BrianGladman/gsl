void TYPE(test_sort_vector)(size_t N, size_t stride);
void FUNCTION(my,initialize) (TYPE(gsl_vector) * v);
void FUNCTION(my,randomize) (TYPE(gsl_vector) * v);
int FUNCTION(my,check) (TYPE(gsl_vector) * data, TYPE(gsl_vector) * orig);

void 
TYPE(test_sort_vector)(size_t N, size_t stride) 
{
  int status;

  TYPE(gsl_block) * b1 = FUNCTION(gsl_block, calloc) (N*stride) ;
  TYPE(gsl_block) * b2 = FUNCTION(gsl_block, calloc) (N*stride) ;

  TYPE(gsl_vector) * orig = FUNCTION(gsl_vector, alloc_from_block) (b1,0,N,stride) ;
  TYPE(gsl_vector) * data = FUNCTION(gsl_vector, alloc_from_block) (b2,0,N,stride) ;
  
  FUNCTION(my,initialize) (orig);

  /* Already sorted */
  FUNCTION(gsl_vector,cpy) (data, orig);
  TYPE(gsl_sort_vector) (data);
  status = FUNCTION(my,check) (data, orig);

  gsl_test (status, "sorting " NAME(gsl_vector) ", ordered, n = %u, stride = %u", N, stride) ;  

  /* Reverse the data */

  FUNCTION(gsl_vector,cpy) (data, orig);
  FUNCTION(gsl_vector,reverse) (data) ;
  TYPE(gsl_sort_vector) (data);
  status = FUNCTION(my,check) (data, orig);

  gsl_test (status, "sorting " NAME(gsl_vector) ", reversed, n = %u, stride = %u", N, stride) ;

  /* Perform some shuffling */

  FUNCTION(gsl_vector,cpy) (data, orig);
  FUNCTION(my,randomize) (data) ;
  TYPE(gsl_sort_vector) (data);
  status = FUNCTION(my,check) (data, orig);

  gsl_test (status, "sorting " NAME(gsl_vector) ", randomized, n = %u, stride = %u", N, stride) ;

  FUNCTION(gsl_vector,free) (orig);
  FUNCTION(gsl_vector,free) (data);
  FUNCTION(gsl_block,free) (b1);
  FUNCTION(gsl_block,free) (b2);
}


void
FUNCTION(my,initialize) (TYPE(gsl_vector) * v)
{
  size_t i;

  for (i = 0; i < v->size; i++) 
    {
      FUNCTION(gsl_vector,set)(v,i,(ATOMIC)i) ;
    }
}

void
FUNCTION(my,randomize) (TYPE(gsl_vector) * v)
{
  size_t i;
  
  for (i = 0; i < v->size; i++) 
    {
      size_t j = urand() * v->size;
      FUNCTION(gsl_vector,swap)(v,i,j);
    }
}

int
FUNCTION(my,check) (TYPE(gsl_vector) * data, TYPE(gsl_vector) * orig)
{
  size_t i;

  for (i = 0; i < data->size; i++) 
    {
      if (FUNCTION(gsl_vector,get)(data,i) != FUNCTION(gsl_vector,get)(orig,i)) 
        {
          return GSL_FAILURE;
        }
    }

  return GSL_SUCCESS;
}
