static inline void FUNCTION(my,downheap) (BASE *data, const size_t stride, const size_t N, size_t k);

static inline void
FUNCTION(my,downheap) (BASE *data, const size_t stride, const size_t N, size_t k)
{
  BASE v = data[k*stride] ;

  while (k <= N / 2)
    {
      size_t j = 2 * k;

      if (j < N && data[j*stride] < data[(j+1)*stride])
        {
          j++;
        }
      
      if (v >= data[j*stride])
        {
          break ;
        }

      data[k*stride] = data[j*stride] ;

      k = j;
    }

  data[k*stride] = v;
}

void
TYPE (gsl_sort_vector) (TYPE (gsl_vector) * v)
{
  BASE * data = v->data ;
  const size_t n = v->size;
  const size_t stride = v->stride ;
  
  size_t N;
  size_t k;

  if (n == 0)
    {
      return ; /* No data to sort */
    }

  /* We have n_data elements, last element is at 'n_data-1', first at
     '0' Set N to the last element number. */

  N = n - 1;

  k = N / 2;
  k++;                          /* Compensate the first use of 'k--' */
  do
    {
      k--;
      FUNCTION(my,downheap) (data, stride, N, k);
    }
  while (k > 0);

  while (N > 0)
    {
      /* first swap the elements */
      BASE tmp = data[0*stride] ;
      data[0*stride] = data[N*stride] ;
      data[N*stride] = tmp ;

      /* then process the heap */
      N--;

      FUNCTION(my,downheap) (data, stride, N, 0);
    }
}





