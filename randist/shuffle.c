#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

/* Inline swap and copy functions for moving objects around */

static inline 
void swap (void * base, size_t size, size_t i, size_t j)
{
  register char * a = size * i + (char *) base ;
  register char * b = size * j + (char *) base ;
  register size_t s = size ;

  if (i == j)
    return ;
  
  do						
    {						
      char tmp = *a;				
      *a++ = *b;				
      *b++ = tmp;				
    } 
  while (--s > 0);				
}

static inline void 
copy (void * dest, size_t i, void * src, size_t j, size_t size)
{
  register char * a = size * i + (char *) dest ;
  register char * b = size * j + (char *) src ;
  register size_t s = size ;
  
  do						
    {						
      *a++ = *b++;				
    } 
  while (--s > 0);				
}

/* Randomly permute (shuffle) N indices

   Supply an array x[N] with nmemb members, each of size size and on
   return it will be shuffled into a random order.  The algorithm is
   from Knuth, SemiNumerical Algorithms, v2, p139, who cites Moses and
   Oakford, and Durstenfeld */

void
gsl_ran_shuffle (const gsl_rng * r, void * base, size_t n, size_t size)
{
  size_t i ;

  for (i = n - 1; i > 0; i--)
    {
      size_t j = (i + 1) * gsl_rng_uniform (r); 

      swap (base, size, i, j) ;
    }
}

void *
gsl_ran_choose (const gsl_rng * r, void * dest, size_t k, void * src, 
		 size_t n, size_t size)
{
  size_t i, j = 0;

  /* Choose k out of n items, return an array x[] of the k items.
     These items will prevserve the relative order of the original
     input -- you can use shuffle() to randomize the output if you
     wish */

  if (k > n)
    {
      GSL_ERROR_RETURN ("k is greater than n, cannot sample more than n items",
			GSL_EINVAL, 0) ;
    }

  for (i = 0; i < n && j < k; i++)
    {
      if ((n - i) * gsl_rng_uniform (r) < k - j)
	{
	  copy (dest, j, src, i, size) ;
	  j++ ;
	}
    }

  return dest;
}

void *
gsl_ran_sample (const gsl_rng * r, void * dest, size_t k, void * src, 
		size_t n, size_t size)
{
  size_t i, j = 0;

  /* Choose k out of n items, with replacement */

  if (k > n)
    {
      GSL_ERROR_RETURN ("k is greater than n, cannot sample more than n items",
			GSL_EINVAL, 0) ;
    }

  for (i = 0; i < k; i++)
    {
      j = n * gsl_rng_uniform (r) ;

      copy (dest, i, src, j, size) ;
    }

  return dest;
}
