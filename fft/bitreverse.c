#include <config.h>
#include <gsl_fft.h>

#include "bitreverse.h"

int gsl_fft_complex_bitreverse_order (gsl_complex data[], 
				      const size_t n,
				      const size_t logn)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      size_t j = 0;
      size_t i_tmp = i;
      size_t bit;

      for (bit = 0; bit < logn; bit++)
	{
	  j <<= 1;		/* reverse shift i into j */
	  j |= i_tmp & 1;
	  i_tmp >>= 1;
	}
      
      if (i < j)
	{
	  const gsl_complex data_tmp = data[i];
	  data[i] = data[j];
	  data[j] = data_tmp;
	}
    }
  return 0;
}

int gsl_fft_real_bitreverse_order (double data[], 
				   const size_t n,
				   const size_t logn)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      size_t j = 0;
      size_t i_tmp = i;
      size_t bit;

      for (bit = 0; bit < logn; bit++)
	{
	  j <<= 1;		/* reverse shift i into j */
	  j |= i_tmp & 1;
	  i_tmp >>= 1;
	}
      
      if (i < j)
	{
	  const double data_tmp = data[i];
	  data[i] = data[j];
	  data[j] = data_tmp;
	}
    }
  return 0;
}


int gsl_fft_complex_goldrader_bitreverse_order (gsl_complex data[], 
						const size_t n)
{
  size_t i;
  size_t j = 0;

  for (i = 0; i < n - 1; i++)
    {
      size_t k = n / 2 ;

      if (i < j)
	{
	  const gsl_complex data_tmp = data[i];
	  data[i] = data[j];
	  data[j] = data_tmp;
	}

      while (k <= j) 
	{
	  j = j - k ;
	  k = k / 2 ;
	}
      
      j += k ;

    }

  return 0;
}


int gsl_fft_complex_rodriguez_bitreverse_order (gsl_complex data[], 
						const size_t n,
						const size_t logn)
{
  size_t i;
  size_t j = 0 ;
  unsigned last = n - ( 1 << ((logn + 1)/2) ) ;

  for (i = 1; i < last ; i++)
    {
      size_t k = n / 2 ;

      while (k <= j) 
	{
	  j = j - k ;
	  k = k / 2 ;
	}

      j += k ;

      if (i < j)
	{
	  const gsl_complex data_tmp = data[i];
	  data[i] = data[j];
	  data[j] = data_tmp;
	}

    }

  return 0;
}







