#include <config.h>
#include <gsl_fft.h>

#include "bitreverse.h"

int 
bitreverse_complex (double data[], 
		    const size_t stride,
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
	  const double tmp_real = REAL(data,stride,i);
	  const double tmp_imag = IMAG(data,stride,i);
	  REAL(data,stride,i) = REAL(data,stride,j);
	  IMAG(data,stride,i) = IMAG(data,stride,j);
	  REAL(data,stride,j) = tmp_real;
	  IMAG(data,stride,j) = tmp_imag;
	}
    }
  return 0;
}

int bitreverse_real (double data[], 
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


int 
goldrader_bitreverse_order_complex (gsl_complex data[], 
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


int 
rodriguez_bitreverse_order_complex (gsl_complex data[], 
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







