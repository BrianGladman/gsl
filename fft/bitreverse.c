#include <gsl_fft.h>

int gsl_fft_complex_bitreverse_order (complex data[], 
				      const unsigned int n,
				      const unsigned int logn)
{
  unsigned int i;

  for (i = 0; i < n; i++)
    {
      unsigned int j = 0;
      unsigned int i_tmp = i;
      unsigned int bit;

      for (bit = 0; bit < logn; bit++)
	{
	  j <<= 1;		/* reverse shift i into j */
	  j |= i_tmp & 1;
	  i_tmp >>= 1;
	}
      
      if (i < j)
	{
	  const complex data_tmp = data[i];
	  data[i] = data[j];
	  data[j] = data_tmp;
	}
    }
  return 0;
}

int gsl_fft_real_bitreverse_order (double data[], 
				   const unsigned int n,
				   const unsigned int logn)
{
  unsigned int i;

  for (i = 0; i < n; i++)
    {
      unsigned int j = 0;
      unsigned int i_tmp = i;
      unsigned int bit;

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






