

#include <stdio.h>

int
  FUNCTION (gsl_block, fprintf) (FILE * stream, const ATOMIC * data, const size_t n,
				 const size_t stride, const char *format)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
#if MULTIPLICITY == 1
      int status = fprintf (stream,
			    format,
			    data[i * stride]);
#elif MULTIPLICITY == 2
      int status = fprintf (stream,
			    format,
			    data[2 * i * stride],
			    data[2 * i * stride + 1]);
#else
#error Unsupported multiplicity > 2
#endif

      if (status < 0)
	{
	  GSL_ERROR ("fprintf failed", GSL_EFAILED);
	}

      status = putc ('\n', stream);
      if (status == EOF)
	{
	  GSL_ERROR ("putc failed", GSL_EFAILED);
	}
    }

  return 0;
}

int
  FUNCTION (gsl_block, fscanf) (FILE * stream, ATOMIC * data, const size_t n, const size_t stride)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
#if MULTIPLICITY == 1
      int status = fscanf (stream, IN_FORMAT, data + i * stride);
#elif MULTIPLICITY == 2
      int status = fscanf (stream,
			   IN_FORMAT,
			   data + 2 * i * stride,
			   data + 2 * i * stride + 1);
#else
#error Unsupported multiplicity > 2
#endif

      if (status != MULTIPLICITY)
	GSL_ERROR ("fscanf failed", GSL_EFAILED);
    }

  return GSL_SUCCESS;
}
