#include <config.h>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_histogram.h>

int
gsl_histogram_fread (FILE * stream, gsl_histogram * h)
{
  int status = gsl_block_raw_fread (stream, h->range, h->n + 1, 1);

  if (status)
    return status;

  status = gsl_block_raw_fread (stream, h->bin, h->n, 1);
  return status;
}

int
gsl_histogram_fwrite (FILE * stream, const gsl_histogram * h)
{
  int status = gsl_block_raw_fwrite (stream, h->range, h->n + 1, 1);

  if (status)
    return status;

  status = gsl_block_raw_fwrite (stream, h->bin, h->n, 1);
  return status;
}

int
gsl_histogram_fprintf (FILE * stream, const gsl_histogram * h,
		       const char *range_format, const char *bin_format)
{
  size_t i;
  const size_t n = h->n;

  for (i = 0; i < n; i++)
    {
      int status = fprintf (stream, range_format, h->range[i]);

      if (status < 0)
	{
	  GSL_ERROR ("fprintf failed", GSL_EFAILED);
	}

      status = putc (' ', stream);

      if (status == EOF)
	{
	  GSL_ERROR ("putc failed", GSL_EFAILED);
	}

      status = fprintf (stream, range_format, h->range[i + 1]);

      if (status < 0)
	{
	  GSL_ERROR ("fprintf failed", GSL_EFAILED);
	}

      status = putc (' ', stream);

      if (status == EOF)
	{
	  GSL_ERROR ("putc failed", GSL_EFAILED);
	}

      status = fprintf (stream, bin_format, h->bin[i]);

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
gsl_histogram_fscanf (FILE * stream, gsl_histogram * h)
{
  size_t i;
  const size_t n = h->n;
  double upper;

  for (i = 0; i < n; i++)
    {
      int status = fscanf (stream,
			   "%lg %lg %lg", h->range + i, &upper,
			   h->bin + i);

      if (status != 3)
	{
	  GSL_ERROR ("fscanf failed", GSL_EFAILED);
	}
    }

  h->range[n] = upper;

  return 0;
}
