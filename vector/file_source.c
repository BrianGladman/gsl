#include <stdio.h>
#include <unistd.h>
#include <gsl_errno.h>

int
FUNCTION (gsl_vector, fread) (FILE * stream, TYPE (gsl_vector) * v)
{
  int status = FUNCTION (gsl_block, fread) (stream,
					    v->data,
					    v->size,
					    v->stride);
  return status;
}

int
FUNCTION (gsl_vector, fwrite) (FILE * stream, const TYPE (gsl_vector) * v)
{
  int status = FUNCTION (gsl_block, fwrite) (stream,
					     v->data,
					     v->size,
					     v->stride);
  return status;
}

#if ! (defined(BASE_LONG_DOUBLE) && !defined(HAVE_PRINTF_LONGDOUBLE))
int
FUNCTION (gsl_vector, fprintf) (FILE * stream, const TYPE (gsl_vector) * v,
				const char *format)
{
  int status = FUNCTION (gsl_block, fprintf) (stream,
					      v->data,
					      v->size,
					      v->stride,
					      format);
  return status;
}

int
FUNCTION (gsl_vector, fscanf) (FILE * stream, TYPE (gsl_vector) * v)
{
  int status = FUNCTION (gsl_block, fscanf) (stream,
					     v->data,
					     v->size,
					     v->stride);
  return status;
}
#endif

int
FUNCTION (gsl_block, fread) (FILE * stream, ATOMIC * data, const size_t n,
			     const size_t stride)
{
  if (stride == 1)
    {
      size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

      if (items != n)
	{
	  GSL_ERROR ("fread failed", GSL_EFAILED);
	}
    }
  else
    {
      size_t i;

      for (i = 0; i < n; i++)
	{
	  size_t item = fread (data + MULTIPLICITY * i * stride,
			       MULTIPLICITY * sizeof (ATOMIC), 1, stream);
	  if (item != 1)
	    {
	      GSL_ERROR ("fread failed", GSL_EFAILED);
	    }
	}
    }

  return 0;
}

int
FUNCTION (gsl_block, fwrite) (FILE * stream, const ATOMIC * data,
			      const size_t n, const size_t stride)
{

  if (stride == 1)
    {
      size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

      if (items != n)
	{
	  GSL_ERROR ("fwrite failed", GSL_EFAILED);
	}
    }
  else
    {
      size_t i;

      for (i = 0; i < n; i++)
	{
	  size_t item = fwrite (data + MULTIPLICITY * i * stride,
				MULTIPLICITY * sizeof (ATOMIC),
				1, stream);
	  if (item != 1)
	    {
	      GSL_ERROR ("fwrite failed", GSL_EFAILED);
	    }
	}
    }

  return 0;
}
