#if !(defined(USES_LONGDOUBLE) && !defined(HAVE_PRINTF_LONGDOUBLE))

int
FUNCTION (gsl_block, fprintf) (FILE * stream, const TYPE(gsl_block) * b, const char *format)
{
  size_t n = b->size ;
  
  ATOMIC * data = b->data ;
  
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      int status;

      for (k = 0; k < MULTIPLICITY; k++)
	{
	  if (k > 0)
	    {
	      status = putc (' ', stream);

	      if (status == EOF)
		{
		  GSL_ERROR ("putc failed", GSL_EFAILED);
		}
	    }
	  status = fprintf (stream,
			    format,
			    data[MULTIPLICITY * i + k]);
	  if (status < 0)
	    {
	      GSL_ERROR ("fprintf failed", GSL_EFAILED);
	    }
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
FUNCTION (gsl_block, fscanf) (FILE * stream, TYPE(gsl_block) * b)
{
  size_t n = b->size ;
  
  ATOMIC * data = b->data ;

  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      for (k = 0; k < MULTIPLICITY; k++)
	{
	  int status = fscanf (stream,
			       IN_FORMAT,
			       data + MULTIPLICITY * i + k);
	  if (status != 1)
            {
              GSL_ERROR ("fscanf failed", GSL_EFAILED);
            }
	}
    }

  return GSL_SUCCESS;
}


int
FUNCTION (gsl_block, raw_fprintf) (FILE * stream, 
                                   const ATOMIC * data,
                                   const size_t n,
                                   const size_t stride, 
                                   const char *format)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      int status;

      for (k = 0; k < MULTIPLICITY; k++)
	{
	  if (k > 0)
	    {
	      status = putc (' ', stream);

	      if (status == EOF)
		{
		  GSL_ERROR ("putc failed", GSL_EFAILED);
		}
	    }
	  status = fprintf (stream,
			    format,
			    data[MULTIPLICITY * i * stride + k]);
	  if (status < 0)
	    {
	      GSL_ERROR ("fprintf failed", GSL_EFAILED);
	    }
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
FUNCTION (gsl_block, raw_fscanf) (FILE * stream, 
                                  ATOMIC * data,
                                  const size_t n, 
                                  const size_t stride)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      for (k = 0; k < MULTIPLICITY; k++)
	{
	  int status = fscanf (stream,
			       IN_FORMAT,
			       data + MULTIPLICITY * i * stride + k);
	  if (status != 1)
	    GSL_ERROR ("fscanf failed", GSL_EFAILED);
	}
    }

  return GSL_SUCCESS;
}

#endif
