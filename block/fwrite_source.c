int
FUNCTION (gsl_block, fread) (FILE * stream, TYPE(gsl_block) * b)
{
  size_t n = b->size ;

  ATOMIC * data = b->data ;

  size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fread failed", GSL_EFAILED);
    }
      
return 0;
}

int
FUNCTION (gsl_block, fwrite) (FILE * stream, const TYPE(gsl_block) * b)
{
  size_t n = b->size ;

  ATOMIC * data = b->data ;
  
  size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fwrite failed", GSL_EFAILED);
    }

  return 0;
}

int
FUNCTION (gsl_block, raw_fread) (FILE * stream, ATOMIC * data, 
                                 const size_t n, const size_t stride)
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

  return GSL_SUCCESS;
}

int
FUNCTION (gsl_block, raw_fwrite) (FILE * stream, const ATOMIC * data,
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

  return GSL_SUCCESS;
}
