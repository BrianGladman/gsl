int
FUNCTION(gsl_block,fprintf)(FILE * stream, const BASE * data, const size_t n,
			     const char * format)
{
  size_t i ;

  for (i = 0 ; i < n ; i++) 
    {
      int status = fprintf(stream, format, data[i]) ;
      
      if (status < 0) 
	{
	  GSL_ERROR ("fprintf failed", GSL_EFAILED) ;
	}

      status = putc('\n', stream) ;

      if (status == EOF) 
	{
	  GSL_ERROR ("putc failed", GSL_EFAILED) ;
	}
    }

  return 0 ;
}

int
FUNCTION(gsl_block,fscanf)(FILE * stream, BASE * data, const size_t n)
{
  size_t i ;

  for (i = 0 ; i < n ; i++) 
    {
      int status = fscanf(stream, IN_FORMAT, data + i) ;

      if (status != 1) 
	{
	  GSL_ERROR ("fscanf failed", GSL_EFAILED) ;
	}
    }

  return 0;
} 


