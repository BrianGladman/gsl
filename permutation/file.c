#include <config.h>
#include <stdio.h>
#include <gsl_errno.h>
#include <gsl_permutation.h>

#define IN_FORMAT "%lu"

int
gsl_permutation_fread (FILE * stream, gsl_permutation * p)
{
  size_t n = p->size ;

  size_t * data = p->data ;

  size_t items = fread (data, sizeof (size_t), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fread failed", GSL_EFAILED);
    }
      
  return GSL_SUCCESS;
}

int
gsl_permutation_fwrite (FILE * stream, const gsl_permutation * p)
{
  size_t n = p->size ;

  size_t * data = p->data ;
  
  size_t items = fwrite (data, sizeof (size_t), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fwrite failed", GSL_EFAILED);
    }

  return GSL_SUCCESS;
}

int
gsl_permutation_fprintf (FILE * stream, const gsl_permutation * p, const char *format)
{
  size_t n = p->size ;
  
  size_t * data = p->data ;
  
  size_t i;

  for (i = 0; i < n; i++)
    {
      int status = fprintf (stream, format, data[i]);

      if (status < 0)
        {
          GSL_ERROR ("fprintf failed", GSL_EFAILED);
        }
    }

  return GSL_SUCCESS;
}

int
gsl_permutation_fscanf (FILE * stream, gsl_permutation * p)
{
  size_t n = p->size ;
  
  size_t * data = p->data ;

  size_t i;

  for (i = 0; i < n; i++)
    {
      unsigned long j ;  

      /* FIXME: what if size_t != unsigned long ??? 

         want read in size_t but have to read in unsigned long to avoid
         error from compiler */

      int status = fscanf (stream, IN_FORMAT, &j);  

      if (status != 1)
        {
          GSL_ERROR ("fscanf failed", GSL_EFAILED);
        }

      data[i] = j;
    }

  return GSL_SUCCESS;
}

