int
FUNCTION (gsl_matrix, fread) (FILE * stream, TYPE (gsl_matrix) * m)
{
  int status = 0;

  const size_t size1 = m->size1; 
  const size_t size2 = m->size2;

  const size_t tda = m->tda;

  if (tda == size2) /* the rows are contiguous */
    {
      status = FUNCTION (gsl_block, raw_fread) (stream, 
                                                m->data, 
                                                size1 * size2, 1);
    }
  else
    {
      size_t i;

      for (i = 0 ; i < size1 ; i++)  /* read each row separately */
        {
          status = FUNCTION (gsl_block, raw_fread) (stream, 
                                                    m->data + i * tda, 
                                                    size2, 1);
          if (status)
            break;
        }
    }

  return status;
}

int
FUNCTION (gsl_matrix, fwrite) (FILE * stream, const TYPE (gsl_matrix) * m)
{
  int status = 0;

  const size_t size1 = m->size1; 
  const size_t size2 = m->size2;

  const size_t tda = m->tda;

  if (tda == size2) /* the rows are contiguous */
    {
      status = FUNCTION (gsl_block, raw_fwrite) (stream, 
                                                 m->data, 
                                                 size1 * size2, 1);
    }
  else
    {
      size_t i;

      for (i = 0 ; i < size1 ; i++)  /* write each row separately */
        {
          status = FUNCTION (gsl_block, raw_fwrite) (stream, 
                                                     m->data + i * tda, 
                                                     size2, 1);
          if (status)
            break;
        }
    }

  return status;

}

#if !(defined(USES_LONGDOUBLE) && !defined(HAVE_PRINTF_LONGDOUBLE))
int
FUNCTION (gsl_matrix, fprintf) (FILE * stream, const TYPE (gsl_matrix) * m,
				const char *format)
{
  int status = 0;

  const size_t size1 = m->size1; 
  const size_t size2 = m->size2;

  const size_t tda = m->tda;

  if (tda == size2) /* the rows are contiguous */
    {
      status = FUNCTION (gsl_block, raw_fprintf) (stream, 
                                                  m->data, 
                                                  size1 * size2, 1,
                                                  format);
    }
  else
    {
      size_t i;

      for (i = 0 ; i < size1 ; i++)  /* print each row separately */
        {
          status = FUNCTION (gsl_block, raw_fprintf) (stream, 
                                                      m->data + i * tda, 
                                                      size2, 1,
                                                      format);
          if (status)
            break;
        }
    }

  return status;
}

int
FUNCTION (gsl_matrix, fscanf) (FILE * stream, TYPE (gsl_matrix) * m)
{
  int status = 0;

  const size_t size1 = m->size1; 
  const size_t size2 = m->size2;

  const size_t tda = m->tda;

  if (tda == size2)  /* the rows are contiguous */
    {
      status = FUNCTION (gsl_block, raw_fscanf) (stream, 
                                                 m->data, 
                                                 size1 * size2, 1);
    }
  else
    {
      size_t i;

      for (i = 0 ; i < size1 ; i++)  /* scan each row separately */
        {
          status = FUNCTION (gsl_block, raw_fscanf) (stream, 
                                                     m->data + i * tda, 
                                                     size2, 1);
          if (status)
            break;
        }
    }

  return status;
}
#endif

