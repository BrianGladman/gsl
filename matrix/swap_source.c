int
FUNCTION (gsl_matrix, swap_rows) (TYPE (gsl_matrix) * m,
				 const size_t i, const size_t j)
{
  const size_t size1 = m->size1;
  const size_t size2 = m->size2;

  if (i >= size1)
    {
      GSL_ERROR ("first row index is out of range", GSL_EINVAL);
    }

  if (j >= size1)
    {
      GSL_ERROR ("second row index is out of range", GSL_EINVAL);
    }

  if (i != j)
    {
      ATOMIC *row1 = m->data + MULTIPLICITY * i * m->tda;
      ATOMIC *row2 = m->data + MULTIPLICITY * j * m->tda;
      
      size_t k;
      
      for (k = 0; k < MULTIPLICITY * size2; k++)
        {
          ATOMIC tmp = row1[k] ;
          row1[k] = row2[k] ;
          row2[k] = tmp ;
        }
    }

  return GSL_SUCCESS;
}

int
FUNCTION (gsl_matrix, swap_columns) (TYPE (gsl_matrix) * m,
                                     const size_t i, const size_t j)
{
  const size_t size1 = m->size1;
  const size_t size2 = m->size2;

  if (i >= size2)
    {
      GSL_ERROR ("first column index is out of range", GSL_EINVAL);
    }

  if (j >= size2)
    {
      GSL_ERROR ("second column index is out of range", GSL_EINVAL);
    }

  if (i != j)
    {
      ATOMIC *col1 = m->data + MULTIPLICITY * i;
      ATOMIC *col2 = m->data + MULTIPLICITY * j;
      
      size_t p;
      
      for (p = 0; p < size1; p++)
        {
          size_t k;
          size_t n = p * MULTIPLICITY * m->tda;
 
          for (k = 0; k < MULTIPLICITY; k++)
            {
              ATOMIC tmp = col1[n+k] ;
              col1[n+k] = col2[n+k] ;
              col2[n+k] = tmp ;
            }
        }
    }

  return GSL_SUCCESS;
}


int
FUNCTION (gsl_matrix, swap_rowcol) (TYPE (gsl_matrix) * m,
                                    const size_t i, const size_t j)
{
  const size_t size1 = m->size1;
  const size_t size2 = m->size2;

  if (size1 != size2)
    {
      GSL_ERROR ("matrix must be square to swap row and column", GSL_ENOTSQR);
    }

  if (i >= size1)
    {
      GSL_ERROR ("row index is out of range", GSL_EINVAL);
    }

  if (j >= size2)
    {
      GSL_ERROR ("column index is out of range", GSL_EINVAL);
    }

  {
    ATOMIC *row = m->data + MULTIPLICITY * i * m->tda;
    ATOMIC *col = m->data + MULTIPLICITY * j;
      
    size_t p;
    
    for (p = 0; p < size1; p++)
      {
        size_t k;

        size_t r = p * MULTIPLICITY;
        size_t c = p * MULTIPLICITY * m->tda;
        
          for (k = 0; k < MULTIPLICITY; k++)
            {
              ATOMIC tmp = col[c+k] ;
              col[c+k] = row[r+k] ;
              row[r+k] = tmp ;
            }
        }
    }

  return GSL_SUCCESS;
}

