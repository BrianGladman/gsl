int
FUNCTION (gsl_matrix, isnull) (const TYPE (gsl_matrix) * m)
{
  const size_t size1 = m->size1;
  const size_t size2 = m->size2;
  const size_t tda = m->tda ;
  
  size_t i, j, k;

  for (i = 0; i < size1 ; i++)
    {
      for (j = 0; j < size2; j++)
        {
          for (k = 0; k < MULTIPLICITY; k++) 
            {
              if (m->data[(i * tda + j) * MULTIPLICITY + k] != 0.0)
                {
                  return 0;
                }
            }
        }
    }
      
  return 1;
}

