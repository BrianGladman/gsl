int
FUNCTION (gsl_vector, isnull) (TYPE (gsl_vector) * v)
{
  const size_t n = v->size;
  const size_t stride = v->stride ;
  
  size_t j;

  for (j = 0; j < n; j++)
    {
      size_t k;
      
      for (k = 0; k < MULTIPLICITY; k++) 
        {
          if (v->data[stride * j + k] != 0.0)
            {
              return 0;
            }
        }
    }

  return 1;
}

