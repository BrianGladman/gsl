TYPE (gsl_matrix) *
FUNCTION (gsl_matrix, alloc_from_block) (TYPE(gsl_block) * block, 
                                         const size_t offset,
                                         const size_t n1, 
                                         const size_t n2,
                                         const size_t d2)
{
  TYPE (gsl_matrix) * m;

  if (n1 == 0)
    {
      GSL_ERROR_RETURN ("matrix dimension n1 must be positive integer",
			GSL_EDOM, 0);
    }
  else if (n2 == 0)
    {
      GSL_ERROR_RETURN ("matrix dimension n2 must be positive integer",
			GSL_EDOM, 0);
    }
  else if (d2 < n2)
    {
      GSL_ERROR_RETURN ("matrix dimension d2 must be greater than n2",
			GSL_EDOM, 0);
    }
  else if (block->size < offset + n1 * d2)
    {
      GSL_ERROR_RETURN ("matrix size exceeds available block size",
			GSL_EDOM, 0);
    }

  m = (TYPE (gsl_matrix) *) malloc (sizeof (TYPE (gsl_matrix)));

  if (m == 0)
    {
      GSL_ERROR_RETURN ("failed to allocate space for matrix struct",
			GSL_ENOMEM, 0);
    }

  m->data = block->data + offset;
  m->size1 = n1;
  m->size2 = n2;
  m->dim2 = d2;

  return m;
}


TYPE (gsl_matrix) *
FUNCTION (gsl_matrix, alloc_from_matrix) (TYPE(gsl_matrix) * mm, 
                                          const size_t k1,
                                          const size_t k2,
                                          const size_t n1, 
                                          const size_t n2)
{
  TYPE (gsl_matrix) * m;

  if (n1 == 0)
    {
      GSL_ERROR_RETURN ("matrix dimension n1 must be positive integer",
			GSL_EDOM, 0);
    }
  else if (n2 == 0)
    {
      GSL_ERROR_RETURN ("matrix dimension n2 must be positive integer",
			GSL_EDOM, 0);
    }
  else if (k1 + n1 > mm->size1)
    {
      GSL_ERROR_RETURN ("submatrix dimension 1 exceeds size of original",
			GSL_EDOM, 0);
    }
  else if (k2 + n2 > mm->size2)
    {
      GSL_ERROR_RETURN ("submatrix dimension 2 exceeds size of original",
			GSL_EDOM, 0);
    }

  m = (TYPE (gsl_matrix) *) malloc (sizeof (TYPE (gsl_matrix)));

  if (m == 0)
    {
      GSL_ERROR_RETURN ("failed to allocate space for matrix struct",
			GSL_ENOMEM, 0);
    }

  m->data = mm->data + k1 * mm-> dim2 + k2 ;
  m->size1 = n1;
  m->size2 = n2;
  m->dim2 = mm->dim2;

  return m;
}


void
FUNCTION (gsl_matrix, free) (TYPE (gsl_matrix) * m)
{
  free (m);
}
