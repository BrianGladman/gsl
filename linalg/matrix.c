gsl_matrix
submatrix (gsl_matrix * m, size_t i1, size_t j1, size_t i2, size_t j2)
{
  gsl_matrix subm ;

  subm.data = m->data + i1 * m->tda + j1 ;
  subm.size1 = (i2 - i1) + 1 ;
  subm.size2 = (j2 - j1) + 1 ;
  subm.tda = m->tda ;
  subm.block = 0 ;

  return subm;
}

gsl_vector
col (gsl_matrix *m, size_t i1, size_t j1, size_t i2) 
{
  gsl_vector v ;

  v.data = m->data + i1 * m->tda + j1 ;
  v.size = (i2 - i1) + 1;
  v.stride = m->tda ;
  v.block = 0 ;

  return v;
}

gsl_vector
row (gsl_matrix *m, size_t i1, size_t j1, size_t j2) 
{
  gsl_vector v ;

  v.data = m->data + i1 * m->tda + j1 ;
  v.size = (j2 - j1) + 1;
  v.stride = 1 ;
  v.block = 0 ;

  return v;
}

