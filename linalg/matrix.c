static inline gsl_matrix
submatrix (gsl_matrix * m, size_t i1, size_t j1, size_t i2, size_t j2)
{
  gsl_matrix subm = {0,0,0,0,0} ;

  if (i1 >= m->size1 || i2 >= m->size1 || j1 >= m->size2 || j2 >= m->size2) 
    return subm;

  if (i1 > i2 || j1 > j2) 
    return subm;

  subm.data = m->data + i1 * m->tda + j1 ;
  subm.size1 = (i2 - i1) + 1 ;
  subm.size2 = (j2 - j1) + 1 ;
  subm.tda = m->tda ;
  subm.block = 0 ;

  return subm;
}

static inline gsl_vector
subvector (gsl_vector * v, size_t i1, size_t i2)
{
  gsl_vector subv = {0,0,0,0} ;

  if (i1 >= m->size1 || i2 >= m->size1 || i1 > i2) 
    return subv;

  subv.data = v->data + i1 * v->stride ;
  subv.size = (i2 - i1) + 1 ;
  subv.stride = v->stride ;
  subv.block = 0 ;

  return subv;
}

static inline gsl_vector
col (gsl_matrix *m, size_t i1, size_t j1, size_t i2) 
{
  gsl_vector v = {0,0,0,0} ;

  if (i1 >= m->size1 || i2 >= m->size1 || j1 >= m->size2 || i1 > i2) 
    return v;

  v.data = m->data + i1 * m->tda + j1 ;
  v.size = (i2 - i1) + 1;
  v.stride = m->tda ;
  v.block = 0 ;

  return v;
}

static inline gsl_vector
row (gsl_matrix *m, size_t i1, size_t j1, size_t j2) 
{
  gsl_vector v = {0,0,0,0} ;

  if (i1 >= m->size1 || j1 >= m->size2  || j2 >= m->size2 || j1 > j2) 
    return v;

  v.data = m->data + i1 * m->tda + j1 ;
  v.size = (j2 - j1) + 1;
  v.stride = 1 ;
  v.block = 0 ;

  return v;
}

