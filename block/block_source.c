size_t
FUNCTION(gsl_block,size) (const TYPE(gsl_block) * b)
{
  return b->size ;
}

ATOMIC *
FUNCTION(gsl_block,data) (const TYPE(gsl_block) * b)
{
  return b->data ;
}
