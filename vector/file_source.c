int
FUNCTION (gsl_vector, fread) (FILE * stream, TYPE (gsl_vector) * v)
{
  int status = FUNCTION (gsl_block, raw_fread) (stream,
                                                v->data,
                                                v->size,
                                                v->stride);
  return status;
}

int
FUNCTION (gsl_vector, fwrite) (FILE * stream, const TYPE (gsl_vector) * v)
{
  int status = FUNCTION (gsl_block, raw_fwrite) (stream,
                                                 v->data,
                                                 v->size,
                                                 v->stride);
  return status;
}

#if !(defined(USES_LONGDOUBLE) && !defined(HAVE_PRINTF_LONGDOUBLE))
int
FUNCTION (gsl_vector, fprintf) (FILE * stream, const TYPE (gsl_vector) * v,
				const char *format)
{
  int status = FUNCTION (gsl_block, raw_fprintf) (stream,
                                                  v->data,
                                                  v->size,
                                                  v->stride,
                                                  format);
  return status;
}

int
FUNCTION (gsl_vector, fscanf) (FILE * stream, TYPE (gsl_vector) * v)
{
  int status = FUNCTION (gsl_block, raw_fscanf) (stream,
                                                 v->data,
                                                 v->size,
                                                 v->stride);
  return status;
}
#endif

