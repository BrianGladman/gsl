inline BASE
FUNCTION (gsl_vector, get) (const TYPE (gsl_vector) * v, const size_t i)
{
  if (gsl_check_range)
    {
      if (i >= v->size)		/* size_t is unsigned, can't be negative */
	{
	  const BASE zero = ZERO;
	  GSL_ERROR_RETURN ("index out of range", GSL_EINVAL, zero);
	}
    }

  /* The following line is a generalization of return v->data[i] */

  return *(BASE *) (v->data + MULTIPLICITY * i * v->stride);
}

inline void
FUNCTION (gsl_vector, set) (TYPE (gsl_vector) * v, const size_t i, BASE x)
{
  if (gsl_check_range)
    {
      if (i >= v->size)		/* size_t is unsigned, can't be negative */
	{
	  GSL_ERROR_RETURN_NOTHING ("index out of range", GSL_EINVAL);
	}
    }

  /* The following line is a generalization of v->data[i] = x */

  *(BASE *) (v->data + MULTIPLICITY * i * v->stride) = x;
}

void
FUNCTION (gsl_vector, set_all) (TYPE (gsl_vector) * v, BASE x)
{
  ATOMIC * const data = v->data;
  const size_t n = v->size;
  const size_t stride = v->stride;

  size_t i;

  for (i = 0; i < n; i++)
    {
      *(BASE *) (data + MULTIPLICITY * i * stride) = x;
    }
}

