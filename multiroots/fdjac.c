#include <config.h>
#include <gsl/gsl_multiroots.h>

int
gsl_multiroot_fdjacobian (gsl_multiroot_function * F,
                           const gsl_vector * x, const gsl_vector * f,
                           double epsrel, gsl_matrix * jacobian)
{
  const size_t n = x->size;
  const size_t m = f->size;
  const size_t n1 = jacobian->size1;
  const size_t n2 = jacobian->size2;

  if (m != n1 || n != n2)
    {
      GSL_ERROR ("function and jacobian are not conformant", GSL_EBADLEN);
    }

  {
    size_t i,j;
    gsl_vector *x1, *f1;

    x1 = gsl_vector_alloc (n);

    if (x1 == 0)
      {
	GSL_ERROR ("failed to allocate space for x1 workspace", GSL_ENOMEM);
      }

    f1 = gsl_vector_alloc (m);

    if (f1 == 0)
      {
	gsl_vector_free (x1);

	GSL_ERROR ("failed to allocate space for f1 workspace", GSL_ENOMEM);
      }

    gsl_vector_cpy (x1, x);	/* copy x into x1 */

    for (j = 0; j < n; j++)
      {
	double xj = gsl_vector_get (x, j);
	double dx = epsrel * fabs (xj);

	if (dx == 0)
	  {
	    dx = epsrel;
	  }

	gsl_vector_set (x1, j, xj + dx);
	GSL_MULTIROOT_FN_EVAL (F, x1, f1);
	gsl_vector_set (x1, j, xj);

	for (i = 0; i < m; i++)
	  {
	    double g1 = gsl_vector_get (f1, i);
	    double g0 = gsl_vector_get (f, i);
	    gsl_matrix_set (jacobian, i, j, (g1 - g0) / dx);
	  }
      }

    gsl_vector_free (x1);
    gsl_vector_free (f1);
  }
  
  return GSL_SUCCESS;
}
