#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

int
main ()
{
  double data[] = { 1.0  , 1/2.0, 1/3.0, 1/4.0,
                    1/2.0, 1/3.0, 1/4.0, 1/5.0,
                    1/3.0, 1/4.0, 1/5.0, 1/6.0,
                    1/4.0, 1/5.0, 1/6.0, 1/7.0 };

  gsl_matrix m = gsl_matrix_view(data, 4, 4);

  gsl_vector *eval = gsl_vector_alloc (4);
  gsl_matrix *evec = gsl_matrix_alloc (4, 4);

  gsl_eigen_symmv_workspace * w = 
    gsl_eigen_symmv_alloc (4);
  
  gsl_eigen_symmv (&m, eval, evec, w);

  gsl_eigen_symmv_free(w);

  gsl_eigen_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);
  
  {
    int i;

    for (i = 0; i < 4; i++)
      {
        gsl_vector v = gsl_matrix_column(evec, i);
        printf("eigenvalue = %g\n", gsl_vector_get(eval, i));
        printf("eigenvector = \n");
        gsl_vector_fprintf(stdout, &v, "%g");
      }
  }
}
