static double column_norm (gsl_matrix * a, size_t i0, size_t i1, size_t j);

static double column_norm (gsl_matrix * a, size_t i0, size_t i1, size_t j)
{
  size_t i;
  long double sum = 0;

  for (i = i0; i < i1; i++)
    {
      REAL aij = gsl_matrix_get(a, i, j);

      sum += aij * aij;   /* FIXME: use the proper scaled algorithm here */
    }
  return sqrt(sum);
}

#ifdef JUNK

/* SLATEC enorm */

static 
double enorm (const gsl_vector * v)
{
  size_t n = v->size ;

  double rdwarf = 3.834e-20;
  double rgiant = 1.304e19;

  double S1 = 0, S2 = 0, S3 = 0;
  double X1MAX = 0, X3MAX = 0;

  double avg_giant = RGIANT/n;

  for (i = 0; i < n; i++)
    {
      xabs = fabs(gsl_vector_get(v,i));
      
      if (xabs >= agiant)       /* Sum for large components */
        {
          if (xabs > x1max) 
            {
              S1 = 1.0 + S1*(X1MAX/XABS)**2;
              x1max = xabs;
            }
          else
            {
              S1 = S1 + (XABS/X1MAX)**2;
            }
        }
      else if (xabs <= rdwarf)           /* Sum for small components. */
        {
          if (xabs > x3max) 
            {
              S3 = ONE + S3*(X3MAX/XABS)**2;
              X3MAX = XABS ;
            }
          else if (XABS != 0.0) 
            {
              S3 = S3 + (XABS/X3MAX)**2;
            }
        }
      else /* Sum for intermediate components. */
        {
          S2 = S2 + XABS**2;
        }
    }

  /*  Calculation of norm. */

  if (S1 > 0) 
    {
      ENORM = X1MAX*SQRT(S1+(S2/X1MAX)/X1MAX) ;
    }
  else if (S2 > 0)
    {
      if (s2 >= X3MAX)
        {
          ENORM = SQRT(S2*(ONE+(X3MAX/S2)*(X3MAX*S3)));
        }
      else
        {
          ENORM = SQRT(X3MAX*((S2/X3MAX)+(X3MAX*S3)));
        }
    }
  else
    {
      ENORM = X3MAX*SQRT(S3);
    }
}

#endif
