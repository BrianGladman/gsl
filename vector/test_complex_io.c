void FUNCTION (test, text) (void);


void
FUNCTION (test, text) (void)
{
  TYPE (gsl_block) * bv = FUNCTION (gsl_block, alloc) (N);
  TYPE (gsl_block) * bw = FUNCTION (gsl_block, alloc) (N);
  TYPE (gsl_vector) * v = FUNCTION (gsl_vector, alloc) (bv,0,N,1);
  TYPE (gsl_vector) * w = FUNCTION (gsl_vector, alloc) (bw,0,N,1);

  size_t i;

  {
    FILE *f = fopen ("test.txt", "w");

    for (i = 0; i < N; i++)
      {
	BASE x;
	GSL_REAL (x) = i;
	GSL_IMAG (x) = i + 1;
	FUNCTION (gsl_vector, set) (v, i, x);
      };

    FUNCTION (gsl_vector, fprintf) (f, v, OUT_FORMAT);

    fclose (f);
  }

  {
    FILE *f = fopen ("test.txt", "r");

    FUNCTION (gsl_vector, fscanf) (f, w);

    status = 0;
    for (i = 0; i < N; i++)
      {
	if (w->data[2 * i] != (ATOMIC) i || w->data[2 * i + 1] != (ATOMIC) (i + 1))
	  status = 1;
      };
    fclose (f);
  }

  FUNCTION (gsl_vector, free) (v);
  FUNCTION (gsl_vector, free) (w);
  FUNCTION (gsl_block, free) (bv);
  FUNCTION (gsl_block, free) (bw);

  gsl_test (status, NAME (gsl_vector) "_fprintf and fscanf work correctly");
}
