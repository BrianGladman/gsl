void FUNCTION (test, text) (void);


void
FUNCTION (test, text) (void)
{
  TYPE (gsl_vector) * w, *v = FUNCTION (gsl_vector, alloc) (N);

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

  w = FUNCTION (gsl_vector, calloc) (N);

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

  gsl_test (status, NAME (gsl_vector) "_fprintf and fscanf work correctly");
}
