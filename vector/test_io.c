void FUNCTION (test, text) (void);

void
FUNCTION (test, text) (void)
{
  TYPE (gsl_vector) * v, *w;
  size_t i;

  v = FUNCTION (gsl_vector, alloc) (N);

  {
    FILE *f = fopen ("test.txt", "w");

    for (i = 0; i < N; i++)
      {
	FUNCTION (gsl_vector, set) (v, i, (ATOMIC) i);
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
	if (w->data[i] != (ATOMIC) i)
	  status = 1;
      };

    gsl_test (status, NAME (gsl_vector) "_fprintf and fscanf work correctly");

    fclose (f);
  }

  FUNCTION (gsl_vector, free) (v);
  FUNCTION (gsl_vector, free) (w);
}


