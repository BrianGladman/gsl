void FUNCTION (test, text) (void);

void
FUNCTION (test, text) (void)
{
  size_t i;
  
  {
    TYPE (gsl_block) * v = FUNCTION (gsl_block, alloc) (N);

    FILE *f = fopen ("test.txt", "w");

    for (i = 0; i < N; i++)
      {
	v->data[i] = (ATOMIC) i;
      };

    FUNCTION (gsl_block, fprintf) (f, v, OUT_FORMAT);

    fclose (f);

    FUNCTION (gsl_block, free) (v);
  }

  {
    TYPE (gsl_block) * w = FUNCTION (gsl_block, calloc) (N);

    FILE *f = fopen ("test.txt", "r");

    FUNCTION (gsl_block, fscanf) (f, w);

    status = 0;
    for (i = 0; i < N; i++)
      {
	if (w->data[i] != (ATOMIC) i)
	  status = 1;
      };

    gsl_test (status, NAME (gsl_block) "_fprintf and fscanf work correctly");

    fclose (f);

    FUNCTION (gsl_block, free) (w);
  }
}


