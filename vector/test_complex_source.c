int FUNCTION (test, func) (void);

int
  FUNCTION (test, func) (void)
{
  TYPE (gsl_vector) * v, *w;
  size_t i;

  v = FUNCTION (gsl_vector, alloc) (N);

  gsl_test (v->data == 0, NAME (gsl_vector) "_alloc returns valid pointer");
  gsl_test (v->size != N, NAME (gsl_vector) "_alloc returns valid size");

  for (i = 0; i < N; i++)
    {
      BASE x = {{i,i+1}} ;
      FUNCTION (gsl_vector, set) (v, i, x);
    };

  {
    int status = 0;

    for (i = 0; i < N; i++)
      {
	if (v->data[2*i] != (ATOMIC) i || v->data[2*i+1] != (ATOMIC) (i+1))
	  status = 1;
      };

    gsl_test (status, NAME (gsl_vector) "_set writes into array correctly");
  }

  {
    int status = 0;

    for (i = 0; i < N; i++)
      {
	BASE x = {{i,i+1}} ;
	BASE y = FUNCTION (gsl_vector, get) (v, i); 
	if (!GSL_COMPLEX_EQ(x,y) )
	  status = 1;
      };
    gsl_test (status, NAME (gsl_vector) "_get reads from array correctly");
  }

  FUNCTION (gsl_vector, free) (v);	/* free whatever is in v */

  v = FUNCTION (gsl_vector, calloc) (N);

  gsl_test (v->data == 0, NAME (gsl_vector) "_calloc returns valid pointer");
  gsl_test (v->size != N, NAME (gsl_vector) "_calloc returns valid size");

  {
    int status = 0;

    for (i = 0; i < N; i++)
      {
	if (v->data[2*i] != 0.0 || v->data[2*i+1] != 0.0)
	  status = 1;
      };

    gsl_test (status, NAME (gsl_vector) "_calloc initializes array to zero");
  }

  {
    FILE *f = fopen ("test.txt", "w");

    for (i = 0; i < N; i++)
      {
	BASE x = {{i,i+1}} ;
	FUNCTION (gsl_vector, set) (v, i, x);
      };

    FUNCTION (gsl_vector, fprintf) (f, v, OUT_FORMAT);

    fclose (f);
  }

  w = FUNCTION (gsl_vector, calloc) (N);

  {
    int status = 0;
    FILE *f = fopen ("test.txt", "r");

    FUNCTION (gsl_vector, fscanf) (f, w);

    for (i = 0; i < N; i++)
      {
	if (w->data[2*i] != (ATOMIC) i || w->data[2*i+1] != (ATOMIC) (i+1))
	  status = 1;
      };

    gsl_test (status, NAME (gsl_vector) "_fprintf and fscanf work correctly");

    fclose (f);
  }


  {
    FILE *f = fopen ("test.dat", "w");

    for (i = 0; i < N; i++)
      {
	BASE x = {{N-i,N-i+1}} ;
	FUNCTION (gsl_vector, set) (v, i, x);
      };

    FUNCTION (gsl_vector, fwrite) (f, v);

    fclose (f);
  }

  {
    int status = 0;
    FILE *f = fopen ("test.dat", "r");

    FUNCTION (gsl_vector, fread) (f, w);

    for (i = 0; i < N; i++)
      {
	if (w->data[2*i] != (ATOMIC) (N - i) || w->data[2*i+1] != (ATOMIC) (N - i + 1)  )
	  status = 1;
      };

    gsl_test (status, NAME (gsl_vector) "_write and read work correctly");

    fclose (f);
  }

  return gsl_test_summary ();
}

