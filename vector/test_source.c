void FUNCTION (test, func) (void);
void FUNCTION (test, text) (void);
void FUNCTION (test, binary) (void);

void
FUNCTION (test, func) (void)
{
  TYPE (gsl_vector) * v;
  size_t i;

  v = FUNCTION (gsl_vector, alloc) (N);

  gsl_test (v->data == 0, NAME (gsl_vector) "_alloc returns valid pointer");
  gsl_test (v->size != N, NAME (gsl_vector) "_alloc returns valid size");
  gsl_test (v->stride != 1, NAME (gsl_vector) "_alloc returns unit stride");

  for (i = 0; i < N; i++)
    {
      FUNCTION (gsl_vector, set) (v, i, (ATOMIC) i);
    }


  status = 0;
  for (i = 0; i < N; i++)
    {
      if (v->data[i] != (ATOMIC) i)
	status = 1;
    };
  
  gsl_test (status,
	    NAME (gsl_vector) "_set" DESC " writes into array correctly");


  status = 0;
  for (i = 0; i < N; i++)
    {
      if (FUNCTION (gsl_vector, get) (v, i) != (ATOMIC) i)
	status = 1;
    };
  gsl_test (status,
	    NAME (gsl_vector) "_get" DESC " reads from array correctly");

  /* Now set stride to 2 */

  v->stride = 2;
  
  status = 0;
  for (i = 0; i < N / 2; i++)
    {
      if (FUNCTION (gsl_vector, get) (v, i) != (ATOMIC) (2 * i))
	status = 1;
    };
  gsl_test (status, NAME (gsl_vector) "_get" DESC " reads correctly with stride");
  
  for (i = 0; i < N / 2; i++)
    {
      FUNCTION (gsl_vector, set) (v, i, (ATOMIC) (i + 1000));
    };
  
  status = 0;

  for (i = 0; i < N / 2; i++)
    {
      if (v->data[2 * i] != (ATOMIC) (i + 1000))
	status = 1;
    };
  
  gsl_test (status, NAME (gsl_vector) "_set" DESC " writes correctly with stride");

  FUNCTION (gsl_vector, free) (v);	/* free whatever is in v */

  v = FUNCTION (gsl_vector, calloc) (N);

  gsl_test (v->data == 0, NAME (gsl_vector) "_calloc returns valid pointer");
  gsl_test (v->size != N, NAME (gsl_vector) "_calloc returns valid size");

  status = 0;

  for (i = 0; i < N; i++)
    {
      if (v->data[i] != 0.0)
	status = 1;
    };
  
  gsl_test (status, NAME (gsl_vector) "_calloc initializes array to zero");

}

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
}

void
FUNCTION (test, binary) (void)
{
  TYPE (gsl_vector) * v = FUNCTION (gsl_vector, alloc) (N);
  TYPE (gsl_vector) * w = FUNCTION (gsl_vector, alloc) (N);

  size_t i;

  {
    FILE *f = fopen ("test.dat", "w");

    for (i = 0; i < N; i++)
      {
	FUNCTION (gsl_vector, set) (v, i, (ATOMIC) (N - i));
      };

    FUNCTION (gsl_vector, fwrite) (f, v);

    fclose (f);
  }

  {
    FILE *f = fopen ("test.dat", "r");

    FUNCTION (gsl_vector, fread) (f, w);

    status = 0;
    for (i = 0; i < N; i++)
      {
	if (w->data[i] != (ATOMIC) (N - i))
	  status = 1;
      };

    gsl_test (status, NAME (gsl_vector) "_write and read work correctly");

    fclose (f);
  }

}

void FUNCTION (test, trap) (void);

void
FUNCTION (test, trap) (void)
{
  TYPE (gsl_vector) * v = FUNCTION (gsl_vector, alloc) (N);

  size_t j = 0;
  double x;


  status = 0;
  FUNCTION (gsl_vector, set) (v, j - 1, 1.2);
  gsl_test (!status, NAME (gsl_vector) "_set traps index below lower bound");

  status = 0;
  FUNCTION (gsl_vector, set) (v, N + 1, 1.2);
  gsl_test (!status, NAME (gsl_vector) "_set traps index above upper bound");

  status = 0;
  FUNCTION (gsl_vector, set) (v, N, 1.2);
  gsl_test (!status, NAME (gsl_vector) "_set traps index at upper bound");

  status = 0;
  x = FUNCTION (gsl_vector, get) (v, j - 1);
  gsl_test (!status, NAME (gsl_vector) "_get traps index below lower bound");
  gsl_test (x != 0,
	 NAME (gsl_vector) "_get returns zero for index below lower bound");

  status = 0;
  x = FUNCTION (gsl_vector, get) (v, N + 1);
  gsl_test (!status, NAME (gsl_vector) "_get traps index above upper bound");
  gsl_test (x != 0,
	 NAME (gsl_vector) "_get returns zero for index above upper bound");

  status = 0;
  x = FUNCTION (gsl_vector, get) (v, N);
  gsl_test (!status, NAME (gsl_vector) "_get traps index at upper bound");
  gsl_test (x != 0,
	    NAME (gsl_vector) "_get returns zero for index at upper bound");
}
