void FUNCTION (test, func) (void);
void FUNCTION (test, binary) (void);
void FUNCTION (test, trap) (void);

void
FUNCTION (test, func) (void)
{
  TYPE (gsl_block) * b;
  TYPE (gsl_vector) * v;
  size_t i;

  b = FUNCTION (gsl_block, alloc) (N);
  v = FUNCTION (gsl_vector, alloc_from_block) (b, 0, N, 1);

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
  FUNCTION (gsl_block, free) (b);

}


void
FUNCTION (test, binary) (void)
{
  TYPE (gsl_block) * bv = FUNCTION (gsl_block, alloc) (N);
  TYPE (gsl_block) * bw = FUNCTION (gsl_block, alloc) (N);

  TYPE (gsl_vector) * v = FUNCTION (gsl_vector, alloc_from_block) (bv,0,N,1);
  TYPE (gsl_vector) * w = FUNCTION (gsl_vector, alloc_from_block) (bw,0,N,1);

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

  FUNCTION (gsl_vector, free) (v);	/* free whatever is in v */
  FUNCTION (gsl_vector, free) (w);	/* free whatever is in w */
  FUNCTION (gsl_block, free) (bv);
  FUNCTION (gsl_block, free) (bw);
}

void
FUNCTION (test, trap) (void)
{
  TYPE (gsl_block) * b = FUNCTION (gsl_block, alloc) (N);
  TYPE (gsl_vector) * v = FUNCTION (gsl_vector, alloc_from_block) (b,0,N,1);

  size_t j = 0;
  double x;
  BASE * y;

  status = 0;
  FUNCTION (gsl_vector, set) (v, j - 1, (ATOMIC)0);
  gsl_test (!status, NAME (gsl_vector) "_set traps index below lower bound");

  status = 0;
  FUNCTION (gsl_vector, set) (v, N + 1, (ATOMIC)0);
  gsl_test (!status, NAME (gsl_vector) "_set traps index above upper bound");

  status = 0;
  FUNCTION (gsl_vector, set) (v, N, (ATOMIC)0);
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

  FUNCTION (gsl_vector, free) (v);	/* free whatever is in v */
  FUNCTION (gsl_block, free) (b);
}





