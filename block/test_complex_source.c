void FUNCTION (test, func) (void);
void FUNCTION (test, binary) (void);
void FUNCTION (test, trap) (void);

void
FUNCTION (test, func) (void)
{
  TYPE (gsl_block) * b;
  ATOMIC * data;
  size_t i, size;

  b = FUNCTION (gsl_block, alloc) (N);

  gsl_test (b->data == 0, NAME (gsl_block) "_alloc returns valid pointer");
  gsl_test (b->size != N, NAME (gsl_block) "_alloc returns valid size");

  data = FUNCTION(gsl_block, data) (b);
  size = FUNCTION(gsl_block, size) (b);

  gsl_test (data == 0, NAME (gsl_block) "_data returns valid pointer");
  gsl_test (size != N, NAME (gsl_block) "_size returns valid size");
  
  FUNCTION (gsl_block, free) (b);	/* free whatever is in v */

  b = FUNCTION (gsl_block, calloc) (N);

  gsl_test (b->data == 0, NAME (gsl_block) "_calloc returns valid pointer");
  gsl_test (b->size != N, NAME (gsl_block) "_calloc returns valid size");

  data = FUNCTION(gsl_block, data) (b);
  size = FUNCTION(gsl_block, size) (b);

  gsl_test (data == 0, NAME (gsl_block) "_data returns valid pointer from calloc");
  gsl_test (size != N, NAME (gsl_block) "_size returns valid size from calloc");

  status = 0;
  
  for (i = 0; i < N; i++)
    {
      if (b->data[2 * i] != 0.0 || b->data[2 * i + 1] != 0.0)
	status = 1;
    };

  gsl_test (status, NAME (gsl_block) "_calloc initializes array to zero");

  FUNCTION (gsl_block, free) (b);
}

void
FUNCTION (test, binary) (void)
{
  size_t i;

  {
    TYPE (gsl_block) * v = FUNCTION (gsl_block, calloc) (N);

    FILE *f = fopen ("test.dat", "w");

    for (i = 0; i < N; i++)
      {
	v->data[2*i] = N - i;
        v->data[2*i + 1] = 10*(N-i) + 1;
      };

    FUNCTION (gsl_block, fwrite) (f, v);

    fclose (f);
    
    FUNCTION (gsl_block, free) (v);
  }

  {
    TYPE (gsl_block) * w = FUNCTION (gsl_block, calloc) (N);

    FILE *f = fopen ("test.dat", "r");

    FUNCTION (gsl_block, fread) (f, w);

    status = 0;
    for (i = 0; i < N; i++)
      {
	if (w->data[2 * i] != (ATOMIC) (N - i) || w->data[2 * i + 1] != (ATOMIC) (10*(N - i) + 1))
	  status = 1;
      };
    fclose (f);

    FUNCTION (gsl_block, free) (w);
  }

  gsl_test (status, NAME (gsl_block) "_write and read work correctly");

}


void
FUNCTION (test, trap) (void)
{
  TYPE (gsl_block) * b = FUNCTION (gsl_block, alloc) (0);

  gsl_test (b != 0, NAME (gsl_block) "_alloc traps zero length");
}




