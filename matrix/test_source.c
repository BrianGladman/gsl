void FUNCTION(test,func) (void);
void FUNCTION(test,trap) (void);
void FUNCTION(test,text) (void);
void FUNCTION(test,binary) (void);

 
void 
FUNCTION(test,func) (void) 
{
  TYPE(gsl_matrix) * m;
  size_t i, j;
  int k = 0;

  m = FUNCTION(gsl_matrix,alloc) (N,M) ;

  gsl_test(m->data == 0, NAME(gsl_matrix) "_alloc returns valid pointer") ;
  gsl_test(m->size1 != N, NAME(gsl_matrix) "_alloc returns valid size1") ;
  gsl_test(m->size2 != M, NAME(gsl_matrix) "_alloc returns valid size2") ;

  for (i = 0 ; i < N ; i++) {
    for (j = 0 ; j < M ; j++) {
      k++ ; FUNCTION(gsl_matrix,set)(m,i,j,(BASE)k) ;
    }
  }

  { 
    int status = 0 ;
    k = 0 ;
    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	k++ ;
	if(m->data[i*M + j] != (BASE)k) 
	  status = 1 ;
      } ;
    } ;
    
    gsl_test(status, NAME(gsl_matrix) "_set writes into array correctly") ;
  }

  {
    int status = 0 ;
    k = 0 ;
    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	k++ ;
	if(FUNCTION(gsl_matrix,get)(m,i,j) != (BASE)k) 
	  status = 1 ;
      } ;
    } ;
    gsl_test(status, NAME(gsl_matrix) "_get reads from array correctly") ;
  }

  FUNCTION(gsl_matrix,free)(m);  /* free whatever is in m */
  
  m = FUNCTION(gsl_matrix,calloc)(N,M) ;

  gsl_test(m->data == 0, NAME(gsl_matrix) "_calloc returns valid pointer") ;
  gsl_test(m->size1 != N, NAME(gsl_matrix) "_calloc returns valid size1") ;
  gsl_test(m->size2 != M, NAME(gsl_matrix) "_calloc returns valid size2") ;

  { 
    int status = 0 ;

    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	if(m->data[i*M + j] != 0.0) 
	  status = 1 ;
      } 
    }
    gsl_test(status, NAME(gsl_matrix) "_calloc initializes array to zero") ;
  }
  FUNCTION(gsl_matrix,free) (m) ;
}

void 
FUNCTION(test,text) (void) 
{
  TYPE(gsl_matrix) * m = FUNCTION(gsl_matrix,alloc) (N,M) ; 
  size_t i, j;
  int k = 0 ;

  {
    FILE * f = fopen("test.txt","w") ;

    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	k++ ; FUNCTION(gsl_matrix,set)(m,i,j,(BASE)k) ;
      }
    }
    
    FUNCTION(gsl_matrix,fprintf)(f, m, OUT_FORMAT) ;
    fclose(f) ;
  }

  {
    int status = 0 ;
    FILE * f = fopen("test.txt","r") ;
    TYPE(gsl_matrix) * mm = FUNCTION(gsl_matrix,calloc)(N,M) ;

    FUNCTION(gsl_matrix,fscanf)(f, mm) ;
    k = 0 ;
    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	k++ ; 
	if (mm->data[i*M + j] != (BASE)k)
	  status = 1 ;
      }
    }
    
    gsl_test(status, NAME(gsl_matrix) "_fprintf and fscanf work correctly") ;
    
    fclose (f) ;
    FUNCTION(gsl_matrix,free) (mm) ;
  }
  FUNCTION(gsl_matrix,free) (m) ;
}

void 
FUNCTION(test,binary) (void) 
{
  TYPE(gsl_matrix) * m = FUNCTION(gsl_matrix,alloc) (N,M) ; 
  size_t i, j;
  int k = 0 ;

  {
    FILE * f = fopen("test.dat", "w") ;
    k = 0 ;
    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	k++ ; FUNCTION(gsl_matrix,set)(m,i,j,(BASE)k) ;
      }
    }

    FUNCTION(gsl_matrix,fwrite)(f, m) ;
    fclose (f) ;
  }

  {
    int status = 0 ;
    FILE * f = fopen("test.dat","r") ;
    TYPE(gsl_matrix) * mm = FUNCTION(gsl_matrix,calloc)(N,M) ;

    FUNCTION(gsl_matrix,fread)(f, mm) ;
    k = 0 ;
    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	k++ ; 
	if (mm->data[i*M + j] != (BASE)k)
	    status = 1 ;
      }
    }

    gsl_test(status, NAME(gsl_matrix) "_write and read work correctly") ;

    fclose (f) ;
    FUNCTION(gsl_matrix,free) (mm) ;
  }

  FUNCTION(gsl_matrix,free) (m) ;
}

#undef GSL_RANGE_CHECK_OFF
#define GSL_WARNINGS_OFF

void
FUNCTION(test,trap) (void)
{
  TYPE(gsl_matrix) * m = FUNCTION(gsl_matrix,alloc) (N,M) ;

  size_t i = 0, j = 0;
  double x; 

  status = 0 ;
  FUNCTION(gsl_matrix,set)(m, N+1, 0, 1.2) ;
  gsl_test(!status, 
	   NAME(gsl_matrix) "_set traps 1st index above upper bound") ;

  status = 0 ;
  FUNCTION(gsl_matrix,set)(m, 0, M+1, 1.2) ;
  gsl_test(!status, 
	   NAME(gsl_matrix) "_set traps 2nd index above upper bound") ;

  status = 0 ;
  FUNCTION(gsl_matrix,set)(m, N, 0, 1.2) ;
  gsl_test(!status, 
	   NAME(gsl_matrix) "_set traps 1st index at upper bound") ;

  status = 0 ;
  FUNCTION(gsl_matrix,set)(m, 0, M, 1.2) ;
  gsl_test(!status, 
	   NAME(gsl_matrix) "_set traps 2nd index at upper bound") ;

  status = 0 ;
  x = FUNCTION(gsl_matrix,get)(m, i-1, 0) ;
  gsl_test(!status,
	   NAME(gsl_matrix) "_get traps 1st index below lower bound") ;
  gsl_test(x != 0, 
	   NAME(gsl_matrix) "_get returns zero for 1st index below lower bound") ;

  status = 0 ;
  x = FUNCTION(gsl_matrix,get)(m, 0, j-1) ;
  gsl_test(!status,
	   NAME(gsl_matrix) "_get traps 2nd index below lower bound") ;
  gsl_test(x != 0, 
	   NAME(gsl_matrix) "_get returns zero for 2nd index below lower bound") ;

  status = 0 ;
  x = FUNCTION(gsl_matrix,get)(m, N+1, 0) ;
  gsl_test(!status, 
	   NAME(gsl_matrix) "_get traps 1st index above upper bound") ;
  gsl_test(x != 0, 
	   NAME(gsl_matrix) "_get returns zero for 1st index above upper bound") ;

  status = 0 ;
  x = FUNCTION(gsl_matrix,get)(m, 0, M+1) ;
  gsl_test(!status, 
	   NAME(gsl_matrix) "_get traps 2nd index above upper bound") ;
  gsl_test(x != 0, 
	   NAME(gsl_matrix) "_get returns zero for 2nd index above upper bound") ;

  status = 0 ;
  x = FUNCTION(gsl_matrix,get)(m, N, 0) ;
  gsl_test(!status, 
	   NAME(gsl_matrix) "_get traps 1st index at upper bound") ;
  gsl_test(x != 0, 
	   NAME(gsl_matrix) "_get returns zero for 1st index at upper bound") ;

  status = 0 ;
  x = FUNCTION(gsl_matrix,get)(m, 0, M) ;
  gsl_test(!status, 
	   NAME(gsl_matrix) "_get traps 2nd index at upper bound") ;
  gsl_test(x != 0, 
	   NAME(gsl_matrix) "_get returns zero for 2nd index at upper bound") ;

}

