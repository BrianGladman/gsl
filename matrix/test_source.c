int main (void) 
{
  gsl_matrix * m;
  size_t i, j;
  int k = 0;

  m = gsl_matrix_alloc(N,M) ;

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
	if(m->data[i*M + j] != k) 
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
	if(FUNCTION(gsl_matrix,get)(m,i,j) != k) 
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

  {
    FILE * f = fopen("test.txt","w") ;
    k = 0 ;
    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	k++ ; FUNCTION(gsl_matrix,set)(m,i,j,(double)k) ;
      }
    }
    
    FUNCTION(gsl_matrix,fprintf)(f, m, "%.4f") ;

    fclose(f) ;
  }


  {
    int status = 0 ;
    FILE * f = fopen("test.txt","r") ;
    BASE * mm = FUNCTION(gsl_matrix,calloc)(N,M) ;

    FUNCTION(gsl_matrix,fscanf)(f, mm) ;
    k = 0 ;
    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	k++ ; 
	if (mm->data[i*M + j] != k)
	  status = 1 ;
      }
    }
    
    gsl_test(status, NAME(gsl_matrix) "_fprintf and fscanf work correctly") ;
    
    fclose (f) ;
    gsl_matrix_free (mm) ;
  }


  {
    FILE * f = fopen("test.dat", "w") ;
    k = 0 ;
    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	k++ ; FUNCTION(gsl_matrix,set)(m,i,j,(double)k) ;
      }
    }

    FUNCTION(gsl_matrix,fwrite)(f, m) ;
    fclose (f) ;
  }

  {
    int status = 0 ;
    FILE * f = fopen("test.dat","r") ;
    BASE * mm = FUNCTION(gsl_matrix,calloc)(N,M) ;

    FUNCTION(gsl_matrix,fread)(f, mm) ;
    k = 0 ;
    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	k++ ; 
	if (mm->data[i*M + j] != k)
	    status = 1 ;
      }
    }

    gsl_test(status, NAME(gsl_matrix) "_write and read work correctly") ;

    fclose (f) ;
    FUNCTION(gsl_matrix,free) (mm) ;
  }

  return gsl_test_summary ();
}
