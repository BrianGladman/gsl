#include <config.h>
#include <fcntl.h>
#include <stdio.h>
#include <gsl_histogram2d.h>
#include <gsl_test.h>

#define N 107
#define M 239
#define N1 17
#define M1 23

int main (void) 
{
  gsl_histogram2d * h;
  size_t i, j, k;
  
  h = gsl_histogram2d_alloc(N,M) ;
  
  gsl_test(h->xrange == 0, "gsl_histogram2d_alloc returns valid yrange pointer") ;
  gsl_test(h->yrange == 0, "gsl_histogram2d_alloc returns valid yrange pointer") ;
  gsl_test(h->bin == 0, "gsl_histogram2d_alloc returns valid bin pointer") ;
  
  gsl_test(h->nx != N, "gsl_histogram2d_alloc returns valid nx") ;
  gsl_test(h->ny != M, "gsl_histogram2d_alloc returns valid ny") ;
  
  k = 0 ;
  for (i = 0 ; i < N ; i++) 
    {
      for (j = 0; j < M; j++)
	{
	  k++ ;
	  gsl_histogram2d_accumulate (h, (double)i, (double)j, (double)k) ;
	} ;
    }

  { 
    int status = 0 ;
    k = 0 ;
    for (i = 0 ; i < N ; i++) 
      {
	for (j = 0 ; j < M ; j++) 
	  {
	    k++ ;
	    if(h->bin[i * M + j] != (double)k) 
	      {
		status = 1 ;
	      }
	  }
      }
    gsl_test(status,"gsl_histogram2d_accumulate writes into array correctly") ;
  }
  
  {
    int status = 0 ;
    k = 0 ;
    for (i = 0 ; i < N ; i++) {
      for (j = 0 ; j < M ; j++) {
	k++ ;
	if(gsl_histogram2d_get (h, i, j) != (double)k ) 
	  status = 1 ;
      } ;
    }
    gsl_test(status, "gsl_histogram2d_get reads from array correctly") ;
  }

  gsl_histogram2d_reset (h) ;
  
  { 
    int status = 0 ;

    for (i = 0 ; i < N*M ; i++) {
      if (h->bin[i] != 0)
	status = 1 ;
    }
    gsl_test(status, "gsl_histogram2d_reset zeros array correctly") ;
  }
  
  gsl_histogram2d_free (h) ;
  h = gsl_histogram2d_alloc(N1,M1) ;

  {

    int status = 0 ;

    for (i = 0 ; i < N1 ; i++) 
      {
	for (j = 0 ; j < M1 ; j++) 
	  {
	    gsl_histogram2d_increment (h, (double)i, (double)j) ;
	
	    for (k = 0 ; k <= i * M1 + j ; k++) {
	      if (h->bin[k] != 1) {
		status = 1;
	      }
	    }
	    
	    for (k = i * M1 + j + 1 ; k < N1 * M1 ; k++) {
	      if (h->bin[k] != 0) {
		status = 1;
	      }
	    }
	  }
      }
    gsl_test(status, "gsl_histogram2d_increment works correctly") ;
  }

  gsl_histogram2d_free (h) ;
  h = gsl_histogram2d_alloc(N,M) ;

  {
    int status = 0 ;
    for (i = 0; i < N; i++)
      {
	double x0 = gsl_histogram2d_get_xlowerlimit (h, i);
	double x1 = gsl_histogram2d_get_xupperlimit (h, i);

	if (x0 != i || x1 != i + 1)
	  {
	    status = 1 ;
	  }
      }
    gsl_test(status,"gsl_histogram2d_get_xlowerlimit and xupperlimit works correctly") ;
  }


  {
    int status = 0 ;
    for (i = 0; i < M; i++)
      {
	double y0 = gsl_histogram2d_get_ylowerlimit (h, i);
	double y1 = gsl_histogram2d_get_yupperlimit (h, i);

	if (y0 != i || y1 != i + 1)
	  {
	    status = 1 ;
	  }
      }
    gsl_test(status,"gsl_histogram2d_get_ylowerlimit and yupperlimit works correctly") ;
  }


  {
    int status = 0 ;
    if (gsl_histogram2d_xmax (h) != N) 
      status = 1 ;
    gsl_test(status, "gsl_histogram2d_xmax works correctly") ;
  }

  {
    int status = 0 ;
    if (gsl_histogram2d_xmin (h) != 0) 
      status = 1 ;
    gsl_test(status, "gsl_histogram2d_xmin works correctly") ;
  }

  {
    int status = 0 ;
    if (gsl_histogram2d_nx (h) != N) 
      status = 1 ;
    gsl_test(status, "gsl_histogram2d_nx works correctly") ;
  }

  {
    int status = 0 ;
    if (gsl_histogram2d_ymax (h) != M) 
      status = 1 ;
    gsl_test(status, "gsl_histogram2d_ymax works correctly") ;
  }

  {
    int status = 0 ;
    if (gsl_histogram2d_ymin (h) != 0) 
      status = 1 ;
    gsl_test(status, "gsl_histogram2d_ymin works correctly") ;
  }

  {
    int status = 0 ;
    if (gsl_histogram2d_ny (h) != M) 
      status = 1 ;
    gsl_test(status, "gsl_histogram2d_ny works correctly") ;
  }


  gsl_histogram2d_free(h);  /* free whatever is in h */

  h = gsl_histogram2d_alloc_uniform (N1, M1, 0.0, 5.0, 0.0, 5.0) ;

  gsl_test(h->xrange == 0,
	   "gsl_histogram2d_alloc_uniform returns valid range pointer") ;
  gsl_test(h->yrange == 0,
	   "gsl_histogram2d_alloc_uniform returns valid range pointer") ;
  gsl_test(h->bin == 0,
	   "gsl_histogram2d_alloc_uniform returns valid bin pointer") ;
  gsl_test(h->nx != N1,
	   "gsl_histogram2d_alloc_uniform returns valid nx") ;
  gsl_test(h->ny != M1,
	   "gsl_histogram2d_alloc_uniform returns valid ny") ;

  gsl_histogram2d_accumulate (h, 0.0, 3.01, 1.0) ;
  gsl_histogram2d_accumulate (h, 0.1, 2.01, 2.0) ;
  gsl_histogram2d_accumulate (h, 0.2, 1.01, 3.0) ;
  gsl_histogram2d_accumulate (h, 0.3, 0.01, 4.0) ;

  { 
    size_t i1, i2, i3, i4;
    size_t j1, j2, j3, j4;
    double expected ;
    int status ; 
    status = gsl_histogram2d_find (h, 0.0, 3.01, &i1, &j1) ;
    status = gsl_histogram2d_find (h, 0.1, 2.01, &i2, &j2) ;
    status = gsl_histogram2d_find (h, 0.2, 1.01, &i3, &j3) ;
    status = gsl_histogram2d_find (h, 0.3, 0.01, &i4, &j4) ;

    for (i = 0; i < N1; i++)
      {
	for (j = 0; j < M1; j++)
	  {
	    if (i == i1 && j == j1) {
	      expected = 1.0 ;
	    } else if (i == i2 && j == j2) {
	      expected = 2.0 ;
	    } else if (i == i3 && j == j3) {
	      expected = 3.0 ;
	    } else if (i == i4 && j == j4) {
	      expected = 4.0 ;
	    } else {
	      expected = 0.0 ;
	    }

	    if (h->bin[i * M1 + j] != expected) {
	      status = 1 ;
	    }
	  }
      }
    gsl_test(status, "gsl_histogram2d_find works correctly") ;
  }

  {
    FILE * f = fopen("test.dat","w") ;
    gsl_histogram2d_fprintf(f, h, "%.18g") ;
    fclose(f) ;
  }

  {
    FILE * f = fopen("test.dat","r") ;
    gsl_histogram2d * hh = gsl_histogram2d_alloc (N1,M1);
    int status = 0 ;

    gsl_histogram2d_fscanf(f, hh) ;

    for (i = 0; i <= N1; i++) 
      {
	if (h->xrange[i] != hh->xrange[i]) 
	  {
	    printf("xrange[%d] : %g orig vs %g\n", i, h->xrange[i],hh->xrange[i]) ;
	    status = 1 ;
	  }
      }

    for (j = 0; j <= M1; j++) 
      {
	if (h->yrange[j] != hh->yrange[j]) 
	  {
	    printf("yrange[%d] : %g orig vs %g\n", j, h->yrange[j],hh->yrange[j]) ;
	    status = 1 ;
	  }
      }

    for (i = 0; i < N1*M1; i++) 
      {
	if (h->bin[i] != hh->bin[i]) 
	  {
	    printf("bin[%d] : %g orig vs %g\n", i, h->bin[i],hh->bin[i]) ;
	    status = 1 ;
	  }
      }

    gsl_test(status, "gsl_histogram2d_fprintf and fscanf work correctly") ;

    gsl_histogram2d_free (hh);
    fclose(f) ;
  }

  {
    FILE * f = fopen("test.dat","w") ;
    gsl_histogram2d_fwrite (f, h) ;
    fclose(f) ;
  }

  {
    FILE * f = fopen("test.dat","r") ;
    gsl_histogram2d * hh = gsl_histogram2d_alloc (N1,M1);
    int status = 0 ;

    gsl_histogram2d_fread(f, hh) ;

    for (i = 0; i <= N1; i++) 
      {
	if (h->xrange[i] != hh->xrange[i]) 
	  {
	    printf("xrange[%d] : %g orig vs %g\n", i, h->xrange[i],hh->xrange[i]) ;
	    status = 1 ;
	  }
      }

    for (j = 0; j <= M1; j++) 
      {
	if (h->yrange[j] != hh->yrange[j]) 
	  {
	    printf("yrange[%d] : %g orig vs %g\n", j, h->yrange[j],hh->yrange[j]) ;
	    status = 1 ;
	  }
      }

    for (i = 0; i < N1*M1; i++) 
      {
	if (h->bin[i] != hh->bin[i]) 
	  {
	    printf("bin[%d] : %g orig vs %g\n", i, h->bin[i],hh->bin[i]) ;
	    status = 1 ;
	  }
      }

    gsl_test(status, "gsl_histogram2d_fwrite and fread work correctly") ;

    gsl_histogram2d_free (hh);
    fclose(f) ;
  }


  return gsl_test_summary ();
}


