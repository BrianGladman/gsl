#include <stdlib.h>

#include <gsl_errno.h>

typedef struct
{
  size_t n;
  double * data ;
} gsl_vector ;

double
gsl_vector_get(gsl_vector v, size_t i)
{
#ifdef GSL_RANGE_CHECK
  if (i < 0 || i >= v.n)   /* if size_t is unsigned then i<0 is impossible! */
    {
      abort() ;
    }
#endif
  return v.data[i] ;
}

void
gsl_vector_set(gsl_vector * v, size_t i, double x)
{
#ifdef GSL_RANGE_CHECK
  if (i < 0 || i >= v.n) 
    {
      abort() ;
    }
#endif
  v.data[i] = x ;
}

gsl_vector * 
gsl_vector_alloc (size_t n)
{
  gsl_vector * v ;

  if (n == 0)
    {
      GSL_ERROR ("vector length n must be positive integer", GSL_EDOM);
    }
  
  
  v = (gsl_vector *) malloc(sizeof(gsl_vector)) ;

  if (v == 0) 
    {
      GSL_ERROR ("failed to allocate space for vector", GSL_ENOMEM);
    }
  else 
    
  v->data = malloc(n * sizeof(double)) ;

  if (v->data == 0) 
    {
      GSL_ERROR ("failed to allocate space for vector", GSL_ENOMEM);
    }
  
  v->n = n ;

  return v ;
}

int
gsl_vector_calloc (gsl_vector * v, size_t n)
{
  size_t i ;

  int status = gsl_vector_alloc (v, n) ;
  
  if (status) 
    return status ;

  for (i = 0 ; i < n; i++)  /* initialize to zero */
    {
      v->data[i] = 0.0 ;
    }

  return 0 ;
}


int
gsl_vector_free (gsl_vector * v)
{
  free(v->data) ;
  v->n = 0 ;
  return 0 ;
}


/* Maximum and minimum elements of a vector.
 */
void vector_min_max_i(const int    * v, unsigned long size, int    * min, int    * max);
void vector_min_max_l(const long   * v, unsigned long size, long   * min, long   * max);
void vector_min_max_f(const float  * v, unsigned long size, float  * min, float  * max);
void vector_min_max_d(const double * v, unsigned long size, double * min, double * max);


/* Multiply vector by a scalar.
 */
void vector_multiply_scalar_i(int    * v, unsigned long size, int    s);
void vector_multiply_scalar_l(long   * v, unsigned long size, long   s);
void vector_multiply_scalar_f(float  * v, unsigned long size, float  s);
void vector_multiply_scalar_d(double * v, unsigned long size, double s);


/* Inner product.
 */
int    vector_dot_i(const int    * v1, const int    * v2, unsigned long size);
long   vector_dot_l(const long   * v1, const long   * v2, unsigned long size);
float  vector_dot_f(const float  * v1, const float  * v2, unsigned long size);
double vector_dot_d(const double * v1, const double * v2, unsigned long size);


/* Euclidean length of a vector.
 */
float  vector_euclen_f(const float  * v, unsigned long size);
double vector_euclen_d(const double * v, unsigned long size);


/* Distance between two vectors, sup-norm.
 */
int    vector_dist_sup_i(const int    * v1, const int    * v2, unsigned long size);
long   vector_dist_sup_l(const long   * v1, const long   * v2, unsigned long size);
float  vector_dist_sup_f(const float  * v1, const float  * v2, unsigned long size);
double vector_dist_sup_d(const double * v1, const double * v2, unsigned long size);


/* Distance between two vectors, l_2-norm.
 */
float  vector_dist_l2_f(const float  * v1, const float  * v2, unsigned long size);
double vector_dist_l2_d(const double * v1, const double * v2, unsigned long size);


/* Elements below the given threshold in absolute value,
 * relative to the largest element, are set to zero.
 */
void vector_lower_cut_f(float  * v, unsigned long size, float  cut);
void vector_lower_cut_d(double * v, unsigned long size, double cut);


/* Dump a vector to a stream using a given separator string.
 */
void vector_dump_i(FILE *, const int    * v, unsigned long size, const char * format, const char * sep);
void vector_dump_l(FILE *, const long   * v, unsigned long size, const char * format, const char * sep);
void vector_dump_f(FILE *, const float  * v, unsigned long size, const char * format, const char * sep);
void vector_dump_d(FILE *, const double * v, unsigned long size, const char * format, const char * sep);
void vector_dump_c(FILE *, const char   * v, unsigned long size, const char * format, const char * sep);



enum jpeg_scale  {jpegLinearScale,  jpegLogScale};
enum jpeg_scheme {jpegSignedScheme, jpegAbsScheme};

/* Create a jpeg image from a vector, using a given scan line length.
 * Values are quantized in the given range (lo, hi), and
 * values outside this range are color clipped; lo and hi
 * should be _positive_ numbers.
 *     scale  = jpegLinearScale  | jpegLogScale
 *     scheme = jpegSignedScheme | jpegAbsScheme
 * Available only for libraries compiled with the jpeg option.
 */
void vector_jpeg_d(FILE *fd, const double * v, unsigned long size,
		   int scanlength,
		   double lo, double hi,
		   enum jpeg_scale  scale,
		   enum jpeg_scheme scheme
		   );
		   
