/* Created: [GJ] Sun Apr 20 05:12:22 EDT 1997
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "constants.h"
#include "sorting.h"
#include "vectors.h"


/* yibdee bibdee pazulaboo... */
#define GEN_ALLOC(FUNC, TYPE)                                         \
TYPE * FUNC(unsigned long size)                                       \
{                                                                     \
  TYPE * p = (TYPE *) malloc((size_t)(size * sizeof(TYPE)));          \
  if(p == 0) {                                                        \
    char buff[100];                                                   \
    sprintf(buff,"FUNC: alloc of size %ld failed", size);             \
    push_error(buff, Error_Alloc_);                                   \
  }                                                                   \
  return p;                                                           \
}                                                                     \

GEN_ALLOC(new_vector_i, int)
GEN_ALLOC(new_vector_l, long)
GEN_ALLOC(new_vector_d, double)
GEN_ALLOC(new_vector_f, float)
GEN_ALLOC(new_vector_c, char)

GEN_ALLOC(new_vector_ip, int *)
GEN_ALLOC(new_vector_lp, long *)
GEN_ALLOC(new_vector_dp, double *)
GEN_ALLOC(new_vector_fp, float *)
GEN_ALLOC(new_vector_cp, char *)


void free_vector(void * p)
{
  if(p != 0) free(p);
}


/* yibdee bibdee pazulaboo indeed... */
#define GEN_MATALLOC(FUNC, PALC, ALC, TYPE)                           \
TYPE ** FUNC(unsigned long size)                                      \
{                                                                     \
  TYPE ** m = PALC(size);                                             \
  if(m == 0) {                                                        \
    char buff[100];                                                   \
    sprintf(buff,"FUNC: alloc of size %ld failed", size);             \
    push_error(buff, Error_Alloc_);                                   \
    return 0;                                                         \
  }                                                                   \
  else {                                                              \
    unsigned long i;                                                  \
    for(i=0; i<size; i++) {                                           \
      m[i] = ALC(size);                                               \
    }                                                                 \
    return m;                                                         \
  }                                                                   \
}                                                                     \

#define GEN_MATFREE(FUNC, TYPE)                                     \
void FUNC(TYPE ** m, unsigned long size)                            \
{                                                                   \
  unsigned long i;                                                  \
  if(m != 0) {                                                      \
    for(i=0; i<size; i++) {                                         \
      if(m[i] != 0) free_vector(m[i]);                              \
    }                                                               \
    free_vector(m);                                                 \
  }                                                                 \
}                                                                   \


GEN_MATALLOC(new_matrix_i, new_vector_ip, new_vector_i, int)
GEN_MATALLOC(new_matrix_l, new_vector_lp, new_vector_l, long)
GEN_MATALLOC(new_matrix_d, new_vector_dp, new_vector_d, double)
GEN_MATALLOC(new_matrix_f, new_vector_fp, new_vector_f, float)

GEN_MATFREE(free_matrix_i, int)
GEN_MATFREE(free_matrix_l, long)
GEN_MATFREE(free_matrix_d, double)
GEN_MATFREE(free_matrix_f, float)


#define GEN_COPYVEC(FUNC, ALLOC, TYPE)                              \
TYPE * FUNC(TYPE * v, unsigned long size)                           \
{                                                                   \
  if(v != 0) {                                                      \
    TYPE * nv = ALLOC(size);                                        \
    if(nv==0) return 0;                                             \
    memcpy(nv, v, size*sizeof(TYPE));                               \
    return nv;                                                      \
  }                                                                 \
  else {                                                            \
    return 0;                                                       \
  }                                                                 \
}                                                                   \

GEN_COPYVEC(copy_vector_i, new_vector_i, int)
GEN_COPYVEC(copy_vector_l, new_vector_l, long)
GEN_COPYVEC(copy_vector_d, new_vector_d, double)
GEN_COPYVEC(copy_vector_f, new_vector_f, float)
GEN_COPYVEC(copy_vector_c, new_vector_c, char)



#define GEN_COPYMAT(FUNC, ALLOC, COPYVEC, TYPE)                     \
TYPE ** FUNC(TYPE ** m, unsigned long size)                         \
{                                                                   \
  if(m != 0) {                                                      \
    unsigned long i;                                                \
    TYPE ** nm = ALLOC(size);                                       \
    if(nm==0) return 0;                                             \
    for(i=0; i<size; i++) {                                         \
      nm[i] = COPYVEC(m[i], size);                                  \
    }                                                               \
    return nm;                                                      \
  }                                                                 \
  else {                                                            \
    return 0;                                                       \
  }                                                                 \
}                                                                   \

GEN_COPYMAT(copy_matrix_i, new_matrix_i, copy_vector_i, int)
GEN_COPYMAT(copy_matrix_l, new_matrix_l, copy_vector_l, long)
GEN_COPYMAT(copy_matrix_d, new_matrix_d, copy_vector_d, double)
GEN_COPYMAT(copy_matrix_f, new_matrix_f, copy_vector_f, float)


#define GEN_MINMAX(FUNC, TYPE, INIT)                                       \
void FUNC(const TYPE * v, unsigned long size, TYPE * min, TYPE * max)      \
{                                                                          \
  unsigned long i;                                                         \
  TYPE min_t = INIT;                                                       \
  TYPE max_t = -INIT;                                                      \
  for(i=0; i<size; i++) {                                                  \
    min_t = Min(v[i], min_t);                                              \
    max_t = Max(v[i], max_t);                                              \
  }                                                                        \
  *min = min_t;                                                            \
  *max = max_t;                                                            \
}                                                                          \

GEN_MINMAX(vector_min_max_i, int,    max_int)
GEN_MINMAX(vector_min_max_l, long,   max_long)
GEN_MINMAX(vector_min_max_f, float,  max_float)
GEN_MINMAX(vector_min_max_d, double, max_double)



#define GEN_MULT_SCALAR(FUNC, TYPE)                                       \
void FUNC(TYPE * v, unsigned long size, TYPE s)                           \
{                                                                         \
  unsigned long i;                                                        \
  for(i=0; i<size; i++) v[i] *= s;                                        \
}                                                                         \

GEN_MULT_SCALAR(vector_multiply_scalar_i, int)
GEN_MULT_SCALAR(vector_multiply_scalar_l, long)
GEN_MULT_SCALAR(vector_multiply_scalar_f, float)
GEN_MULT_SCALAR(vector_multiply_scalar_d, double)


#define GEN_DOT(FUNC, TYPE)                                               \
TYPE FUNC(const TYPE * v1, const TYPE * v2, unsigned long size)           \
{                                                                         \
  unsigned long i;                                                        \
  TYPE result = 0.;                                                       \
  for(i=0; i<size; i++) result += v1[i]*v2[i];                            \
  return result;                                                          \
}                                                                         \

GEN_DOT(vector_dot_i, int)
GEN_DOT(vector_dot_l, long)
GEN_DOT(vector_dot_f, float)
GEN_DOT(vector_dot_d, double)


float vector_euclen_f(const float * v, unsigned long size)
{
  unsigned long i;
  float sumsq = 0.;
  for(i=0; i<size; i++) sumsq += v[i]*v[i];
  return sqrt(sumsq);
}

double vector_euclen_d(const double * v, unsigned long size)
{
  unsigned long i;
  double sumsq = 0.;
  for(i=0; i<size; i++) sumsq += v[i]*v[i];
  return sqrt(sumsq);
}


#define GEN_DIST_SUP(FUNC, TYPE, ABS)                                   \
TYPE FUNC(const TYPE * v1, const TYPE * v2, unsigned long size)         \
{                                                                       \
  unsigned long i;                                                      \
  TYPE dist = 0;                                                        \
  for(i=0; i<size; i++) {                                               \
    dist = Max(dist, ABS(v1[i]-v2[i]));                                 \
  }                                                                     \
  return dist;                                                          \
}                                                                       \

GEN_DIST_SUP(vector_dist_sup_i, int,    abs)
GEN_DIST_SUP(vector_dist_sup_l, long,   abs)
GEN_DIST_SUP(vector_dist_sup_f, float,  fabs)
GEN_DIST_SUP(vector_dist_sup_d, double, fabs)


float vector_dist_l2_f(const float * v1, const float * v2, unsigned long size)
{
  unsigned long i;
  float sumsq = 0.;
  for(i=0; i<size; i++) sumsq += (v1[i]-v2[i])*(v1[i]-v2[i]);
  return sqrt(sumsq);
}

double vector_dist_l2_d(const double * v1, const double * v2, unsigned long size)
{
  unsigned long i;
  double sumsq = 0.;
  for(i=0; i<size; i++) sumsq += (v1[i]-v2[i])*(v1[i]-v2[i]);
  return sqrt(sumsq);
}


#define GEN_LOWER_CUT(FUNC, TYPE, ABS, MINMAX)                          \
void FUNC(TYPE * v, unsigned long size, TYPE cut)                       \
{                                                                       \
  unsigned long i;                                                      \
  TYPE min, max, abs_max;                                               \
  MINMAX(v, size, &min, &max);                                          \
  abs_max = Max(ABS(min),ABS(max));                                     \
  for(i=0; i<size; i++) {                                               \
    if(ABS(v[i])/abs_max < cut) v[i] = (TYPE) 0;                        \
  }                                                                     \
}                                                                       \

GEN_LOWER_CUT(vector_lower_cut_f, float,  fabs, vector_min_max_f)
GEN_LOWER_CUT(vector_lower_cut_d, double, fabs, vector_min_max_d)


#define GEN_DUMP(FUNC, TYPE)                                            \
void FUNC(FILE * fp, const TYPE * v, unsigned long size,                \
		 const char * format, const char * sep)                 \
{                                                                       \
  unsigned long i;                                                      \
  for(i=0; i<size; i++) {                                               \
    fprintf(fp, format, v[i]);                                          \
    if(i < size-1) fprintf(fp, "%s", sep);                              \
  }                                                                     \
}                                                                       \

GEN_DUMP(vector_dump_i, int)
GEN_DUMP(vector_dump_l, long)
GEN_DUMP(vector_dump_f, float)
GEN_DUMP(vector_dump_d, double)
GEN_DUMP(vector_dump_c, char)


#ifdef JPEG_OPTION
#include <jpeglib.h>

static void get_rgb_color(double val,
			  unsigned char * r,
			  unsigned char * g,
			  unsigned char * b,
			  double lo, double hi,
			  enum jpeg_scale scale,
			  enum jpeg_scheme scheme)
{
  /* enforce positivity of thresholds */
  hi = fabs(hi);
  lo = fabs(lo);

  /* swap them if necessary */
  if(hi < lo) {
    double t = hi; hi = lo; lo = t;
  }

  if(hi == lo) {
    *r = 255; *g = 255; *b = 255;
  }
  else {
    double intensity;             /* intensity in (0,1) */
    double abs_val = fabs(val);

    abs_val = Min(abs_val, hi);   /* value clipping */
    abs_val = Max(abs_val, lo);

    if(scale == jpegLinearScale) {
      intensity = (abs_val - lo)/(hi-lo);
    }
    else if(scale == jpegLogScale) {
      lo += 2.e-111;
      hi += 2.e-111;
      intensity = (log(abs_val)-log(lo))/(log(hi)-log(lo));
      intensity = Max(intensity, 0.);
      intensity = Min(intensity, 1.);
    }

    if(scheme == jpegAbsScheme) {
      *r = Min(255, (int)(255.*intensity));
      *g = 0;
      *b = Min(255, (int)(255.*intensity));
    }
    else if(scheme == jpegSignedScheme) {
      *r = Min(255, (int)(255.*intensity));
      *g = Min(50,  (int)(50.*intensity));
      *b = 0;
      if(val < 0.) {
	unsigned char t = *r; *r = *b; *b = t; /* swap red and blue */
      }
    }
  }
}
			  

static void make_rgb_row_buffer(JSAMPLE * row_buff, int width,
				const double * v,
				double norm,
				double lo, double hi,
				enum jpeg_scale scale,
				enum jpeg_scheme scheme)
{
  int i_v;
  for(i_v=0; i_v<width; i_v++) {
    unsigned char r, g, b;
    get_rgb_color(v[i_v] / norm, &r, &g, &b, lo, hi, scale, scheme);
    row_buff[3*i_v]     = r;
    row_buff[3*i_v + 1] = g;
    row_buff[3*i_v + 2] = b;
  }
}


void vector_jpeg_d(FILE * fd, const double * v, unsigned long size,
		   int scanlength,
		   double lo, double hi,
		   enum jpeg_scale  scale,
		   enum jpeg_scheme scheme
		   )
{
  unsigned long i_input = 0;
  int image_width  = scanlength;
  int image_height = size / scanlength;
  int quality = 75;
  int row_stride;

  double min, max, abs_max;

  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  JSAMPROW row_pointer[1];

  vector_min_max_d(v, size, &min, &max);
  abs_max = Max(fabs(min), fabs(max));

  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo, fd);

  cinfo.image_width  = image_width;
  cinfo.image_height = image_height;
  cinfo.input_components = 3;
  cinfo.in_color_space = JCS_RGB;

  jpeg_set_defaults(&cinfo);
  jpeg_set_quality(&cinfo, quality, TRUE /* limit to baseline-JPEG values */);
  
  jpeg_start_compress(&cinfo, TRUE);

  row_stride = image_width * 3;	/* JSAMPLEs per row */
  row_pointer[0] = (JSAMPLE *) malloc((size_t) row_stride);
  while (cinfo.next_scanline < cinfo.image_height) {
    /* jpeg_write_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could pass
     * more than one scanline at a time if that's more convenient.
     */
    make_rgb_row_buffer(row_pointer[0], image_width, v + i_input,
			abs_max, lo, hi, scale, scheme
			);
    i_input += image_width;
    jpeg_write_scanlines(&cinfo, row_pointer, 1);
  }
  free(row_pointer[0]);

  jpeg_finish_compress(&cinfo);
  fflush(fd);

  jpeg_destroy_compress(&cinfo);
}


/*
void jpeg_test(void)
{
  int i, j;
  double data[512][512];

  for(i=0; i<512; i++) {
    for(j=0; j<512; j++) {
      data[i][j] = exp( -(fabs(i-256.) + fabs(j-256.))/512.);
      if(i < j) data[i][j] *= -1.;
    }
  }

  vector_jpeg_d(stdout, &(data[0][0]), 512*512, 512,
		3.e-1, 1.,
		jpegLogScale,
		jpegSignedScheme
		);
}
*/

#endif /* !JPEG_OPTION */
