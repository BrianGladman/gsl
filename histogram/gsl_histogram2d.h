#ifndef GSL_HISTOGRAM2D_H 
#define GSL_HISTOGRAM2D_H 

#include <stdlib.h>
#include <stdio.h>

typedef struct {
  size_t nx, ny ;
  double * xrange ;
  double * yrange ;
  double * bin ;
} gsl_histogram2d ;

typedef struct {
  size_t nx, ny ;
  double * xrange ;
  double * yrange ;
  double * sum ;
} gsl_histogram2d_pdf ;

gsl_histogram2d * gsl_histogram2d_calloc (size_t nx, size_t ny);
gsl_histogram2d * gsl_histogram2d_calloc_uniform (size_t nx, size_t ny,
					     double xmin, double xmax,
					     double ymin, double ymax);

void gsl_histogram2d_free (gsl_histogram2d * h);

int gsl_histogram2d_increment (gsl_histogram2d * h, double x, double y);
int gsl_histogram2d_accumulate (gsl_histogram2d * h, 
				double x, double y, double weight);
int gsl_histogram2d_find (const gsl_histogram2d * h, 
			  double x, double y, size_t * i, size_t * j);
int gsl_histogram2d_find_impl (const gsl_histogram2d * h, 
			       double x, double y, size_t * i, size_t * j);

double gsl_histogram2d_get (const gsl_histogram2d * h, size_t i, size_t j);
double gsl_histogram2d_get_xlowerlimit (const gsl_histogram2d * h, size_t i);
double gsl_histogram2d_get_xupperlimit (const gsl_histogram2d * h, size_t i);
double gsl_histogram2d_get_ylowerlimit (const gsl_histogram2d * h, size_t j);
double gsl_histogram2d_get_yupperlimit (const gsl_histogram2d * h, size_t j);
				     
double gsl_histogram2d_xmax (const gsl_histogram2d * h);
double gsl_histogram2d_xmin (const gsl_histogram2d * h);
size_t gsl_histogram2d_nx (const gsl_histogram2d * h);

double gsl_histogram2d_ymax (const gsl_histogram2d * h);
double gsl_histogram2d_ymin (const gsl_histogram2d * h);
size_t gsl_histogram2d_ny (const gsl_histogram2d * h);

void gsl_histogram2d_reset (gsl_histogram2d * h);

int gsl_histogram2d_fread (FILE * stream, gsl_histogram2d * h);
int gsl_histogram2d_fwrite (FILE * stream, const gsl_histogram2d * h) ;
int gsl_histogram2d_fprintf (FILE * stream, gsl_histogram2d * h, 
			   const char * format);
int gsl_histogram2d_fscanf (FILE * stream, gsl_histogram2d * h);

gsl_histogram2d_pdf * gsl_histogram2d_pdf_alloc (const gsl_histogram2d * h);
void gsl_histogram2d_pdf_free (gsl_histogram2d_pdf * p);
int gsl_histogram2d_pdf_sample (const gsl_histogram2d_pdf * p, 
				   double r1, double r2, 
				   double * x, double * y);

#endif /* GSL_HISTOGRAM2D_H */

