#include <config.h>
#include <math.h>
#include <gsl/gsl_ntuple.h>
#include <gsl/gsl_histogram.h>

struct data
{
  double x;
  double y;
  double z;
};

int sel_func (void *ntuple_data, void *params);
double val_func (void *ntuple_data, void *params);

int
main (void)
{
  struct data ntuple_row;
  int i;

  gsl_ntuple *ntuple = gsl_ntuple_open ("test.dat", &ntuple_row,
                                        sizeof (ntuple_row));

  gsl_histogram *h = gsl_histogram_calloc_uniform (100, 0., 10.);

  gsl_ntuple_select_fn S;
  gsl_ntuple_value_fn V;

  double scale = 1.5;

  S.function = &sel_func;
  S.params = &scale;

  V.function = &val_func;
  V.params = 0;

  gsl_ntuple_project (h, ntuple, &V, &S);

  gsl_histogram_fprintf (stdout, h, "%f", "%f");

  gsl_histogram_free (h);

  gsl_ntuple_close (ntuple);
}

int
sel_func (void *ntuple_data, void *params)
{
  double x, y, z, E, scale;
  scale = *(double *) params;

  x = ((struct data *) ntuple_data)->x;
  y = ((struct data *) ntuple_data)->y;
  z = ((struct data *) ntuple_data)->z;

  E = x * x + y * y + z * z;

  return E / scale > 1;
}

double
val_func (void *ntuple_data, void *params)
{
  double x, y, z;

  x = ((struct data *) ntuple_data)->x;
  y = ((struct data *) ntuple_data)->y;
  z = ((struct data *) ntuple_data)->z;

  return x * x + y * y + z * z;
}
