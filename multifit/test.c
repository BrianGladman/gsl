/* These tests are based on the NIST Statistical Reference Datasets
   See http://www.nist.gov/itl/div898/strd/index.html for more
   information. */

#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_multifit.h>

#include <gsl/gsl_ieee_utils.h>

int longley_n = 16;
int longley_p = 7;

double longley_x [] = {
  1,  83.0,   234289,   2356,     1590,    107608,  1947,
  1,  88.5,   259426,   2325,     1456,    108632,  1948,
  1,  88.2,   258054,   3682,     1616,    109773,  1949,
  1,  89.5,   284599,   3351,     1650,    110929,  1950,
  1,  96.2,   328975,   2099,     3099,    112075,  1951,
  1,  98.1,   346999,   1932,     3594,    113270,  1952,
  1,  99.0,   365385,   1870,     3547,    115094,  1953,
  1, 100.0,   363112,   3578,     3350,    116219,  1954,
  1, 101.2,   397469,   2904,     3048,    117388,  1955,
  1, 104.6,   419180,   2822,     2857,    118734,  1956,
  1, 108.4,   442769,   2936,     2798,    120445,  1957,
  1, 110.8,   444546,   4681,     2637,    121950,  1958,
  1, 112.6,   482704,   3813,     2552,    123366,  1959,
  1, 114.2,   502601,   3931,     2514,    125368,  1960,
  1, 115.7,   518173,   4806,     2572,    127852,  1961,
  1, 116.9,   554894,   4007,     2827,    130081,  1962 } ;

double longley_y[] = {60323, 61122, 60171, 61187, 63221, 63639, 64989, 63761,
                       66019, 67857, 68169, 66513, 68655, 69564, 69331, 70551};

int
main (void)
{
  size_t i;

  gsl_ieee_env_setup();

  {
    gsl_matrix X = gsl_matrix_view (longley_x, longley_n, longley_p);
    gsl_vector y = gsl_vector_view (longley_y, longley_n);
    gsl_vector * w = gsl_vector_alloc (longley_n);
    gsl_vector * c = gsl_vector_alloc (longley_p);
    gsl_matrix * cov = gsl_matrix_alloc (longley_p, longley_p);
    double chisq;

    double expected_c[7] = {  -3482258.63459582,
                              15.0618722713733,
                              -0.358191792925910E-01,
                              -2.02022980381683,
                              -1.03322686717359,
                              -0.511041056535807E-01,
                              1829.15146461355 };

    gsl_vector_set_all (w, 1.0);

    gsl_multifit_wlinear (&X, w, &y, c, cov, &chisq);

    gsl_test_rel (gsl_vector_get(c,0), expected_c[0], 1e-10, "longley gsl_fit_wmultilinear c0") ;
    gsl_test_rel (gsl_vector_get(c,1), expected_c[1], 1e-10, "longley gsl_fit_wmultilinear c1") ;
    gsl_test_rel (gsl_vector_get(c,2), expected_c[2], 1e-10, "longley gsl_fit_wmultilinear c2") ;
    gsl_test_rel (gsl_vector_get(c,3), expected_c[3], 1e-10, "longley gsl_fit_wmultilinear c3") ;
    gsl_test_rel (gsl_vector_get(c,4), expected_c[4], 1e-10, "longley gsl_fit_wmultilinear c4") ;
    gsl_test_rel (gsl_vector_get(c,5), expected_c[5], 1e-10, "longley gsl_fit_wmultilinear c5") ;
    gsl_test_rel (gsl_vector_get(c,6), expected_c[6], 1e-10, "longley gsl_fit_wmultilinear c6") ;


    gsl_matrix_fprintf (stdout, cov, "%g");

  }

  /* now summarize the results */

  return gsl_test_summary ();
}
