/* These tests are based on the NIST Statistical Reference Datasets
   See http://www.nist.gov/itl/div898/strd/index.html for more
   information. */

#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_fit.h>

#include <gsl/gsl_ieee_utils.h>

int norris_n = 36;

double norris_x[] = { 0.2, 337.4, 118.2, 884.6, 10.1, 226.5, 666.3, 996.3,
                      448.6, 777.0, 558.2, 0.4, 0.6, 775.5, 666.9, 338.0, 
                      447.5, 11.6, 556.0, 228.1, 995.8, 887.6, 120.2, 0.3, 
                      0.3, 556.8, 339.1, 887.2, 999.0, 779.0, 11.1, 118.3,
                      229.2, 669.1, 448.9, 0.5 } ;

double norris_y[] = { 0.1, 338.8, 118.1, 888.0, 9.2, 228.1, 668.5, 998.5,
                      449.1, 778.9, 559.2, 0.3, 0.1, 778.1, 668.8, 339.3, 
                      448.9, 10.8, 557.7, 228.3, 998.0, 888.8, 119.6, 0.3, 
                      0.6, 557.6, 339.3, 888.0, 998.5, 778.9, 10.2, 117.6,
                      228.9, 668.4, 449.2, 0.2};

int noint1_n = 11;
double noint1_x[] = { 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70 };
double noint1_y[] = { 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140};

int noint2_n = 3;
double noint2_x[] = { 4, 5, 6 } ;
double noint2_y[] = { 3, 4, 4 } ;

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


  double x[1000], y[1000], w[1000];

  size_t xstride = 2, wstride = 3, ystride = 5;
  size_t i;

  for (i = 0; i < norris_n; i++) 
    {
      x[i*xstride] = norris_x[i];
      w[i*wstride] = 1.0;
      y[i*ystride] = norris_y[i];
    }

  gsl_ieee_env_setup();

  {
    double c0, c1, cov00, cov01, cov11, sumsq;
       
    double expected_c0 = -0.262323073774029;
    double expected_c1 =  1.00211681802045; 
    double expected_cov00 = pow(0.232818234301152, 2.0);
    double expected_cov01 = -7.74327536339570e-05;  /* computed from octave */
    double expected_cov11 = pow(0.429796848199937E-03, 2.0);
    double expected_sumsq = 26.6173985294224;
    
    gsl_fit_linear (x, xstride, y, ystride, norris_n, 
                    &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
    
    /* gsl_fit_wlinear (x, xstride, w, wstride, y, ystride, norris_n, 
                     &c0, &c1, &cov00, &cov01, &cov11, &sumsq); */
  
    gsl_test_rel (c0, expected_c0, 1e-10, "norris gsl_fit_linear c0") ;
    gsl_test_rel (c1, expected_c1, 1e-10, "norris gsl_fit_linear c1") ;
    gsl_test_rel (cov00, expected_cov00, 1e-10, "norris gsl_fit_linear cov00") ;
    gsl_test_rel (cov01, expected_cov01, 1e-10, "norris gsl_fit_linear cov01") ;
    gsl_test_rel (cov11, expected_cov11, 1e-10, "norris gsl_fit_linear cov11") ;
    gsl_test_rel (sumsq, expected_sumsq, 1e-10, "norris gsl_fit_linear sumsq") ;
  }

  {
    double c0, c1, cov00, cov01, cov11, sumsq;
       
    double expected_c0 = -0.262323073774029;
    double expected_c1 =  1.00211681802045; 
    double expected_cov00 = 6.92384428759429e-02;  /* computed from octave */
    double expected_cov01 = -9.89095016390515e-05; /* computed from octave */
    double expected_cov11 = 2.35960747164148e-07;  /* computed from octave */
    double expected_sumsq = 26.6173985294224;
    
    gsl_fit_wlinear (x, xstride, w, wstride, y, ystride, norris_n, 
                     &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
  
    gsl_test_rel (c0, expected_c0, 1e-10, "norris gsl_fit_wlinear c0") ;
    gsl_test_rel (c1, expected_c1, 1e-10, "norris gsl_fit_wlinear c1") ;
    gsl_test_rel (cov00, expected_cov00, 1e-10, "norris gsl_fit_wlinear cov00") ;
    gsl_test_rel (cov01, expected_cov01, 1e-10, "norris gsl_fit_wlinear cov01") ;
    gsl_test_rel (cov11, expected_cov11, 1e-10, "norris gsl_fit_wlinear cov11") ;
    gsl_test_rel (sumsq, expected_sumsq, 1e-10, "norris gsl_fit_wlinear sumsq") ;
  }

  for (i = 0; i < noint1_n; i++) 
    {
      x[i*xstride] = noint1_x[i];
      w[i*wstride] = 1.0;
      y[i*ystride] = noint1_y[i];
    }

  {
    double c0, c1, cov00, cov01, cov11, sumsq;
       
    double expected_c1 = 2.07438016528926; 
    double expected_cov11 = pow(0.165289256198347E-01, 2.0);  
    double expected_sumsq = 127.272727272727;
    
    gsl_fit_mul (x, xstride, y, ystride, noint1_n, &c1, &cov11, &sumsq);
  
    gsl_test_rel (c1, expected_c1, 1e-10, "noint1 gsl_fit_mul c1") ;
    gsl_test_rel (cov11, expected_cov11, 1e-10, "noint1 gsl_fit_mul cov11") ;
    gsl_test_rel (sumsq, expected_sumsq, 1e-10, "noint1 gsl_fit_mul sumsq") ;
  }

  {
    double c0, c1, cov00, cov01, cov11, sumsq;
       
    double expected_c1 = 2.07438016528926; 
    double expected_cov11 = 2.14661371686165e-05; /* computed from octave */
    double expected_sumsq = 127.272727272727;
    
    gsl_fit_wmul (x, xstride, w, wstride, y, ystride, noint1_n, &c1, &cov11, &sumsq);

    gsl_test_rel (c1, expected_c1, 1e-10, "noint1 gsl_fit_wmul c1") ;
    gsl_test_rel (cov11, expected_cov11, 1e-10, "noint1 gsl_fit_wmul cov11") ;
    gsl_test_rel (sumsq, expected_sumsq, 1e-10, "noint1 gsl_fit_wmul sumsq") ;
  }


  for (i = 0; i < noint2_n; i++) 
    {
      x[i*xstride] = noint2_x[i];
      w[i*wstride] = 1.0;
      y[i*ystride] = noint2_y[i];
    }

  {
    double c0, c1, cov00, cov01, cov11, sumsq;
       
    double expected_c1 = 0.727272727272727; 
    double expected_cov11 = pow(0.420827318078432E-01, 2.0);  
    double expected_sumsq = 0.272727272727273;
    
    gsl_fit_mul (x, xstride, y, ystride, noint2_n, &c1, &cov11, &sumsq);
  
    gsl_test_rel (c1, expected_c1, 1e-10, "noint2 gsl_fit_mul c1") ;
    gsl_test_rel (cov11, expected_cov11, 1e-10, "noint2 gsl_fit_mul cov11") ;
    gsl_test_rel (sumsq, expected_sumsq, 1e-10, "noint2 gsl_fit_mul sumsq") ;
  }

  {
    double c0, c1, cov00, cov01, cov11, sumsq;
       
    double expected_c1 = 0.727272727272727; 
    double expected_cov11 = 1.29870129870130e-02 ; /* computed from octave */
    double expected_sumsq = 0.272727272727273;
    
    gsl_fit_wmul (x, xstride, w, wstride, y, ystride, noint2_n, &c1, &cov11, &sumsq);

    gsl_test_rel (c1, expected_c1, 1e-10, "noint2 gsl_fit_wmul c1") ;
    gsl_test_rel (cov11, expected_cov11, 1e-10, "noint2 gsl_fit_wmul cov11") ;
    gsl_test_rel (sumsq, expected_sumsq, 1e-10, "noint2 gsl_fit_wmul sumsq") ;
  }


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

    gsl_fit_wmultilinear (&X, w, &y, c, cov, &chisq);

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
