/* These tests are based on the NIST Statistical Reference Datasets
   See http://www.nist.gov/itl/div898/strd/index.html for more
   information. */

#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_fit.h>

#include <gsl/gsl_ieee_utils.h>

int
main (void)
{
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
    double expected_cov00 = 0.232818234301152;
    double expected_cov01 = 0.0;
    double expected_cov11 = 0.429796848199937E-03;
    double expected_sumsq = 26.6173985294224;
    
    gsl_fit_linear (x, xstride, y, ystride, norris_n, 
                    &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
    
    /* gsl_fit_wlinear (x, xstride, w, wstride, y, ystride, norris_n, 
                     &c0, &c1, &cov00, &cov01, &cov11, &sumsq); */
  
    gsl_test_rel (c0, expected_c0, 1e-10, "norris gsl_fit_linear c0") ;
    gsl_test_rel (c1, expected_c1, 1e-10, "norris gsl_fit_linear c1") ;
    gsl_test_rel (sqrt(cov00), expected_cov00, 1e-10, "norris gsl_fit_linear cov00") ;
    gsl_test_rel (cov01, expected_cov00, 1e-10, "norris gsl_fit_linear cov01") ;
    gsl_test_rel (sqrt(cov11), expected_cov11, 1e-10, "norris gsl_fit_linear cov11") ;
    gsl_test_rel (sumsq, expected_sumsq, 1e-10, "norris gsl_fit_linear sumsq") ;

  printf("c0 = %.18e\n", c0);
  printf("c1 = %.18e\n", c1);

  printf("cov00 = %.18e   s0 = %.18e\n", cov00, sqrt(cov00));
  printf("cov01 = %.18e\n", cov01);
  printf("cov11 = %.18e   s1 = %.18e\n", cov11, sqrt(cov11));

  printf("sumsq = %.18e\n", sumsq);


  }


  /* now summarize the results */

  return gsl_test_summary ();
}
