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

int filip_n = 82;
int filip_p = 11;

double filip_x[] = { -6.860120914, -4.324130045, -4.358625055,
-4.358426747, -6.955852379, -6.661145254, -6.355462942, -6.118102026,
-7.115148017, -6.815308569, -6.519993057, -6.204119983, -5.853871964,
-6.109523091, -5.79832982, -5.482672118, -5.171791386, -4.851705903,
-4.517126416, -4.143573228, -3.709075441, -3.499489089, -6.300769497,
-5.953504836, -5.642065153, -5.031376979, -4.680685696, -4.329846955,
-3.928486195, -8.56735134, -8.363211311, -8.107682739, -7.823908741,
-7.522878745, -7.218819279, -6.920818754, -6.628932138, -6.323946875,
-5.991399828, -8.781464495, -8.663140179, -8.473531488, -8.247337057,
-7.971428747, -7.676129393, -7.352812702, -7.072065318, -6.774174009,
-6.478861916, -6.159517513, -6.835647144, -6.53165267, -6.224098421,
-5.910094889, -5.598599459, -5.290645224, -4.974284616, -4.64454848,
-4.290560426, -3.885055584, -3.408378962, -3.13200249, -8.726767166,
-8.66695597, -8.511026475, -8.165388579, -7.886056648, -7.588043762,
-7.283412422, -6.995678626, -6.691862621, -6.392544977, -6.067374056,
-6.684029655, -6.378719832, -6.065855188, -5.752272167, -5.132414673,
-4.811352704, -4.098269308, -3.66174277, -3.2644011};

double filip_y[] = { 0.8116, 0.9072, 0.9052, 0.9039, 0.8053, 0.8377,
0.8667, 0.8809, 0.7975, 0.8162, 0.8515, 0.8766, 0.8885, 0.8859,
0.8959, 0.8913, 0.8959, 0.8971, 0.9021, 0.909, 0.9139, 0.9199, 0.8692,
0.8872, 0.89, 0.891, 0.8977, 0.9035, 0.9078, 0.7675, 0.7705, 0.7713,
0.7736, 0.7775, 0.7841, 0.7971, 0.8329, 0.8641, 0.8804, 0.7668,
0.7633, 0.7678, 0.7697, 0.77, 0.7749, 0.7796, 0.7897, 0.8131, 0.8498,
0.8741, 0.8061, 0.846, 0.8751, 0.8856, 0.8919, 0.8934, 0.894, 0.8957,
0.9047, 0.9129, 0.9209, 0.9219, 0.7739, 0.7681, 0.7665, 0.7703,
0.7702, 0.7761, 0.7809, 0.7961, 0.8253, 0.8602, 0.8809, 0.8301,
0.8664, 0.8834, 0.8898, 0.8964, 0.8963, 0.9074, 0.9119, 0.9228 } ;

int
main (void)
{
  size_t i, j;

  gsl_ieee_env_setup();

  {
    gsl_matrix X = gsl_matrix_view (longley_x, longley_n, longley_p);
    gsl_vector y = gsl_vector_view (longley_y, longley_n);
    gsl_vector * c = gsl_vector_alloc (longley_p);
    gsl_matrix * cov = gsl_matrix_alloc (longley_p, longley_p);
    gsl_vector diag;

    double chisq;

    double expected_c[7] = {  -3482258.63459582,
                              15.0618722713733,
                              -0.358191792925910E-01,
                              -2.02022980381683,
                              -1.03322686717359,
                              -0.511041056535807E-01,
                              1829.15146461355 };

    double expected_sd[7]  = {  890420.383607373,      
                                84.9149257747669,      
                                0.334910077722432E-01, 
                                0.488399681651699,     
                                0.214274163161675,     
                                0.226073200069370,     
                                455.478499142212 } ;  

    double expected_chisq = 836424.055505915;

    gsl_multifit_linear (&X, &y, c, cov, &chisq);

    gsl_test_rel (gsl_vector_get(c,0), expected_c[0], 1e-10, "longley gsl_fit_multilinear c0") ;
    gsl_test_rel (gsl_vector_get(c,1), expected_c[1], 1e-10, "longley gsl_fit_multilinear c1") ;
    gsl_test_rel (gsl_vector_get(c,2), expected_c[2], 1e-10, "longley gsl_fit_multilinear c2") ;
    gsl_test_rel (gsl_vector_get(c,3), expected_c[3], 1e-10, "longley gsl_fit_multilinear c3") ;
    gsl_test_rel (gsl_vector_get(c,4), expected_c[4], 1e-10, "longley gsl_fit_multilinear c4") ;
    gsl_test_rel (gsl_vector_get(c,5), expected_c[5], 1e-10, "longley gsl_fit_multilinear c5") ;
    gsl_test_rel (gsl_vector_get(c,6), expected_c[6], 1e-10, "longley gsl_fit_multilinear c6") ;

    diag = gsl_matrix_diagonal (cov);

    gsl_test_rel (gsl_vector_get(&diag,0), pow(expected_sd[0],2.0), 1e-10, "longley gsl_fit_multilinear cov00") ;
    gsl_test_rel (gsl_vector_get(&diag,1), pow(expected_sd[1],2.0), 1e-10, "longley gsl_fit_multilinear cov11") ;
    gsl_test_rel (gsl_vector_get(&diag,2), pow(expected_sd[2],2.0), 1e-10, "longley gsl_fit_multilinear cov22") ;
    gsl_test_rel (gsl_vector_get(&diag,3), pow(expected_sd[3],2.0), 1e-10, "longley gsl_fit_multilinear cov33") ;
    gsl_test_rel (gsl_vector_get(&diag,4), pow(expected_sd[4],2.0), 1e-10, "longley gsl_fit_multilinear cov44") ;
    gsl_test_rel (gsl_vector_get(&diag,5), pow(expected_sd[5],2.0), 1e-10, "longley gsl_fit_multilinear cov55") ;
    gsl_test_rel (gsl_vector_get(&diag,6), pow(expected_sd[6],2.0), 1e-10, "longley gsl_fit_multilinear cov66") ;

    gsl_test_rel (chisq, expected_chisq, 1e-10, "longley gsl_fit_multilinear chisq") ;

  }


  {
    gsl_matrix X = gsl_matrix_view (longley_x, longley_n, longley_p);
    gsl_vector y = gsl_vector_view (longley_y, longley_n);
    gsl_vector * w = gsl_vector_alloc (longley_n);
    gsl_vector * c = gsl_vector_alloc (longley_p);
    gsl_matrix * cov = gsl_matrix_alloc (longley_p, longley_p);
    gsl_vector diag;

    double chisq;

    double expected_c[7] = {  -3482258.63459582,
                              15.0618722713733,
                              -0.358191792925910E-01,
                              -2.02022980381683,
                              -1.03322686717359,
                              -0.511041056535807E-01,
                              1829.15146461355 };

    double expected_sd[7]  = {  890420.383607373,      
                                84.9149257747669,      
                                0.334910077722432E-01, 
                                0.488399681651699,     
                                0.214274163161675,     
                                0.226073200069370,     
                                455.478499142212 } ;  


    double expected_cov[7][7] = { 8531122.56783558, -166.727799925578,
0.261873708176346, 3.91188317230983, 1.1285582054705,
-0.889550869422687, -4362.58709870581, -166.727799925578,
0.0775861253030891, -1.98725210399982e-05, -0.000247667096727256,
-6.82911920718824e-05, 0.000136160797527761, 0.0775255245956248,
0.261873708176346, -1.98725210399982e-05, 1.20690316701888e-08,
1.66429546772984e-07, 3.61843600487847e-08, -6.78805814483582e-08,
-0.00013158719037715, 3.91188317230983, -0.000247667096727256,
1.66429546772984e-07, 2.56665052544717e-06, 6.96541409215597e-07,
-9.00858307771567e-07, -0.00197260370663974, 1.1285582054705,
-6.82911920718824e-05, 3.61843600487847e-08, 6.96541409215597e-07,
4.94032602583969e-07, -9.8469143760973e-08, -0.000576921112208274,
-0.889550869422687, 0.000136160797527761, -6.78805814483582e-08,
-9.00858307771567e-07, -9.8469143760973e-08, 5.49938542664952e-07,
0.000430074434198215, -4362.58709870581, 0.0775255245956248,
-0.00013158719037715, -0.00197260370663974, -0.000576921112208274,
0.000430074434198215, 2.23229587481535 } ;

    double expected_chisq = 836424.055505915;

    gsl_vector_set_all (w, 1.0);

    gsl_multifit_wlinear (&X, w, &y, c, cov, &chisq);

    gsl_test_rel (gsl_vector_get(c,0), expected_c[0], 1e-10, "longley gsl_fit_wmultilinear c0") ;
    gsl_test_rel (gsl_vector_get(c,1), expected_c[1], 1e-10, "longley gsl_fit_wmultilinear c1") ;
    gsl_test_rel (gsl_vector_get(c,2), expected_c[2], 1e-10, "longley gsl_fit_wmultilinear c2") ;
    gsl_test_rel (gsl_vector_get(c,3), expected_c[3], 1e-10, "longley gsl_fit_wmultilinear c3") ;
    gsl_test_rel (gsl_vector_get(c,4), expected_c[4], 1e-10, "longley gsl_fit_wmultilinear c4") ;
    gsl_test_rel (gsl_vector_get(c,5), expected_c[5], 1e-10, "longley gsl_fit_wmultilinear c5") ;
    gsl_test_rel (gsl_vector_get(c,6), expected_c[6], 1e-10, "longley gsl_fit_wmultilinear c6") ;

    for (i = 0; i < longley_p; i++) 
      {
        for (j = 0; j < longley_p; j++)
          {
            gsl_test_rel (gsl_matrix_get(cov,i,j), expected_cov[i][j], 1e-7, 
                          "longley gsl_fit_wmultilinear cov(%d,%d)", i, j) ;
          }
      }

    gsl_test_rel (chisq, expected_chisq, 1e-10, "longley gsl_fit_wmultilinear chisq") ;

  }



  {
    gsl_matrix * X = gsl_matrix_alloc (filip_n, filip_p);
    gsl_vector y = gsl_vector_view (filip_y, filip_n);
    gsl_vector * c = gsl_vector_alloc (filip_p);
    gsl_matrix * cov = gsl_matrix_alloc (filip_p, filip_p);
    gsl_vector diag;

    double chisq;

    double expected_c[11] = { -1467.48961422980,      
                              -2772.17959193342,      
                              -2316.37108160893,      
                              -1127.97394098372,      
                              -354.478233703349,      
                              -75.1242017393757,      
                              -10.8753180355343,      
                              -1.06221498588947,      
                              -0.670191154593408E-01, 
                              -0.246781078275479E-02, 
                              -0.402962525080404E-04 };

    double expected_sd[11]  = { 298.084530995537,     
                               559.779865474950,     
                               466.477572127796,     
                               227.204274477751,     
                               71.6478660875927,     
                               15.2897178747400,     
                               2.23691159816033,     
                               0.221624321934227,    
                               0.142363763154724E-01,
                               0.535617408889821E-03,
                               0.896632837373868E-05 };

    double expected_chisq = 0.795851382172941E-03;

    for (i = 0 ; i < filip_n; i++) 
      {
        for (j = 0; j < filip_p; j++) 
          {
            gsl_matrix_set(X, i, j, pow(filip_x[i], j));
          }
      }

    gsl_multifit_linear (X, &y, c, cov, &chisq);

    gsl_test_rel (gsl_vector_get(c,0), expected_c[0], 1e-7, "filip gsl_fit_multilinear c0") ;
    gsl_test_rel (gsl_vector_get(c,1), expected_c[1], 1e-7, "filip gsl_fit_multilinear c1") ;
    gsl_test_rel (gsl_vector_get(c,2), expected_c[2], 1e-7, "filip gsl_fit_multilinear c2") ;
    gsl_test_rel (gsl_vector_get(c,3), expected_c[3], 1e-7, "filip gsl_fit_multilinear c3") ;
    gsl_test_rel (gsl_vector_get(c,4), expected_c[4], 1e-7, "filip gsl_fit_multilinear c4") ;
    gsl_test_rel (gsl_vector_get(c,5), expected_c[5], 1e-7, "filip gsl_fit_multilinear c5") ;
    gsl_test_rel (gsl_vector_get(c,6), expected_c[6], 1e-7, "filip gsl_fit_multilinear c6") ;
    gsl_test_rel (gsl_vector_get(c,7), expected_c[7], 1e-7, "filip gsl_fit_multilinear c7") ;
    gsl_test_rel (gsl_vector_get(c,8), expected_c[8], 1e-7, "filip gsl_fit_multilinear c8") ;
    gsl_test_rel (gsl_vector_get(c,9), expected_c[9], 1e-7, "filip gsl_fit_multilinear c9") ;
    gsl_test_rel (gsl_vector_get(c,10), expected_c[10], 1e-7, "filip gsl_fit_multilinear c10") ;

    diag = gsl_matrix_diagonal (cov);

    gsl_test_rel (gsl_vector_get(&diag,0), pow(expected_sd[0],2.0), 1e-7, "filip gsl_fit_multilinear cov00") ;
    gsl_test_rel (gsl_vector_get(&diag,1), pow(expected_sd[1],2.0), 1e-7, "filip gsl_fit_multilinear cov11") ;
    gsl_test_rel (gsl_vector_get(&diag,2), pow(expected_sd[2],2.0), 1e-7, "filip gsl_fit_multilinear cov22") ;
    gsl_test_rel (gsl_vector_get(&diag,3), pow(expected_sd[3],2.0), 1e-7, "filip gsl_fit_multilinear cov33") ;
    gsl_test_rel (gsl_vector_get(&diag,4), pow(expected_sd[4],2.0), 1e-7, "filip gsl_fit_multilinear cov44") ;
    gsl_test_rel (gsl_vector_get(&diag,5), pow(expected_sd[5],2.0), 1e-7, "filip gsl_fit_multilinear cov55") ;
    gsl_test_rel (gsl_vector_get(&diag,6), pow(expected_sd[6],2.0), 1e-7, "filip gsl_fit_multilinear cov66") ;

    gsl_test_rel (chisq, expected_chisq, 1e-7, "filip gsl_fit_multilinear chisq") ;

  }



  /* now summarize the results */

  return gsl_test_summary ();
}
