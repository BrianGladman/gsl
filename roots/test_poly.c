#include <math.h>
#include <float.h>

#include <gsl_errno.h>
#include <gsl_math.h>
#include <gsl_roots.h>
#include <gsl_test.h>

void 
test_poly(void)
{
  int n = 0 ;
  double x[2], y[3] ;
  gsl_complex z[2], Z[3] ;

  /* Quadratic */
  
  n = gsl_root_solve_quadratic (4.0, -20.0, 26.0, x) ;
  
  gsl_test(n != 0, "gsl_root_solve_quadratic, no roots, (2x - 5)^2 = -1") ;
  
  n = gsl_root_solve_quadratic (4.0, -20.0, 25.0, x) ;
  
  gsl_test(n != 1, "gsl_root_solve_quadratic, one root, (2x - 5)^2 = 0");
  gsl_test_rel(x[0], 2.5, 1e-9, "x1, (2x - 5)^2 = 0");
  
  n = gsl_root_solve_quadratic (4.0, -20.0, 21.0, x);
  
  gsl_test(n != 2, "gsl_root_solve_quadratic, two roots, (2x - 5)^2 = 4") ;
  gsl_test_rel(x[0], 1.5, 1e-9, "x1, (2x - 5)^2 = 4");
  gsl_test_rel(x[1], 3.5, 1e-9, "x2, (2x - 5)^2 = 4");
  
  n = gsl_root_solve_quadratic (4.0, 7.0, 0.0, x);
  
  gsl_test(n != 2, "gsl_root_solve_quadratic, two roots, x(4x + 7) = 0") ;
  gsl_test_rel(x[0], -1.75, 1e-9, "x1, x(4x + 7) = 0");
  gsl_test_rel(x[1], 0.0, 1e-9, "x2, x(4x + 7) = 0");
  
  n = gsl_root_solve_quadratic (5.0, 0.0, -20.0, x);
  
  gsl_test(n != 2, "gsl_root_solve_quadratic, two roots b = 0, 5 x^2 = 20") ;
  gsl_test_rel(x[0], -2.0, 1e-9, "x1, 5 x^2 = 20");
  gsl_test_rel(x[1], 2.0, 1e-9, "x2, 5 x^2 = 20");

  /* Cubic */

  n = gsl_root_solve_cubic (0.0, 0.0, -27.0, y) ;
  
  gsl_test(n != 1, "gsl_root_solve_cubic, one root, x^3 = 27");
  gsl_test_rel(y[0], 3.0, 1e-9, "x1, x^3 = 27");
  
  n = gsl_root_solve_cubic (-51.0, 867.0, -4913.0, y);
  
  gsl_test(n != 3, "gsl_root_solve_cubic, three roots, (x-17)^3=0") ;
  gsl_test_rel(y[0], 17.0, 1e-9, "x1, (x-17)^3=0");
  gsl_test_rel(y[1], 17.0, 1e-9, "x2, (x-17)^3=0");
  gsl_test_rel(y[2], 17.0, 1e-9, "x3, (x-17)^3=0");
  
  n = gsl_root_solve_cubic (-57.0, 1071.0, -6647.0, y);
  
  gsl_test(n != 3, "gsl_root_solve_cubic, three roots, (x-17)(x-17)(x-23)=0") ;
  gsl_test_rel(y[0], 17.0, 1e-9, "x1, (x-17)(x-17)(x-23)=0");
  gsl_test_rel(y[1], 17.0, 1e-9, "x2, (x-17)(x-17)(x-23)=0");
  gsl_test_rel(y[2], 23.0, 1e-9, "x3, (x-17)(x-17)(x-23)=0");
  
  n = gsl_root_solve_cubic (-143.0, 5087.0, -50065.0, y);
  
  gsl_test(n != 3, "gsl_root_solve_cubic, three roots, (x-17)(x-31)(x-95)=0") ;
  gsl_test_rel(y[0], 17.0, 1e-9, "x1, (x-17)(x-31)(x-95)=0");
  gsl_test_rel(y[1], 31.0, 1e-9, "x2, (x-17)(x-31)(x-95)=0");
  gsl_test_rel(y[2], 95.0, 1e-9, "x3, (x-17)(x-31)(x-95)=0");
  
  /* Quadratic with complex roots */

  n = gsl_root_complex_solve_quadratic (4.0, -20.0, 26.0, z) ;
  
  gsl_test(n != 2, "gsl_root_complex_solve_quadratic, 2 roots (2x - 5)^2 = -1");
  gsl_test_rel(GSL_REAL(z[0]), 2.5, 1e-9, "z1.real, (2x - 5)^2 = -1");
  gsl_test_rel(GSL_IMAG(z[0]), -0.5, 1e-9, "z1.imag, (2x - 5)^2 = -1");
  
  gsl_test_rel(GSL_REAL(z[1]), 2.5, 1e-9, "z2.real, (2x - 5)^2 = -1");
  gsl_test_rel(GSL_IMAG(z[1]), 0.5, 1e-9, "z2.imag, (2x - 5)^2 = -1");
  
  n = gsl_root_complex_solve_quadratic (4.0, -20.0, 25.0, z) ;
  
  gsl_test(n != 1, "gsl_root_complex_solve_quadratic, one root, (2x - 5)^2 = 0");
  gsl_test_rel(GSL_REAL(z[0]), 2.5, 1e-9, "z1.real, (2x - 5)^2 = 0");
  gsl_test_rel(GSL_IMAG(z[0]), 0.0, 1e-9, "z1.imag (2x - 5)^2 = 0");
  
  n = gsl_root_complex_solve_quadratic (4.0, -20.0, 21.0, z);
  
  gsl_test(n != 2, "gsl_root_complex_solve_quadratic, two roots, (2x - 5)^2 = 4") ;
  gsl_test_rel(GSL_REAL(z[0]), 1.5, 1e-9, "z1.real, (2x - 5)^2 = 4");
  gsl_test_rel(GSL_IMAG(z[0]), 0.0, 1e-9, "z1.imag, (2x - 5)^2 = 4");
  gsl_test_rel(GSL_REAL(z[1]), 3.5, 1e-9, "z2.real, (2x - 5)^2 = 4");
  gsl_test_rel(GSL_IMAG(z[1]), 0.0, 1e-9, "z2.imag, (2x - 5)^2 = 4");
  
  n = gsl_root_complex_solve_quadratic (4.0, 7.0, 0.0, z);
  
  gsl_test(n != 2, "gsl_root_complex_solve_quadratic, two roots, x(4x + 7) = 0") ;
  gsl_test_rel(GSL_REAL(z[0]), -1.75, 1e-9, "z1.real, x(4x + 7) = 0");
  gsl_test_rel(GSL_IMAG(z[0]), 0.0, 1e-9, "z1.imag, x(4x + 7) = 0");
  gsl_test_rel(GSL_REAL(z[1]), 0.0, 1e-9, "z2.real, x(4x + 7) = 0");
  gsl_test_rel(GSL_IMAG(z[1]), 0.0, 1e-9, "z2.imag, x(4x + 7) = 0");
  
  n = gsl_root_complex_solve_quadratic (5.0, 0.0, -20.0, z);
  
  gsl_test(n != 2, "gsl_root_complex_solve_quadratic, two roots b = 0, 5 x^2 = 20") ;
  gsl_test_rel(GSL_REAL(z[0]), -2.0, 1e-9, "z1.real, 5 x^2 = 20");
  gsl_test_rel(GSL_IMAG(z[0]), 0.0, 1e-9, "z1.imag, 5 x^2 = 20");
  gsl_test_rel(GSL_REAL(z[1]), 2.0, 1e-9, "z2.real, 5 x^2 = 20");
  gsl_test_rel(GSL_IMAG(z[1]), 0.0, 1e-9, "z2.imag, 5 x^2 = 20");
  
  n = gsl_root_complex_solve_quadratic (5.0, 0.0, 20.0, z);
  
  gsl_test(n != 2, "gsl_root_complex_solve_quadratic, two roots b = 0, 5 x^2 = -20") ;
  gsl_test_rel(GSL_REAL(z[0]), 0.0, 1e-9, "z1.real, 5 x^2 = -20");
  gsl_test_rel(GSL_IMAG(z[0]), -2.0, 1e-9, "z1.imag, 5 x^2 = -20");
  gsl_test_rel(GSL_REAL(z[1]), 0.0, 1e-9, "z2.real, 5 x^2 = -20");
  gsl_test_rel(GSL_IMAG(z[1]), 2.0, 1e-9, "z2.imag, 5 x^2 = -20");

  /* Cubic with complex roots */
  
  n = gsl_root_complex_solve_cubic (0.0, 0.0, -27.0, Z) ;
  
  gsl_test(n != 3, "gsl_root_complex_solve_cubic, three root, x^3 = 27");
  gsl_test_rel(GSL_REAL(Z[0]), -1.5, 1e-9, "z1.real, x^3 = 27");
  gsl_test_rel(GSL_IMAG(Z[0]), -1.5 * sqrt(3.0), 1e-9, "z1.imag, x^3 = 27");
  gsl_test_rel(GSL_REAL(Z[1]), -1.5, 1e-9, "z2.real, x^3 = 27");
  gsl_test_rel(GSL_IMAG(Z[1]), 1.5 * sqrt(3.0), 1e-9, "z2.imag, x^3 = 27");
  gsl_test_rel(GSL_REAL(Z[2]), 3.0, 1e-9, "z3.real, x^3 = 27");
  gsl_test_rel(GSL_IMAG(Z[2]), 0.0, 1e-9, "z3.imag, x^3 = 27");
  
  n = gsl_root_complex_solve_cubic (-51.0, 867.0, -4913.0, Z);
  
  gsl_test(n != 3, "gsl_root_complex_solve_cubic, three roots, (x-17)^3=0") ;
  gsl_test_rel(GSL_REAL(Z[0]), 17.0, 1e-9, "z1.real, (x-17)^3=0");
  gsl_test_rel(GSL_IMAG(Z[0]), 0.0, 1e-9, "z1.imag, (x-17)^3=0");
  gsl_test_rel(GSL_REAL(Z[1]), 17.0, 1e-9, "z2.real, (x-17)^3=0");
  gsl_test_rel(GSL_IMAG(Z[1]), 0.0, 1e-9, "z2.imag, (x-17)^3=0");
  gsl_test_rel(GSL_REAL(Z[2]), 17.0, 1e-9, "z3.real, (x-17)^3=0");
  gsl_test_rel(GSL_IMAG(Z[2]), 0.0, 1e-9, "z3.imag, (x-17)^3=0");
  
  n = gsl_root_complex_solve_cubic (-57.0, 1071.0, -6647.0, Z);
  
  gsl_test(n != 3, "gsl_root_complex_solve_cubic, three roots, (x-17)(x-17)(x-23)=0") ;
  gsl_test_rel(GSL_REAL(Z[0]), 17.0, 1e-9, "z1.real, (x-17)(x-17)(x-23)=0");
  gsl_test_rel(GSL_IMAG(Z[0]), 0.0, 1e-9, "z1.imag, (x-17)(x-17)(x-23)=0");
  gsl_test_rel(GSL_REAL(Z[1]), 17.0, 1e-9, "z2.real, (x-17)(x-17)(x-23)=0");
  gsl_test_rel(GSL_IMAG(Z[1]), 0.0, 1e-9, "z2.imag, (x-17)(x-17)(x-23)=0");
  gsl_test_rel(GSL_REAL(Z[2]), 23.0, 1e-9, "z3.real, (x-17)(x-17)(x-23)=0");
  gsl_test_rel(GSL_IMAG(Z[2]), 0.0, 1e-9, "z3.imag, (x-17)(x-17)(x-23)=0");
  
  
  n = gsl_root_complex_solve_cubic (-143.0, 5087.0, -50065.0, Z);
  
  gsl_test(n != 3, "gsl_root_complex_solve_cubic, three roots, (x-17)(x-31)(x-95)=0") ;
  gsl_test_rel(GSL_REAL(Z[0]), 17.0, 1e-9, "z1.real, (x-17)(x-31)(x-95)=0");
  gsl_test_rel(GSL_IMAG(Z[0]), 0.0, 1e-9, "z1.imag, (x-17)(x-31)(x-95)=0");
  gsl_test_rel(GSL_REAL(Z[1]), 31.0, 1e-9, "z2.real, (x-17)(x-31)(x-95)=0");
  gsl_test_rel(GSL_IMAG(Z[1]), 0.0, 1e-9, "z2.imag, (x-17)(x-31)(x-95)=0");
  gsl_test_rel(GSL_REAL(Z[2]), 95.0, 1e-9, "z3.real, (x-17)(x-31)(x-95)=0");
  gsl_test_rel(GSL_IMAG(Z[2]), 0.0, 1e-9, "z3.imag, (x-17)(x-31)(x-95)=0");

}



