#include <math.h>
#include <float.h>

#include <gsl_errno.h>
#include <gsl_math.h>
#include <gsl_poly.h>
#include <gsl_test.h>

void 
test_poly(void)
{
  int n = 0 ;
  double x0, x1, y0, y1, y2 ;
  gsl_complex z0, z1, Z0, Z1, Z2 ;

  /* Quadratic */
  
  n = gsl_poly_solve_quadratic (4.0, -20.0, 26.0, &x0, &x1) ;
  
  gsl_test(n != 0, "gsl_poly_solve_quadratic, no roots, (2x - 5)^2 = -1") ;
  
  n = gsl_poly_solve_quadratic (4.0, -20.0, 25.0, &x0, &x1) ;
  
  gsl_test(n != 1, "gsl_poly_solve_quadratic, one root, (2x - 5)^2 = 0");
  gsl_test_rel(x0, 2.5, 1e-9, "x1, (2x - 5)^2 = 0");
  
  n = gsl_poly_solve_quadratic (4.0, -20.0, 21.0, &x0, &x1);
  
  gsl_test(n != 2, "gsl_poly_solve_quadratic, two roots, (2x - 5)^2 = 4") ;
  gsl_test_rel(x0, 1.5, 1e-9, "x1, (2x - 5)^2 = 4");
  gsl_test_rel(x1, 3.5, 1e-9, "x2, (2x - 5)^2 = 4");
  
  n = gsl_poly_solve_quadratic (4.0, 7.0, 0.0, &x0, &x1);
  
  gsl_test(n != 2, "gsl_poly_solve_quadratic, two roots, x(4x + 7) = 0") ;
  gsl_test_rel(x0, -1.75, 1e-9, "x1, x(4x + 7) = 0");
  gsl_test_rel(x1, 0.0, 1e-9, "x2, x(4x + 7) = 0");
  
  n = gsl_poly_solve_quadratic (5.0, 0.0, -20.0, &x0, &x1);
  
  gsl_test(n != 2, "gsl_poly_solve_quadratic, two roots b = 0, 5 x^2 = 20") ;
  gsl_test_rel(x0, -2.0, 1e-9, "x1, 5 x^2 = 20");
  gsl_test_rel(x1, 2.0, 1e-9, "x2, 5 x^2 = 20");

  /* Cubic */

  n = gsl_poly_solve_cubic (0.0, 0.0, -27.0, &y0, &y1, &y2) ;
  
  gsl_test(n != 1, "gsl_poly_solve_cubic, one root, x^3 = 27");
  gsl_test_rel(y0, 3.0, 1e-9, "x1, x^3 = 27");
  
  n = gsl_poly_solve_cubic (-51.0, 867.0, -4913.0, &y0, &y1, &y2);
  
  gsl_test(n != 3, "gsl_poly_solve_cubic, three roots, (x-17)^3=0") ;
  gsl_test_rel(y0, 17.0, 1e-9, "x1, (x-17)^3=0");
  gsl_test_rel(y1, 17.0, 1e-9, "x2, (x-17)^3=0");
  gsl_test_rel(y2, 17.0, 1e-9, "x3, (x-17)^3=0");
  
  n = gsl_poly_solve_cubic (-57.0, 1071.0, -6647.0, &y0, &y1, &y2);
  
  gsl_test(n != 3, "gsl_poly_solve_cubic, three roots, (x-17)(x-17)(x-23)=0") ;
  gsl_test_rel(y0, 17.0, 1e-9, "x1, (x-17)(x-17)(x-23)=0");
  gsl_test_rel(y1, 17.0, 1e-9, "x2, (x-17)(x-17)(x-23)=0");
  gsl_test_rel(y2, 23.0, 1e-9, "x3, (x-17)(x-17)(x-23)=0");
  
  n = gsl_poly_solve_cubic (-143.0, 5087.0, -50065.0, &y0, &y1, &y2);
  
  gsl_test(n != 3, "gsl_poly_solve_cubic, three roots, (x-17)(x-31)(x-95)=0") ;
  gsl_test_rel(y0, 17.0, 1e-9, "x1, (x-17)(x-31)(x-95)=0");
  gsl_test_rel(y1, 31.0, 1e-9, "x2, (x-17)(x-31)(x-95)=0");
  gsl_test_rel(y2, 95.0, 1e-9, "x3, (x-17)(x-31)(x-95)=0");
  
  /* Quadratic with complex roots */

  n = gsl_poly_complex_solve_quadratic (4.0, -20.0, 26.0, &z0, &z1) ;
  
  gsl_test(n != 2, "gsl_poly_complex_solve_quadratic, 2 roots (2x - 5)^2 = -1");
  gsl_test_rel(GSL_REAL(z0), 2.5, 1e-9, "z1.real, (2x - 5)^2 = -1");
  gsl_test_rel(GSL_IMAG(z0), -0.5, 1e-9, "z1.imag, (2x - 5)^2 = -1");
  
  gsl_test_rel(GSL_REAL(z1), 2.5, 1e-9, "z2.real, (2x - 5)^2 = -1");
  gsl_test_rel(GSL_IMAG(z1), 0.5, 1e-9, "z2.imag, (2x - 5)^2 = -1");
  
  n = gsl_poly_complex_solve_quadratic (4.0, -20.0, 25.0, &z0, &z1) ;
  
  gsl_test(n != 1, "gsl_poly_complex_solve_quadratic, one root, (2x - 5)^2 = 0");
  gsl_test_rel(GSL_REAL(z0), 2.5, 1e-9, "z1.real, (2x - 5)^2 = 0");
  gsl_test_rel(GSL_IMAG(z0), 0.0, 1e-9, "z1.imag (2x - 5)^2 = 0");
  
  n = gsl_poly_complex_solve_quadratic (4.0, -20.0, 21.0, &z0, &z1);
  
  gsl_test(n != 2, "gsl_poly_complex_solve_quadratic, two roots, (2x - 5)^2 = 4") ;
  gsl_test_rel(GSL_REAL(z0), 1.5, 1e-9, "z1.real, (2x - 5)^2 = 4");
  gsl_test_rel(GSL_IMAG(z0), 0.0, 1e-9, "z1.imag, (2x - 5)^2 = 4");
  gsl_test_rel(GSL_REAL(z1), 3.5, 1e-9, "z2.real, (2x - 5)^2 = 4");
  gsl_test_rel(GSL_IMAG(z1), 0.0, 1e-9, "z2.imag, (2x - 5)^2 = 4");
  
  n = gsl_poly_complex_solve_quadratic (4.0, 7.0, 0.0, &z0, &z1);
  
  gsl_test(n != 2, "gsl_poly_complex_solve_quadratic, two roots, x(4x + 7) = 0") ;
  gsl_test_rel(GSL_REAL(z0), -1.75, 1e-9, "z1.real, x(4x + 7) = 0");
  gsl_test_rel(GSL_IMAG(z0), 0.0, 1e-9, "z1.imag, x(4x + 7) = 0");
  gsl_test_rel(GSL_REAL(z1), 0.0, 1e-9, "z2.real, x(4x + 7) = 0");
  gsl_test_rel(GSL_IMAG(z1), 0.0, 1e-9, "z2.imag, x(4x + 7) = 0");
  
  n = gsl_poly_complex_solve_quadratic (5.0, 0.0, -20.0, &z0, &z1);
  
  gsl_test(n != 2, "gsl_poly_complex_solve_quadratic, two roots b = 0, 5 x^2 = 20") ;
  gsl_test_rel(GSL_REAL(z0), -2.0, 1e-9, "z1.real, 5 x^2 = 20");
  gsl_test_rel(GSL_IMAG(z0), 0.0, 1e-9, "z1.imag, 5 x^2 = 20");
  gsl_test_rel(GSL_REAL(z1), 2.0, 1e-9, "z2.real, 5 x^2 = 20");
  gsl_test_rel(GSL_IMAG(z1), 0.0, 1e-9, "z2.imag, 5 x^2 = 20");
  
  n = gsl_poly_complex_solve_quadratic (5.0, 0.0, 20.0, &z0, &z1);
  
  gsl_test(n != 2, "gsl_poly_complex_solve_quadratic, two roots b = 0, 5 x^2 = -20") ;
  gsl_test_rel(GSL_REAL(z0), 0.0, 1e-9, "z1.real, 5 x^2 = -20");
  gsl_test_rel(GSL_IMAG(z0), -2.0, 1e-9, "z1.imag, 5 x^2 = -20");
  gsl_test_rel(GSL_REAL(z1), 0.0, 1e-9, "z2.real, 5 x^2 = -20");
  gsl_test_rel(GSL_IMAG(z1), 2.0, 1e-9, "z2.imag, 5 x^2 = -20");

  /* Cubic with complex roots */
  
  n = gsl_poly_complex_solve_cubic (0.0, 0.0, -27.0, &Z0, &Z1, &Z2) ;
  
  gsl_test(n != 3, "gsl_poly_complex_solve_cubic, three root, x^3 = 27");
  gsl_test_rel(GSL_REAL(Z0), -1.5, 1e-9, "z1.real, x^3 = 27");
  gsl_test_rel(GSL_IMAG(Z0), -1.5 * sqrt(3.0), 1e-9, "z1.imag, x^3 = 27");
  gsl_test_rel(GSL_REAL(Z1), -1.5, 1e-9, "z2.real, x^3 = 27");
  gsl_test_rel(GSL_IMAG(Z1), 1.5 * sqrt(3.0), 1e-9, "z2.imag, x^3 = 27");
  gsl_test_rel(GSL_REAL(Z2), 3.0, 1e-9, "z3.real, x^3 = 27");
  gsl_test_rel(GSL_IMAG(Z2), 0.0, 1e-9, "z3.imag, x^3 = 27");
  
  n = gsl_poly_complex_solve_cubic (-51.0, 867.0, -4913.0, &Z0, &Z1, &Z2);
  
  gsl_test(n != 3, "gsl_poly_complex_solve_cubic, three roots, (x-17)^3=0") ;
  gsl_test_rel(GSL_REAL(Z0), 17.0, 1e-9, "z1.real, (x-17)^3=0");
  gsl_test_rel(GSL_IMAG(Z0), 0.0, 1e-9, "z1.imag, (x-17)^3=0");
  gsl_test_rel(GSL_REAL(Z1), 17.0, 1e-9, "z2.real, (x-17)^3=0");
  gsl_test_rel(GSL_IMAG(Z1), 0.0, 1e-9, "z2.imag, (x-17)^3=0");
  gsl_test_rel(GSL_REAL(Z2), 17.0, 1e-9, "z3.real, (x-17)^3=0");
  gsl_test_rel(GSL_IMAG(Z2), 0.0, 1e-9, "z3.imag, (x-17)^3=0");
  
  n = gsl_poly_complex_solve_cubic (-57.0, 1071.0, -6647.0, &Z0, &Z1, &Z2);
  
  gsl_test(n != 3, "gsl_poly_complex_solve_cubic, three roots, (x-17)(x-17)(x-23)=0") ;
  gsl_test_rel(GSL_REAL(Z0), 17.0, 1e-9, "z1.real, (x-17)(x-17)(x-23)=0");
  gsl_test_rel(GSL_IMAG(Z0), 0.0, 1e-9, "z1.imag, (x-17)(x-17)(x-23)=0");
  gsl_test_rel(GSL_REAL(Z1), 17.0, 1e-9, "z2.real, (x-17)(x-17)(x-23)=0");
  gsl_test_rel(GSL_IMAG(Z1), 0.0, 1e-9, "z2.imag, (x-17)(x-17)(x-23)=0");
  gsl_test_rel(GSL_REAL(Z2), 23.0, 1e-9, "z3.real, (x-17)(x-17)(x-23)=0");
  gsl_test_rel(GSL_IMAG(Z2), 0.0, 1e-9, "z3.imag, (x-17)(x-17)(x-23)=0");
  
  
  n = gsl_poly_complex_solve_cubic (-143.0, 5087.0, -50065.0, &Z0, &Z1, &Z2);
  
  gsl_test(n != 3, "gsl_poly_complex_solve_cubic, three roots, (x-17)(x-31)(x-95)=0") ;
  gsl_test_rel(GSL_REAL(Z0), 17.0, 1e-9, "z1.real, (x-17)(x-31)(x-95)=0");
  gsl_test_rel(GSL_IMAG(Z0), 0.0, 1e-9, "z1.imag, (x-17)(x-31)(x-95)=0");
  gsl_test_rel(GSL_REAL(Z1), 31.0, 1e-9, "z2.real, (x-17)(x-31)(x-95)=0");
  gsl_test_rel(GSL_IMAG(Z1), 0.0, 1e-9, "z2.imag, (x-17)(x-31)(x-95)=0");
  gsl_test_rel(GSL_REAL(Z2), 95.0, 1e-9, "z3.real, (x-17)(x-31)(x-95)=0");
  gsl_test_rel(GSL_IMAG(Z2), 0.0, 1e-9, "z3.imag, (x-17)(x-31)(x-95)=0");

}



