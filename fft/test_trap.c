/* fft/test_trap.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

void FUNCTION(test,trap) (void);

void
FUNCTION(test,trap) (void)
{
  int status ;

  BASE real_x ;
  BASE complex_x ;

  BASE * real_data = &real_x ;
  BASE * complex_data = &complex_x  ; 

  TYPE(gsl_fft_wavetable_complex) * cw;
  TYPE(gsl_fft_wavetable_real) * rw;
  TYPE(gsl_fft_wavetable_halfcomplex) * hcw;

  /* n = 0 in alloc */

  cw = FUNCTION(gsl_fft_complex,alloc) (0);
  gsl_test (cw != 0, "trap for n = 0 in " NAME(gsl_fft_complex) "_alloc");

  rw = FUNCTION(gsl_fft_real,alloc) (0);
  gsl_test (rw != 0, "trap for n = 0 in " NAME(gsl_fft_real) "_alloc" );

  hcw = FUNCTION(gsl_fft_halfcomplex,alloc) (0);
  gsl_test (hcw != 0, "trap for n = 0 in " NAME(gsl_fft_halfcomplex) "_alloc");

  cw = FUNCTION(gsl_fft_complex,alloc) (10);
  hcw = FUNCTION(gsl_fft_halfcomplex,alloc) (10);
  rw = FUNCTION(gsl_fft_real,alloc) (10);

  /* n = 0 in fft forward */

  status = FUNCTION(gsl_fft_complex,forward) (complex_data, 1, 0, cw);
  gsl_test (!status, "trap for n = 0 in " NAME(gsl_fft_complex) "_forward");

  status = FUNCTION(gsl_fft_real,transform) (real_data, 1, 0, rw);
  gsl_test (!status, "trap for n = 0 in " NAME(gsl_fft_real) "_transform");

  status = FUNCTION(gsl_fft_halfcomplex,transform) (real_data, 1, 0, hcw);
  gsl_test (!status, "trap for n = 0 in " NAME(gsl_fft_halfcomplex) "_transform");

  status = FUNCTION(gsl_fft_complex,radix2_forward) (complex_data, 1, 0);
  gsl_test (!status, "trap for n = 0 in " NAME(gsl_fft_complex) "_radix2_forward");

  /* n = 0 in fft backward */

  status = FUNCTION(gsl_fft_complex,backward) (complex_data, 1, 0, cw);
  gsl_test (!status, "trap for n = 0 in " NAME(gsl_fft_complex) "_backward");

  status = FUNCTION(gsl_fft_complex,radix2_backward) (complex_data, 1, 0);
  gsl_test (!status, "trap for n = 0 in " NAME(gsl_fft_complex) "_radix2_backward");

  /* n = 0 in fft inverse */

  status = FUNCTION(gsl_fft_complex,inverse) (complex_data, 1, 0, cw);
  gsl_test (!status, "trap for n = 0 in " NAME(gsl_fft_complex) "_inverse");

  status = FUNCTION(gsl_fft_complex,radix2_inverse) (complex_data, 1, 0);
  gsl_test (!status, "trap for n = 0 in " NAME(gsl_fft_complex) "_radix2_inverse");

  /* n != 2^k in power of 2 routines */

  status = FUNCTION(gsl_fft_complex,radix2_forward) (complex_data, 1, 17);
  gsl_test (!status, "trap for n != 2^k in " NAME(gsl_fft_complex) "_radix2_forward");

  status = FUNCTION(gsl_fft_complex,radix2_backward) (complex_data, 1, 17);
  gsl_test (!status, "trap for n != 2^k in " NAME(gsl_fft_complex) "_radix2_backward");

  status = FUNCTION(gsl_fft_complex,radix2_inverse) (complex_data, 1, 17);
  gsl_test (!status, "trap for n != 2^k in " NAME(gsl_fft_complex) "_radix2_inverse");

  /* n != wavetable.n in mixed radix routines */

  cw->n = 3;
  status = FUNCTION(gsl_fft_complex,forward) (complex_data, 1, 4, cw);
  gsl_test (!status, "trap for n != nw in " NAME(gsl_fft_complex) "_forward");

  cw->n = 3;
  status = FUNCTION(gsl_fft_complex,backward) (complex_data, 1, 4, cw);
  gsl_test (!status, "trap for n != nw in " NAME(gsl_fft_complex) "_backward");

  cw->n = 3;
  status = FUNCTION(gsl_fft_complex,inverse) (complex_data, 1, 4, cw);
  gsl_test (!status, "trap for n != nw in " NAME(gsl_fft_complex) "_inverse");

  rw->n = 3;
  status = FUNCTION(gsl_fft_real,transform) (real_data, 1, 4, rw);
  gsl_test (!status, "trap for n != nw in " NAME(gsl_fft_real) "_transform");

  hcw->n = 3;
  status = FUNCTION(gsl_fft_halfcomplex,transform) (real_data, 1, 4, hcw);
  gsl_test (!status, "trap for n != nw in " NAME(gsl_fft_halfcomplex) "_transform");

  FUNCTION (gsl_fft_halfcomplex,free) (hcw) ;
  FUNCTION (gsl_fft_real,free) (rw) ;
  FUNCTION (gsl_fft_complex,free) (cw) ;

}


