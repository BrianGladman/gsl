#include <math.h>
#include <gsl_complex.h>
#include <gsl_fft_complex.h>

int
gsl_fft_complex_pass_7 (const complex from[],
			complex to[],
			const gsl_fft_direction sign,
			const unsigned int product,
			const unsigned int n,
			const complex twiddle1[],
			const complex twiddle2[],
			const complex twiddle3[],
			const complex twiddle4[],
			const complex twiddle5[],
			const complex twiddle6[])
{

  unsigned int i = 0, j = 0;
  unsigned int k, k1;

  const unsigned int factor = 7;
  const unsigned int m = n / factor;
  const unsigned int q = n / product;
  const unsigned int product_1 = product / factor;
  const unsigned int jump = (factor - 1) * product_1;

  const double c1 = cos(1.0 * 2.0 * M_PI / 7.0) ;
  const double c2 = cos(2.0 * 2.0 * M_PI / 7.0) ;
  const double c3 = cos(3.0 * 2.0 * M_PI / 7.0) ;

  const double s1 = sin(1.0 * 2.0 * M_PI / 7.0) ;
  const double s2 = sin(2.0 * 2.0 * M_PI / 7.0) ;
  const double s3 = sin(3.0 * 2.0 * M_PI / 7.0) ;

  for (k = 0; k < q; k++)
    {
      double w1_real, w1_imag, w2_real, w2_imag, w3_real, w3_imag, w4_real,
        w4_imag, w5_real, w5_imag, w6_real, w6_imag;

      if (k == 0)
	{
	  w1_real = 1.0;
	  w1_imag = 0.0;
	  w2_real = 1.0;
	  w2_imag = 0.0;
	  w3_real = 1.0;
	  w3_imag = 0.0;
	  w4_real = 1.0;
	  w4_imag = 0.0;
	  w5_real = 1.0;
	  w5_imag = 0.0;
	  w6_real = 1.0;
	  w6_imag = 0.0;
	}
      else
	{
	  if (sign == forward)
	    {
	      /* forward tranform */
	      w1_real = twiddle1[k - 1].real;
	      w1_imag = twiddle1[k - 1].imag;
	      w2_real = twiddle2[k - 1].real;
	      w2_imag = twiddle2[k - 1].imag;
	      w3_real = twiddle3[k - 1].real;
	      w3_imag = twiddle3[k - 1].imag;
	      w4_real = twiddle4[k - 1].real;
	      w4_imag = twiddle4[k - 1].imag;
	      w5_real = twiddle5[k - 1].real;
	      w5_imag = twiddle5[k - 1].imag;
	      w6_real = twiddle6[k - 1].real;
	      w6_imag = twiddle6[k - 1].imag;
	    }
	  else
	    {
	      /* backward tranform: w -> conjugate(w) */
	      w1_real = twiddle1[k - 1].real;
	      w1_imag = -twiddle1[k - 1].imag;
	      w2_real = twiddle2[k - 1].real;
	      w2_imag = -twiddle2[k - 1].imag;
	      w3_real = twiddle3[k - 1].real;
	      w3_imag = -twiddle3[k - 1].imag;
	      w4_real = twiddle4[k - 1].real;
	      w4_imag = -twiddle4[k - 1].imag;
	      w5_real = twiddle5[k - 1].real;
	      w5_imag = -twiddle5[k - 1].imag;
	      w6_real = twiddle6[k - 1].real;
	      w6_imag = -twiddle6[k - 1].imag;
	    }
	}

      for (k1 = 0; k1 < product_1; k1++)
	{
	  complex z0, z1, z2, z3, z4, z5, z6;

	  {
	    const unsigned int from0 = i;
	    const unsigned int from1 = from0 + m;
	    const unsigned int from2 = from1 + m;
	    const unsigned int from3 = from2 + m;
	    const unsigned int from4 = from3 + m;
	    const unsigned int from5 = from4 + m;
	    const unsigned int from6 = from5 + m;

	    z0 = from[from0];
	    z1 = from[from1];
	    z2 = from[from2];
	    z3 = from[from3];
	    z4 = from[from4];
	    z5 = from[from5];
	    z6 = from[from6];
	  }

	  /* compute x = W(7) z */

	  {
	    /* t0 = z1 + z6 */
	    const double t0_real = z1.real + z6.real ;
	    const double t0_imag = z1.imag + z6.imag ; 

	    /* t1 = z1 - z6 */
	    const double t1_real = z1.real - z6.real ;
	    const double t1_imag = z1.imag - z6.imag ; 

	    /* t2 = z2 + z5 */
	    const double t2_real = z2.real + z5.real ;
	    const double t2_imag = z2.imag + z5.imag ; 

	    /* t3 = z2 - z5 */
	    const double t3_real = z2.real - z5.real ;
	    const double t3_imag = z2.imag - z5.imag ; 

	    /* t4 = z4 + z3 */
	    const double t4_real = z4.real + z3.real ;
	    const double t4_imag = z4.imag + z3.imag ; 

	    /* t5 = z4 - z3 */
	    const double t5_real = z4.real - z3.real ;
	    const double t5_imag = z4.imag - z3.imag ; 

	    /* t6 = t2 + t0 */
	    const double t6_real = t2_real + t0_real ;
	    const double t6_imag = t2_imag + t0_imag ;

	    /* t7 = t5 + t3 */
	    const double t7_real = t5_real + t3_real ;
	    const double t7_imag = t5_imag + t3_imag ;

	    /* b0 = z0 + t6 + t4 */
	    const double b0_real = z0.real + t6_real + t4_real ;
	    const double b0_imag = z0.imag + t6_imag + t4_imag ;

	    /* b1 = ((cos(2pi/7) + cos(4pi/7) + cos(6pi/7))/3-1) (t6 + t4) */
	    const double b1_real = (((c1 + c2 + c3)/3.0 - 1.0) * (t6_real + t4_real));
	    const double b1_imag = (((c1 + c2 + c3)/3.0 - 1.0) * (t6_imag + t4_imag));

	    /* b2 = ((2*cos(2pi/7) - cos(4pi/7) - cos(6pi/7))/3) (t0 - t4) */
	    const double b2_real = (((2.0 * c1 - c2 - c3)/3.0) * (t0_real - t4_real));
	    const double b2_imag = (((2.0 * c1 - c2 - c3)/3.0) * (t0_imag - t4_imag));

	    /* b3 = ((cos(2pi/7) - 2*cos(4pi/7) + cos(6pi/7))/3) (t4 - t2) */
	    const double b3_real = (((c1 - 2.0*c2 + c3)/3.0) * (t4_real - t2_real));
	    const double b3_imag = (((c1 - 2.0*c2 + c3)/3.0) * (t4_imag - t2_imag));
	    
	    /* b4 = ((cos(2pi/7) + cos(4pi/7) - 2*cos(6pi/7))/3) (t2 - t0) */
	    const double b4_real = (((c1 + c2 - 2.0 * c3)/3.0) * (t2_real - t0_real));
	    const double b4_imag = (((c1 + c2 - 2.0 * c3)/3.0) * (t2_imag - t0_imag));
	    
	    /* b5 = sign * ((sin(2pi/7) + sin(4pi/7) - sin(6pi/7))/3) (t7 + t1) */
	    const double b5_real = (-(int)sign) * ((s1 + s2 - s3)/3.0) * (t7_real + t1_real) ;
	    const double b5_imag = (-(int)sign) * ((s1 + s2 - s3)/3.0) * (t7_imag + t1_imag) ;

	    /* b6 = sign * ((2sin(2pi/7) - sin(4pi/7) + sin(6pi/7))/3) (t1 - t5) */
	    const double b6_real = (-(int)sign) * ((2.0 * s1 - s2 + s3)/3.0) * (t1_real - t5_real) ;
	    const double b6_imag = (-(int)sign) * ((2.0 * s1 - s2 + s3)/3.0) * (t1_imag - t5_imag) ;

	    /* b7 = sign * ((sin(2pi/7) - 2sin(4pi/7) - sin(6pi/7))/3) (t5 - t3) */
	    const double b7_real = (-(int)sign) * ((s1 - 2.0 * s2 - s3)/3.0) * (t5_real - t3_real) ;
	    const double b7_imag = (-(int)sign) * ((s1 - 2.0 * s2 - s3)/3.0) * (t5_imag - t3_imag) ;

	    /* b8 = sign * ((sin(2pi/7) + sin(4pi/7) + 2sin(6pi/7))/3) (t3 - t1) */
	    const double b8_real = (-(int)sign) * ((s1 + s2 + 2.0 * s3)/3.0) * (t3_real - t1_real) ;
	    const double b8_imag = (-(int)sign) * ((s1 + s2 + 2.0 * s3)/3.0) * (t3_imag - t1_imag) ;


	    /* T0 = b0 + b1 */
	    const double T0_real = b0_real + b1_real ;
	    const double T0_imag = b0_imag + b1_imag ;

	    /* T1 = b2 + b3 */
	    const double T1_real = b2_real + b3_real ;
	    const double T1_imag = b2_imag + b3_imag ;

	    /* T2 = b4 - b3 */
	    const double T2_real = b4_real - b3_real ;
	    const double T2_imag = b4_imag - b3_imag ;

	    /* T3 = -b2 - b4 */
	    const double T3_real = -b2_real - b4_real ;
	    const double T3_imag = -b2_imag - b4_imag ;

	    /* T4 = b6 + b7 */
	    const double T4_real = b6_real + b7_real ;
	    const double T4_imag = b6_imag + b7_imag ;

	    /* T5 = b8 - b7 */
	    const double T5_real = b8_real - b7_real ;
	    const double T5_imag = b8_imag - b7_imag ;

	    /* T6 = -b8 - b6 */
	    const double T6_real = -b8_real - b6_real ;
	    const double T6_imag = -b8_imag - b6_imag ;

	    /* T7 = T0 + T1 */
	    const double T7_real = T0_real + T1_real ;
	    const double T7_imag = T0_imag + T1_imag ;

	    /* T8 = T0 + T2 */
	    const double T8_real = T0_real + T2_real ;
	    const double T8_imag = T0_imag + T2_imag ;

	    /* T9 = T0 + T3 */
	    const double T9_real = T0_real + T3_real ;
	    const double T9_imag = T0_imag + T3_imag ;

	    /* T10 = T4 + b5 */
	    const double T10_real = T4_real + b5_real ;
	    const double T10_imag = T4_imag + b5_imag ;

	    /* T11 = T5 + b5 */
	    const double T11_real = T5_real + b5_real ;
	    const double T11_imag = T5_imag + b5_imag ;

	    /* T12 = T6 + b5 */
	    const double T12_real = T6_real + b5_real ;
	    const double T12_imag = T6_imag + b5_imag ;


	    /* x0 = b0 */
	    const double x0_real = b0_real ;
	    const double x0_imag = b0_imag ;
	    
	    /* x1 = T7 - i T10 */
	    const double x1_real = T7_real + T10_imag ;
	    const double x1_imag = T7_imag - T10_real ;
	    
	    /* x2 = T9 - i T12 */
	    const double x2_real = T9_real + T12_imag ;
	    const double x2_imag = T9_imag - T12_real ;

	    /* x3 = T8 + i T11 */
	    const double x3_real = T8_real - T11_imag ;
	    const double x3_imag = T8_imag + T11_real ;
	    
	    /* x4 = T8 - i T11 */
	    const double x4_real = T8_real + T11_imag ;
	    const double x4_imag = T8_imag - T11_real ;
	   
	    /* x5 = T9 + i T12 */
	    const double x5_real = T9_real - T12_imag ;
	    const double x5_imag = T9_imag + T12_real ;

	    /* x6 = T7 + i T10 */
	    const double x6_real = T7_real - T10_imag ;
	    const double x6_imag = T7_imag + T10_real ;

	  /* apply twiddle factors */

	    const unsigned int to0 = j;
	    const unsigned int to1 = to0 + product_1;
	    const unsigned int to2 = to1 + product_1;
	    const unsigned int to3 = to2 + product_1;
	    const unsigned int to4 = to3 + product_1;
	    const unsigned int to5 = to4 + product_1;
	    const unsigned int to6 = to5 + product_1;

	    /* to0 = 1 * x0 */
	    to[to0].real = x0_real;
	    to[to0].imag = x0_imag;

	    /* to1 = w1 * x1 */
	    to[to1].real = w1_real * x1_real - w1_imag * x1_imag;
	    to[to1].imag = w1_real * x1_imag + w1_imag * x1_real;

	    /* to2 = w2 * x2 */
	    to[to2].real = w2_real * x2_real - w2_imag * x2_imag;
	    to[to2].imag = w2_real * x2_imag + w2_imag * x2_real;

	    /* to3 = w3 * x3 */
	    to[to3].real = w3_real * x3_real - w3_imag * x3_imag;
	    to[to3].imag = w3_real * x3_imag + w3_imag * x3_real;

	    /* to4 = w4 * x4 */
	    to[to4].real = w4_real * x4_real - w4_imag * x4_imag;
	    to[to4].imag = w4_real * x4_imag + w4_imag * x4_real;

	    /* to5 = w5 * x5 */
	    to[to5].real = w5_real * x5_real - w5_imag * x5_imag;
	    to[to5].imag = w5_real * x5_imag + w5_imag * x5_real;

	    /* to6 = w6 * x6 */
	    to[to6].real = w6_real * x6_real - w6_imag * x6_imag;
	    to[to6].imag = w6_real * x6_imag + w6_imag * x6_real;

	  }
	  i++;
	  j++;
	}
      j += jump;
    }
  return 0;
}
