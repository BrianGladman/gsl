#include <config.h>
#include <math.h>

#include <gsl_complex.h>
#include <gsl_fft_complex.h>

#include "fft_complex.h"

int
fft_complex_pass_7 (const double in[],
		    const size_t istride,
		    double out[],
		    const size_t ostride,
		    const gsl_fft_direction sign,
		    const size_t product,
		    const size_t n,
		    const gsl_complex twiddle1[],
		    const gsl_complex twiddle2[],
		    const gsl_complex twiddle3[],
		    const gsl_complex twiddle4[],
		    const gsl_complex twiddle5[],
		    const gsl_complex twiddle6[])
{

  size_t i = 0, j = 0;
  size_t k, k1;

  const size_t factor = 7;
  const size_t m = n / factor;
  const size_t q = n / product;
  const size_t p_1 = product / factor;
  const size_t jump = (factor - 1) * p_1;

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
	      w1_real = GSL_REAL(twiddle1[k - 1]);
	      w1_imag = GSL_IMAG(twiddle1[k - 1]);
	      w2_real = GSL_REAL(twiddle2[k - 1]);
	      w2_imag = GSL_IMAG(twiddle2[k - 1]);
	      w3_real = GSL_REAL(twiddle3[k - 1]);
	      w3_imag = GSL_IMAG(twiddle3[k - 1]);
	      w4_real = GSL_REAL(twiddle4[k - 1]);
	      w4_imag = GSL_IMAG(twiddle4[k - 1]);
	      w5_real = GSL_REAL(twiddle5[k - 1]);
	      w5_imag = GSL_IMAG(twiddle5[k - 1]);
	      w6_real = GSL_REAL(twiddle6[k - 1]);
	      w6_imag = GSL_IMAG(twiddle6[k - 1]);
	    }
	  else
	    {
	      /* backward tranform: w -> conjugate(w) */
	      w1_real = GSL_REAL(twiddle1[k - 1]);
	      w1_imag = -GSL_IMAG(twiddle1[k - 1]);
	      w2_real = GSL_REAL(twiddle2[k - 1]);
	      w2_imag = -GSL_IMAG(twiddle2[k - 1]);
	      w3_real = GSL_REAL(twiddle3[k - 1]);
	      w3_imag = -GSL_IMAG(twiddle3[k - 1]);
	      w4_real = GSL_REAL(twiddle4[k - 1]);
	      w4_imag = -GSL_IMAG(twiddle4[k - 1]);
	      w5_real = GSL_REAL(twiddle5[k - 1]);
	      w5_imag = -GSL_IMAG(twiddle5[k - 1]);
	      w6_real = GSL_REAL(twiddle6[k - 1]);
	      w6_imag = -GSL_IMAG(twiddle6[k - 1]);
	    }
	}

      for (k1 = 0; k1 < p_1; k1++)
	{
	  const double z0_real = REAL(in,istride,i);
	  const double z0_imag = IMAG(in,istride,i);
	  const double z1_real = REAL(in,istride,i+m);
	  const double z1_imag = IMAG(in,istride,i+m);
	  const double z2_real = REAL(in,istride,i+2*m);
	  const double z2_imag = IMAG(in,istride,i+2*m);
	  const double z3_real = REAL(in,istride,i+3*m);
	  const double z3_imag = IMAG(in,istride,i+3*m);
	  const double z4_real = REAL(in,istride,i+4*m);
	  const double z4_imag = IMAG(in,istride,i+4*m);
	  const double z5_real = REAL(in,istride,i+5*m);
	  const double z5_imag = IMAG(in,istride,i+5*m);
	  const double z6_real = REAL(in,istride,i+6*m);
	  const double z6_imag = IMAG(in,istride,i+6*m);

	  /* compute x = W(7) z */
	  
	  /* t0 = z1 + z6 */
	  const double t0_real = z1_real + z6_real ;
	  const double t0_imag = z1_imag + z6_imag ; 
	  
	  /* t1 = z1 - z6 */
	  const double t1_real = z1_real - z6_real ;
	  const double t1_imag = z1_imag - z6_imag ; 
	  
	  /* t2 = z2 + z5 */
	  const double t2_real = z2_real + z5_real ;
	  const double t2_imag = z2_imag + z5_imag ; 
	  
	  /* t3 = z2 - z5 */
	  const double t3_real = z2_real - z5_real ;
	  const double t3_imag = z2_imag - z5_imag ; 
	  
	  /* t4 = z4 + z3 */
	  const double t4_real = z4_real + z3_real ;
	  const double t4_imag = z4_imag + z3_imag ; 
	  
	  /* t5 = z4 - z3 */
	  const double t5_real = z4_real - z3_real ;
	  const double t5_imag = z4_imag - z3_imag ; 
	  
	  /* t6 = t2 + t0 */
	  const double t6_real = t2_real + t0_real ;
	  const double t6_imag = t2_imag + t0_imag ;
	  
	  /* t7 = t5 + t3 */
	  const double t7_real = t5_real + t3_real ;
	  const double t7_imag = t5_imag + t3_imag ;
	  
	  /* b0 = z0 + t6 + t4 */
	  const double b0_real = z0_real + t6_real + t4_real ;
	  const double b0_imag = z0_imag + t6_imag + t4_imag ;
	  
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
	  
	  /* to0 = 1 * x0 */
	  REAL(out,ostride,j) = x0_real;
	  IMAG(out,ostride,j) = x0_imag;
	  
	  /* to1 = w1 * x1 */
	  REAL(out,ostride,j+p_1) = w1_real * x1_real - w1_imag * x1_imag;
	  IMAG(out,ostride,j+p_1) = w1_real * x1_imag + w1_imag * x1_real;

	  /* to2 = w2 * x2 */
	  REAL(out,ostride,j+2*p_1) = w2_real * x2_real - w2_imag * x2_imag;
	  IMAG(out,ostride,j+2*p_1) = w2_real * x2_imag + w2_imag * x2_real;

	  /* to3 = w3 * x3 */
	  REAL(out,ostride,j+3*p_1) = w3_real * x3_real - w3_imag * x3_imag;
	  IMAG(out,ostride,j+3*p_1) = w3_real * x3_imag + w3_imag * x3_real;

	  /* to4 = w4 * x4 */
	  REAL(out,ostride,j+4*p_1) = w4_real * x4_real - w4_imag * x4_imag;
	  IMAG(out,ostride,j+4*p_1) = w4_real * x4_imag + w4_imag * x4_real;

	  /* to5 = w5 * x5 */
	  REAL(out,ostride,j+5*p_1) = w5_real * x5_real - w5_imag * x5_imag;
	  IMAG(out,ostride,j+5*p_1) = w5_real * x5_imag + w5_imag * x5_real;
	  
	  /* to6 = w6 * x6 */
	  REAL(out,ostride,j+6*p_1) = w6_real * x6_real - w6_imag * x6_imag;
	  IMAG(out,ostride,j+6*p_1) = w6_real * x6_imag + w6_imag * x6_real;
	  
	  i++; j++;
	}
      j += jump;
    }
  return 0;
}
