int
FUNCTION(fft_complex,pass_4) (const BASE in[],
			      const size_t istride,
			      BASE out[],
			      const size_t ostride,
			      const gsl_fft_direction sign,
			      const size_t product,
			      const size_t n,
			      const gsl_complex twiddle1[],
			      const gsl_complex twiddle2[],
			      const gsl_complex twiddle3[])
{
  size_t i = 0, j = 0;
  size_t k, k1;

  const size_t factor = 4;
  const size_t m = n / factor;
  const size_t q = n / product;
  const size_t p_1 = product / factor;
  const size_t jump = (factor - 1) * p_1;

  for (k = 0; k < q; k++)
    {
      double w1_real, w1_imag, w2_real, w2_imag, w3_real, w3_imag;

      if (k == 0)
	{
	  w1_real = 1.0;
	  w1_imag = 0.0;
	  w2_real = 1.0;
	  w2_imag = 0.0;
	  w3_real = 1.0;
	  w3_imag = 0.0;
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

	  /* compute x = W(4) z */
	  
	  /* t1 = z0 + z2 */
	  const double t1_real = z0_real + z2_real;
	  const double t1_imag = z0_imag + z2_imag;
	  
	  /* t2 = z1 + z3 */
	  const double t2_real = z1_real + z3_real;
	  const double t2_imag = z1_imag + z3_imag;
	  
	  /* t3 = z0 - z2 */
	  const double t3_real = z0_real - z2_real;
	  const double t3_imag = z0_imag - z2_imag;
	  
	  /* t4 = (+/-) (z1 - z3) */
	  const double t4_real = ((int) sign) * (z1_real - z3_real);
	  const double t4_imag = ((int) sign) * (z1_imag - z3_imag);

	    /* x0 = t1 + t2 */
	  const double x0_real = t1_real + t2_real;
	  const double x0_imag = t1_imag + t2_imag;

	    /* x1 = t3 + i t4 */
	  const double x1_real = t3_real - t4_imag;
	  const double x1_imag = t3_imag + t4_real;

	    /* x2 = t1 - t2 */
	  const double x2_real = t1_real - t2_real;
	  const double x2_imag = t1_imag - t2_imag;

	    /* x3 = t3 - i t4 */
	  const double x3_real = t3_real + t4_imag;
	  const double x3_imag = t3_imag - t4_real;

	  /* apply twiddle factors */

	  /* to0 = 1 * x0 */
	  REAL(out,ostride,j) = x0_real;
	  IMAG(out,ostride,j) = x0_imag;

	  /* to1 = w1 * x1 */
	  REAL(out, ostride, j + p_1) = w1_real * x1_real - w1_imag * x1_imag;
	  IMAG(out, ostride, j + p_1) = w1_real * x1_imag + w1_imag * x1_real;
	  
	  /* to2 = w2 * x2 */
	  REAL(out, ostride, j + 2 * p_1) = w2_real * x2_real - w2_imag * x2_imag;
	  IMAG(out, ostride, j + 2 * p_1) = w2_real * x2_imag + w2_imag * x2_real;
	  
	  /* to3 = w3 * x3 */
	  REAL(out, ostride, j + 3 * p_1) = w3_real * x3_real - w3_imag * x3_imag;
	  IMAG(out, ostride, j + 3 * p_1) = w3_real * x3_imag + w3_imag * x3_real;
	  
	  i++;
	  j++;
	}
      j += jump;
    }
  return 0;
}
