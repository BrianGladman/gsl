static int
FUNCTION(fft_complex,pass_n) (BASE in[],
			      const size_t istride,
			      BASE out[],
			      const size_t ostride,
			      const gsl_fft_direction sign,
			      const size_t factor,
			      const size_t product,
			      const size_t n,
			      const gsl_complex twiddle[])
{
  size_t i = 0, j = 0;
  size_t k, k1;

  const size_t m = n / factor;
  const size_t q = n / product;
  const size_t p_1 = product / factor;
  const size_t jump = (factor - 1) * p_1;

  size_t e, e1;

  for (i = 0; i < m; i++)
    {
      REAL(out,ostride,i) = REAL(in,istride,i);
      IMAG(out,ostride,i) = IMAG(in,istride,i);
    }

  for (e = 1; e < (factor - 1) / 2 + 1; e++)
    {
      for (i = 0; i < m; i++)
	{
	  const size_t idx = i + e * m;
	  const size_t idxc = i + (factor - e) * m;
	  REAL(out,ostride,idx) = REAL(in,istride,idx) + REAL(in,istride,idxc);
	  IMAG(out,ostride,idx) = IMAG(in,istride,idx) + IMAG(in,istride,idxc);
	  REAL(out,ostride,idxc) = REAL(in,istride,idx) - REAL(in,istride,idxc);
	  IMAG(out,ostride,idxc) = IMAG(in,istride,idx) - IMAG(in,istride,idxc);
	}
    }

  /* e = 0 */

  for (i=0 ; i<m; i++) 
    {
      REAL(in,istride,i) = REAL(out,ostride,i);
      IMAG(in,istride,i) = IMAG(out,ostride,i);
    }

  for (e1 = 1; e1 < (factor - 1) / 2 + 1; e1++)
    {
      for (i = 0; i < m; i++)
	{
	  REAL(in,istride,i) += REAL(out,ostride,i + e1*m) ;
	  IMAG(in,istride,i) += IMAG(out,ostride,i + e1*m) ;
	}
    }

  for (e = 1; e < (factor-1)/2 + 1; e++)
    {
      size_t idx = e*q ;
      const size_t idx_step = e * q ;
      double w_real, w_imag ;

      const size_t em = e * m ;
      const size_t ecm = (factor - e) * m ;

      for (i = 0; i < m; i++) 
	{
	  REAL(in,istride,i+em) = REAL(out,ostride,i) ;
	  IMAG(in,istride,i+em) = IMAG(out,ostride,i) ;
	  REAL(in,istride,i+ecm) = REAL(out,ostride,i) ;
	  IMAG(in,istride,i+ecm) = IMAG(out,ostride,i) ;
	}

      for (e1 = 1; e1 < (factor - 1) / 2 + 1; e1++)
	{
	  if (idx == 0) {
	    w_real = 1 ;
	    w_imag = 0 ;
	  } else {
	    if (sign == forward) {
	      w_real = GSL_REAL(twiddle[idx - 1]) ;
	      w_imag = GSL_IMAG(twiddle[idx - 1]) ;
	    } else {
	      w_real = GSL_REAL(twiddle[idx - 1]) ;
	      w_imag = -GSL_IMAG(twiddle[idx - 1]) ;
	    }
	  }

	  for (i = 0; i < m; i++) 
	    {
	      const double xp_real = REAL(out,ostride,i + e1 * m);
	      const double xp_imag = IMAG(out,ostride,i + e1 * m);
	      const double xm_real = REAL(out,ostride,i + (factor - e1) *m);
	      const double xm_imag = IMAG(out,ostride,i + (factor - e1) *m);
	
	      const double ap = w_real * xp_real ;
	      const double am = w_imag * xm_imag ; 

	      double sum_real = ap - am;
	      double sumc_real = ap + am;

	      const double bp = w_real * xp_imag ;
	      const double bm = w_imag * xm_real ;

	      double sum_imag = bp + bm;
	      double sumc_imag = bp - bm;

	      REAL(in,istride,i + em) += sum_real;
	      IMAG(in,istride,i + em) += sum_imag;
	      REAL(in,istride,i + ecm) += sumc_real;
	      IMAG(in,istride,i + ecm) += sumc_imag;
	    }
	  idx += idx_step ;
	  idx %= factor * q ;
	}
    }

  i = 0;
  j = 0;

  /* k = 0 */
  for (k1 = 0; k1 < p_1; k1++)
    {
      REAL(out,ostride,k1) = REAL(in,istride,k1);
      IMAG(out,ostride,k1) = IMAG(in,istride,k1);
    }

  for (e1 = 1; e1 < factor; e1++)
    {
      for (k1 = 0; k1 < p_1; k1++)
	{
	  REAL(out,ostride,k1 + e1 * p_1) = REAL(in,istride,k1 + e1 * m) ;
	  IMAG(out,ostride,k1 + e1 * p_1) = IMAG(in,istride,k1 + e1 * m) ;
	}
    }

  i = p_1 ;
  j = product ;

  for (k = 1; k < q; k++)
    {
      for (k1 = 0; k1 < p_1; k1++)
	{
	  REAL(out,ostride,j) = REAL(in,istride,i);
	  IMAG(out,ostride,j) = IMAG(in,istride,i);
	  i++;
	  j++;
	}
      j += jump;
    }

  i = p_1 ;
  j = product ;

  for (k = 1; k < q; k++)
    {
      for (k1 = 0; k1 < p_1; k1++)
	{
	  for (e1 = 1; e1 < factor; e1++)
	    {
	      double x_real = REAL(in, istride,i + e1 * m);
	      double x_imag = IMAG(in, istride,i + e1 * m);

	      double w_real, w_imag ;
	      if (sign == forward) {
		w_real = GSL_REAL(twiddle[(e1-1)*q + k-1]) ;
		w_imag = GSL_IMAG(twiddle[(e1-1)*q + k-1]) ;
	      } else {
		w_real = GSL_REAL(twiddle[(e1-1)*q + k-1]) ;
		w_imag = -GSL_IMAG(twiddle[(e1-1)*q + k-1]) ; 
	      }

	      REAL(out,ostride,j + e1 * p_1) = w_real * x_real - w_imag * x_imag;
	      IMAG(out,ostride,j + e1 * p_1) = w_real * x_imag + w_imag * x_real;
	    }
	  i++;
	  j++;
	}
      j += jump;
    }

  return 0;
}

