/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */

  size_t i, j;
  size_t k = 0;
  const BASE_TYPE aR = REAL0(alpha);
  const BASE_TYPE aI = IMAG0(alpha);
  const BASE_TYPE conj = -1.0;

  if(Uplo == CblasUpper) {
    for(i=0; i<N; i++) {
      const BASE_TYPE XiR = REAL(X, incX, i);
      const BASE_TYPE XiI = IMAG(X, incX, i);
      const BASE_TYPE YiR = REAL(Y, incY, i);
      const BASE_TYPE YiI = IMAG(Y, incY, i);
      for(j=i; j<N; j++) {
	const BASE_TYPE XjR = REAL(X, incX, j);
	const BASE_TYPE XjI = IMAG(X, incX, j);
	const BASE_TYPE YjR = REAL(Y, incY, j);
	const BASE_TYPE YjI = IMAG(Y, incY, j);
        const BASE_TYPE tmpijR = XiR * YjR - conj * XiI * YjI;
	const BASE_TYPE tmpijI = XiI * YjR + conj * XiR * YjI;
	const BASE_TYPE tmpjiR = XjR * YiR - conj * XjI * YiI;
	const BASE_TYPE tmpjiI = XjR * YiI + conj * XjI * YiR;
        REAL(Ap, 1, k) += aR*tmpijR - aI*tmpijI + aR*tmpjiR - conj*aI*tmpjiI;
	IMAG(Ap, 1, k) += aR*tmpijI + aI*tmpijR + aR*tmpjiI + conj*aI*tmpjiR;
	k++;
      }
    }
  }
  else {
    for(i=0; i<N; i++) {
      const BASE_TYPE XiR = REAL(X, incX, i);
      const BASE_TYPE XiI = IMAG(X, incX, i);
      const BASE_TYPE YiR = REAL(Y, incY, i);
      const BASE_TYPE YiI = IMAG(Y, incY, i);
      for(j=0; j<=i; j++) {
	const BASE_TYPE XjR = REAL(X, incX, j);
	const BASE_TYPE XjI = IMAG(X, incX, j);
	const BASE_TYPE YjR = REAL(Y, incY, j);
	const BASE_TYPE YjI = IMAG(Y, incY, j);
        const BASE_TYPE tmpijR = XiR * YjR - conj * XiI * YjI;
	const BASE_TYPE tmpijI = XiI * YjR + conj * XiR * YjI;
	const BASE_TYPE tmpjiR = XjR * YiR - conj * XjI * YiI;
	const BASE_TYPE tmpjiI = XjR * YiI + conj * XjI * YiR;
        REAL(Ap, 1, k) += aR*tmpijR - aI*tmpijI + aR*tmpjiR - conj*aI*tmpjiI;
	IMAG(Ap, 1, k) += aR*tmpijI + aI*tmpijR + aR*tmpjiI + conj*aI*tmpjiR;
	k++;
      }
    }
  }
