##  After running this script with 
##
##     octave test.m
##
## you can trim the size of the file by stripping trailing zeros using
##
## perl -i.bak -p -e 's/(\d\.\d+?)(0+)(\D)/$1$3/g' test_*.c
##
##

rand("seed", 1);
global TEST=0;
global FILE=stdout;
global LIST;
global HEADER;

[LIST] = fopen("tests.c", "w");
[HEADER] = fopen("tests.h", "w");

######################################################################

function S = context (i)
  switch(i)
    case 1
      S.prefix = "s";
      S.precision = "float";
      S.complex = 0;
    case 2
      S.prefix = "d";
      S.precision = "double";
      S.complex = 0;
    case 3
      S.prefix = "c";
      S.precision = "float";
      S.complex = 1;
    case 4
      S.prefix = "z";
      S.precision = "double";
      S.complex = 1;
    otherwise
      error ("unrecognized case %d\n", i);
  endswitch
endfunction

function t = Trans (S)
  if (S.complex == 0)
    t = [111, 112];
    #t = v(1+fix(rand()));
  else
    t = [111, 112, 113];
    #t = v(1+fix(rand()));
  endif
endfunction

function c = coeff (S)
  if (S.complex == 0)
    v = [0, 0.1, 1, -1, -0.3];
    c = v(1+fix(rand()*length(v)));
  else
    v = [0 + 0*I, 0.1 * I, 1, I,  -1, -0.3 + 0.1 * I];
    c = v(1+fix(rand()*length(v)));
  endif
endfunction

function T = test_vectors (S, j)
  T.n = fix(j/4)+1;
  T.s1 = (-1)**(rem(fix(j/2),2));
  T.s2 = (-1)**(rem(j,2)) ;

  if (S.complex == 0)
    T.v1 = random_vector(T.n);
    T.v2 = random_vector(T.n);
  else
    T.v1 = random_vector(T.n) + I * random_vector(T.n);
    T.v2 = random_vector(T.n) + I * random_vector(T.n);
  endif
endfunction

function T = test_vector (S, j)
  T.n = fix(j/2)+1;
  T.s = (-1)**(rem(j,2)) ;

  if (S.complex == 0)
    T.v = random_vector(T.n);
  else
    T.v = random_vector(T.n) + I * random_vector(T.n);
  endif
endfunction

function T = test_matvectors (S, j)
  T.m = fix(j/4)+1;
  T.n = fix(j/4)+1;
  T.lda = fix(j/4)+1;
  T.s1 = (-1)**(rem(fix(j/2),2));
  T.s2 = (-1)**(rem(j,2)) ;

  if (S.complex == 0)
    T.A = random_matrix(T.m, T.n);
    T.v1 = random_vector(T.n);
    T.v2 = random_vector(T.m);
  else
    T.A = random_matrix(T.m, T.n) + I * random_matrix(T.m, T.n);
    T.v1 = random_vector(T.n) + I * random_vector(T.n);
    T.v2 = random_vector(T.m) + I * random_vector(T.m);
  endif
endfunction

function T = test_trmatvector (S, j)
  T.n = fix(j/4)+1;
  T.lda = fix(j/4)+1;
  T.s = (-1)**(rem(fix(j/2),2));

  if (S.complex == 0)
    T.A = random_matrix(T.n, T.n);
    T.v = random_vector(T.n);
  else
    T.A = random_matrix(T.n, T.n) + I * random_matrix(T.n, T.n);
    T.v = random_vector(T.n) + I * random_vector(T.n);
  endif
endfunction

function T = test_trmatvectors (S, j)
  T.n = fix(j/4)+1;
  T.lda = fix(j/4)+1;
  T.s1 = (-1)**(rem(fix(j/2),2));
  T.s2 = (-1)**(rem(j,2)) ;

  if (S.complex == 0)
    T.A = random_matrix(T.n, T.n);
    T.v1 = random_vector(T.n);
    T.v2 = random_vector(T.n);
  else
    T.A = random_matrix(T.n, T.n) + I * random_matrix(T.n, T.n);
    T.v1 = random_vector(T.n) + I * random_vector(T.n);
    T.v2 = random_vector(T.n) + I * random_vector(T.n);
  endif
endfunction


function T = test_bmatvectors (S, j, trans)
  T.kl = fix(j/4)+1;
  T.ku = fix(j/4)+1;
  b = T.kl+T.ku+1;

  T.m =  b;
  T.n =  1+j;
  T.lda = max([T.m,T.n]);
  T.s1 = (-1)**(rem(fix(j/2),2));
  T.s2 = (-1)**(rem(j,2)) ;

  if (S.complex == 0)
    T.A = random_matrix(T.lda, T.lda);
    T.v1 = random_vector(T.n);
    T.v2 = random_vector(T.m);
  else
    T.A = random_matrix(T.lda, T.lda)+ I * random_matrix(T.lda, T.lda);
    T.v1 = random_vector(T.n) + I * random_vector(T.n);
    T.v2 = random_vector(T.m) + I * random_vector(T.m);
  endif
  
  if (trans != 111)
    tmp = T.v1;
    T.v1 = T.v2;
    T.v2 = tmp;
  endif
endfunction


function T = test_tbmatvector (S, j)
  T.k = fix(j/4)+1;
  b = T.k+1;
  T.n =  b+fix(j/2);
  T.lda = T.n;
  T.s = (-1)**(fix(j/2));

  if (S.complex == 0)
    T.A = random_matrix(T.lda, T.lda);
    T.v = random_vector(T.n);
  else
    T.A = random_matrix(T.lda, T.lda)+ I * random_matrix(T.lda, T.lda);
    T.v = random_vector(T.n) + I * random_vector(T.n);
  endif
endfunction

function T = test_tpmatvector (S, j)
  T.n =  1+ fix(j/2);
  T.s = (-1)**(fix(j/2));

  N = T.n * (T.n + 1) / 2;

  if (S.complex == 0)
    T.A = random_vector(N);
    T.v = random_vector(T.n);
  else
    T.A = random_vector(N)+ I * random_vector(N);
    T.v = random_vector(T.n) + I * random_vector(T.n);
  endif
endfunction

function T = test_tpmatvectors (S, j)
  T.n =  1+ fix(j/2);
  T.s1 = (-1)**(rem(fix(j/2),2));
  T.s2 = (-1)**(rem(j,2)) ;

  N = T.n * (T.n + 1) / 2;

  if (S.complex == 0)
    T.A = random_vector(N);
    T.v1 = random_vector(T.n);
    T.v2 = random_vector(T.n);
  else
    T.A = random_vector(N) + I * random_vector(N);
    T.v1 = random_vector(T.n) + I * random_vector(T.n);
    T.v2 = random_vector(T.n) + I * random_vector(T.n);
  endif
endfunction



function T = test_symatvectors (S, j)
  T.n = fix(j/4)+1;
  T.lda = fix(j/4)+1;
  T.s1 = (-1)**(rem(fix(j/2),2));
  T.s2 = (-1)**(rem(j,2)) ;

  if (S.complex == 0)
    T.A = random_matrix(T.n, T.n);
    T.v1 = random_vector(T.n);
    T.v2 = random_vector(T.n);
  else
    T.A = random_matrix(T.n, T.n) + I * random_matrix(T.n, T.n);
    T.v1 = random_vector(T.n) + I * random_vector(T.n);
    T.v2 = random_vector(T.n) + I * random_vector(T.n);
  endif
endfunction

function T = test_sbmatvectors (S, j)
  T.k = fix(j/4)+1;
  b = T.k+1;
  T.n =  b+fix(j/2);
  T.lda = T.n;

  T.s1 = (-1)**(rem(fix(j/2),2));
  T.s2 = (-1)**(rem(j,2)) ;

  if (S.complex == 0)
    T.A = random_matrix(T.n, T.n);
    T.v1 = random_vector(T.n);
    T.v2 = random_vector(T.n);
  else
    T.A = random_matrix(T.n, T.n) + I * random_matrix(T.n, T.n);
    T.v1 = random_vector(T.n) + I * random_vector(T.n);
    T.v2 = random_vector(T.n) + I * random_vector(T.n);
  endif
endfunction

function T = test_spmatvectors (S, j)
  T.n =  1+ fix(j/2);
  N = T.n * (T.n + 1) / 2;

  T.s1 = (-1)**(rem(fix(j/2),2));
  T.s2 = (-1)**(rem(j,2)) ;

  if (S.complex == 0)
    T.A = random_vector(N);
    T.v1 = random_vector(T.n);
    T.v2 = random_vector(T.n);
  else
    T.A = random_vector(N) + I * random_vector(N);
    T.v1 = random_vector(T.n) + I * random_vector(T.n);
    T.v2 = random_vector(T.n) + I * random_vector(T.n);
  endif
endfunction


function T = test_matmat (S, j, order, transA, transB)
  T.m = fix(j/4)+1;
  T.n = fix(j/2)+1;
  T.k = fix(j/1)+1;

  if (order == 101)
    T.ldc = T.n;
  else
    T.ldc = T.m;
  endif

  if ((order == 101 && transA == 111) || (order == 102 && transA != 111))
    T.lda = T.k;
  else
    T.lda = T.m;
  endif

  if ((order == 101 && transB == 111) || (order == 102 && transB != 111))
    T.ldb = T.n;
  else
    T.ldb = T.k;
  endif

  if (S.complex == 0)
    T.C = random_matrix(T.m, T.n);

    T.A = random_matrix(T.m, T.k);
    T.B = random_matrix(T.k, T.n);
  else
    T.C = random_matrix(T.m, T.n) + I * random_matrix(T.m, T.n);

    T.A = random_matrix(T.m, T.k) + I * random_matrix(T.m, T.k);
    T.B = random_matrix(T.k, T.n) + I * random_matrix(T.k, T.n);
  endif
endfunction


function T = test_symatmat (S, j, order, side)
  T.m = fix(j/4)+1;
  T.n = fix(j/2)+1;

  if (side == 141)
    T.k = T.m;
  else
    T.k = T.n;
  endif

  T.lda = T.k;

  if (order == 101)
    T.ldc = T.n;
  else
    T.ldc = T.m;
  endif

  if (order == 101)
    T.ldb = T.n;
  else
    T.ldb = T.m;
  endif

  if (S.complex == 0)
    T.C = random_matrix(T.m, T.n);

    T.A = random_matrix(T.k, T.k);
    T.B = random_matrix(T.m, T.n);
  else
    T.C = random_matrix(T.m, T.n) + I * random_matrix(T.m, T.n);

    T.A = random_matrix(T.k, T.k) + I * random_matrix(T.k, T.k);
    T.B = random_matrix(T.m, T.n) + I * random_matrix(T.m, T.n);
  endif
endfunction


function T = test_syrkmatmat (S, j, order, trans)
  T.k = fix(j/4)+1;
  T.n = fix(j/2)+1;

  T.ldc = T.n;

  if ((order == 101 && trans == 111) || (order != 101 && trans != 111))
    T.lda = T.k;
  else
    T.lda = T.n;
  endif

  if (S.complex == 0)
    T.C = random_matrix(T.n, T.n);

    T.A = random_matrix(T.n, T.k);
  else
    T.C = random_matrix(T.n, T.n) + I * random_matrix(T.n, T.n);

    T.A = random_matrix(T.n, T.k) + I * random_matrix(T.n, T.k);
  endif
endfunction


function T = test_syr2kmatmat (S, j, order, trans)
  T.n = fix(j/4)+1;
  T.k = fix(j/2)+1;

  T.ldc = T.n;

  if (trans == 111)
    size1 = T.n;
    size2 = T.k;
  else
    size1 = T.k;
    size2 = T.n;
  endif

  if (order == 101)
    T.lda = size2;
  else
    T.lda = size1;
  endif

  if (order == 101)
    T.ldb = size2;
  else
    T.ldb = size1;
  endif

  if (S.complex == 0)
    T.C = random_matrix(T.n, T.n);

    T.A = random_matrix(T.n, T.k);
    T.B = random_matrix(T.n, T.k);
  else
    T.C = random_matrix(T.n, T.n) + I * random_matrix(T.n, T.n);

    T.A = random_matrix(T.n, T.k) + I * random_matrix(T.n, T.k);
    T.B = random_matrix(T.n, T.k) + I * random_matrix(T.n, T.k);
  endif
endfunction


function T = test_trmatmat (S, j, order, side, trans)
  T.m = fix(j/2)+1;
  T.n = j;

  if (order == 101)
    T.ldb = T.n;
  else
    T.ldb = T.m;
  endif

  if (side == 141)
    T.k = T.m;
  else
    T.k = T.n;
  endif

  T.lda = T.k;

  if (S.complex == 0)
    T.B = random_matrix(T.m, T.n);

    T.A = random_matrix(T.k, T.k);
  else
    T.B = random_matrix(T.m, T.n) + I * random_matrix(T.m, T.n);

    T.A = random_matrix(T.k, T.k) + I * random_matrix(T.k, T.k);
  endif
endfunction

function r = rnd ()
  r = 10**((rand()-0.5)*10);
  if (rand() < 0.5)
    r = -r;
  endif
endfunction

function v = random_vector(n)
  v = fix((rand(1,n)-0.5)*2000)/1000;
endfunction

function a = random_matrix(m,n)
  a = random_vector(m*n);
endfunction

function v = vector (X, incX, N)
  if (incX > 0)
    v = (X(1+(0:(N-1))*incX)).';
  elseif (incX < 0)
    v = (X(1+((N-1):-1:0)*(-incX))).';
  else
    v = [];
  endif
endfunction

function v = vout (X, incX, N, x)
  v = X;
  if (incX > 0)
    v((1:N)*incX) = x.';
  elseif (incX < 0)
    v((N:-1:1)*(-incX)) = x.';
  else
    v = [];
  endif
endfunction

function m = matrix (order, A, lda, M, N)
#    order
#    A
#   lda
#    M
#    N
  if (order == 102)   # column major
    tmp = reshape(A,lda,length(A)/lda);
    m = tmp(1:M,1:N);
  elseif (order == 101)  # row major
    tmp = reshape(A,lda,length(A)/lda).';
    m = tmp(1:M,1:N);
  endif
endfunction

function a = mout (order, A, lda, M, N, m)
  if (order == 102)  # column major
    a = reshape(A, lda, length(A)/lda);
    a(1:M, 1:N) = m;
    a = reshape(a, 1, length(A));
  elseif (order == 101)        # row major
    a = reshape(A,lda,length(A)/lda).';
    a(1:M, 1:N) = m;
    a = reshape(a.', 1, length(A));
  endif
endfunction

function m = trmatrix (order, uplo, Diag, A, lda, N)
  if (order == 102)   # column major
    tmp = reshape(A,lda,length(A)/lda);
    m = tmp(1:N,1:N);
  elseif (order == 101)  # row major
    tmp = reshape(A,lda,length(A)/lda).';
    m = tmp(1:N,1:N);
  endif

  if (uplo == 121) # upper
    m = triu(m);
  else  # lower
    m = tril(m);
  endif
  
  if (Diag == 132)  # unit diag
    m = m - diag(diag(m)) + eye(size(m));
  endif

endfunction

function m = trmatout (order, uplo, Diag, A, lda, N, a)
  if (order == 102)   # column major
    tmp = reshape(A,lda,length(A)/lda);
    m = tmp(1:N,1:N);
  elseif (order == 101)  # row major
    tmp = reshape(A,lda,length(A)/lda).';
    m = tmp(1:N,1:N);
  endif

  if (Diag == 132)  # unit diag
    if (uplo == 121) # upper
      for i = 1:(N-1)
        m(i,(i+1):N) = a(i,(i+1):N);
      endfor
    else  # lower
      for i = 2:N
        m(i,1:(i-1)) = a(i,1:(i-1));
      endfor
    endif
  else
    if (uplo == 121) # upper
      for i = 1:N
        m(i,i:N) = a(i,i:N);
      endfor
    else  # lower
      for i = 1:N
        m(i,1:i) = a(i,1:i);
      endfor
    endif
  endif

  if (order == 102)   # column major
    m = reshape(m,1,length(A));
  elseif (order == 101)  # row major
    m = reshape(m.',1,length(A));
  endif
endfunction

function m = tbmatrix (order, uplo, Diag, A, lda, N, K)

  if (uplo == 121) # upper
    m = bandmatrix(order, A, lda, N, N, 0, K);
    m = triu(m);
  else  # lower
    m = bandmatrix(order, A, lda, N, N, K, 0);
    m = tril(m);
  endif
  
  if (Diag == 132)  # unit diag
    m = m - diag(diag(m)) + eye(size(m));
  endif
endfunction

function m = tpmatrix (order, uplo, Diag, A, N)
  m = zeros(N,N);
  if (order == 102) # column major 
    if (uplo == 121) # upper
      k = 0;
      for j = 1:N
        m(1:j,j) = A(k + (1:j)).';
        k = k + j;
      endfor
    else
      k = 0;
      for j = 1:N
        m(j:N,j) = A(k + (1:(N-j+1))).';
        k = k + N - j + 1;
      endfor
    endif
  else  # row major
    if (uplo == 121) # upper
      k = 0;
      for j = 1:N
        m(j,j:N) = A(k + (1:(N-j+1)));
        k = k + N - j + 1;
      endfor
    else
      k = 0;
      for j = 1:N
        m(j,1:j) = A(k + (1:j));
        k = k + j;
      endfor
    endif
  endif
  
  if (Diag == 132)  # unit diag
    m = m - diag(diag(m)) + eye(size(m));
  endif
endfunction

function m = tpmatout (order, uplo, Diag, A, N, a)
  if (order == 102) # column major 
    if (uplo == 121) # upper
      k = 0;
      for j = 1:N
        A(k + (1:j)) = a(1:j,j).';
        k = k + j;
      endfor
    else
      k = 0;
      for j = 1:N
        A(k + (1:(N-j+1))) = a(j:N,j).';
        k = k + N - j + 1;
      endfor
    endif
  else  # row major
    if (uplo == 121) # upper
      k = 0;
      for j = 1:N
        A(k + (1:(N-j+1))) = a(j,j:N);
        k = k + N - j + 1;
      endfor
    else
      k = 0;
      for j = 1:N
        A(k + (1:j)) = a(j,1:j);
        k = k + j;
      endfor
    endif
  endif
  
  m = A;
  if (Diag == 132)  # unit diag
    #m = m - diag(diag(m)) + eye(size(m));
  endif
endfunction



function MM = op(M, trans);
  if (trans == 111)
    MM = M ;
  elseif (trans == 112)   # transpose
    MM = M.';
  elseif (trans == 113)   # hermitian conjugate
    MM = M' ;
  endif
endfunction

function m = bandmatrix (order, A, lda, M, N, KL, KU)
  if (order == 102)   # column major
    tmp = reshape(A,lda,length(A)/lda);

    tmp = tmp(1:(KU+KL+1), 1:N);
    m = zeros(M,N);
    for j = 1:min(M+KU,N)
      i0 = min(max(1, j - KU),M);
      i1 = min(M, j + KL);
      k0 = i0 + 1 + KU - j;
      k1 = i1 + 1 + KU - j;
      m(i0:i1, j) = tmp(k0:k1,j);
    endfor

  elseif (order == 101)  # row major
    tmp = reshape(A,lda,length(A)/lda).';

    tmp = tmp(1:M, 1:(KU+KL+1));
    m = zeros(M,N);
    for i = 1:min(M,N+KL)
      j0 = min(max(1, i - KL),N);
      j1 = min(N, i + KU);
      k0 = j0 + 1 + KL - i;
      k1 = j1 + 1 + KL - i;
      m(i, j0:j1) = tmp(i, k0:k1);
    endfor

  endif
endfunction

function a = bmout (order, A, lda, M, N, KL, KU, bm)
  if (order == 102)  # column major
    a = reshape(A, lda, length(A)/lda);

    for j = 1:N
      i0 = max(1, j - KU);
      i1 = min(M, j + KL);
      a((i0:i1)+1+KU-j,j) = bm(i0:i1, j);
    endfor

    a = reshape(a, 1, length(A));
  elseif (order == 101)        # row major
    a = reshape(A,lda,length(A)/lda).';

    for i = 1:M
      j0 = max(1, i - KL);
      j1 = min(M, i + KU);
      a(i,(j0:j1)+1+KL-i) = bm(i,j0:j1);
    endfor

    a = reshape(a.', 1, length(A));
  endif
endfunction


function d = blas_dsdot (N, alpha, X, incX, Y, incY)
  if (N <= 0)
    d = 0;
    return
  endif
  x = vector(X, incX, N);
  y = vector(Y, incY, N);
  d = alpha + (x' * y) ;
endfunction

function d = blas_dot (N, X, incX, Y, incY)
  if (N <= 0)
    d = 0;
    return
  endif
  x = vector(X, incX, N);
  y = vector(Y, incY, N);
  d = (x' * y) ;
endfunction

function d = blas_dotu  (N, X, incX, Y, incY)
  if (N <= 0)
    d = 0;
    return
  endif
  x = vector(X, incX, N);
  y = vector(Y, incY, N);
  d = (x.' * y) ;
endfunction

function d = blas_dotc  (N, X, incX, Y, incY)
  if (N <= 0)
    d = 0;
    return
  endif
  x = vector(X, incX, N);
  y = vector(Y, incY, N);
  d = (x' * y) ;
endfunction

function d = blas_nrm2 (N, X, incX)
  if (N <= 0 || incX <= 0)
    d = 0;
    return;
  endif
  x = vector(X, incX, N);
  d = norm(x);
endfunction

function d = blas_asum (N, X, incX)
  if (N <= 0 || incX <= 0)
    d = 0;
    return;
  endif
  x = vector(X, incX, N);
  d = norm(real(x), 1) + norm(imag(x),1);
endfunction

function k = blas_amax (N, X, incX)
  if (N <= 0 || incX <= 0)
    k = 0;
    return;
  endif
  x = vector(X, incX, N);
  k = 1;
  max = 0;
  for i = 1:length(x)
    v = abs(real(x(i))) + abs(imag(x(i)));
    if (v > max)
      max = v;
      k = i;
    endif
  endfor
  k = k - 1;   #  indices start at 0 in C 
endfunction

function [XX, YY] = blas_swap (N, X, incX, Y, incY)
  if (N <= 0)
    XX=X;
    YY=Y;
    return
  endif
  x = vector(X, incX, N);
  y = vector(Y, incY, N);
  t = x; x = y; y = t;
  XX = vout (X, incX, N, x);
  YY = vout (Y, incY, N, y);
endfunction

function YY = blas_copy (N, X, incX, Y, incY)
  if (N <= 0)
    YY = Y;
    return;
  endif
  x = vector(X, incX, N);
  y = vector(Y, incY, N);
  y = x;
  YY = vout (Y, incY, N, y);
endfunction

function YY = blas_axpy (N, alpha, X, incX, Y, incY)
  if (N <= 0 || alpha == 0.0)
    YY = Y;
    return
  endif
  x = vector(X, incX, N);
  y = vector(Y, incY, N);
  y = alpha * x + y;
  YY = vout (Y, incY, N, y);
endfunction

function [r,z,c,s] = blas_rotg (a, b)
  roe = b;
  if( abs(a) > abs(b) ) 
    roe = a;
  endif
  scale = abs(a) + abs(b);

  if( scale == 0.0 ) 
    c = 1.0;
    s = 0.0;
    r = 0.0;
    z = 0.0;
    return;
  endif
    
  r = scale*sqrt((a/scale)**2 + (b/scale)**2);
  if (roe < 0)
    r = -r;
  endif
  c = a/r;
  s = b/r;
  z = 1.0;
  if( abs(a) > abs(b) ) 
    z = s;
  endif

  if( abs(b) >= abs(a) && c != 0.0) 
    z = 1.0/c;
  endif
endfunction

function [XX, YY] = blas_rot (N, X, incX, Y, incY, c, s)
  if (N <= 0 || (c == 1.0 && s == 0.0))
    XX = X;
    YY = Y;
    return;
  endif
  x = vector(X, incX, N);
  y = vector(Y, incY, N);
  g = [c, s; -s, c];
  M = [x,y] * g';
  x = M(:,1);
  y = M(:,2);
  XX = vout(X, incX, N, x);
  YY = vout(Y, incY, N, y);
endfunction

function [D1, D2, B1, FLAG, H] = blas_rotmg (d1, d2, b1, b2)
  FLAG = NaN;
  H = [ NaN, NaN; NaN, NaN];
  if (d1 < 0)
    H = [0,0;0,0];
    D1 = 0;
    D2 = 0;
    B1 = 0;
    B2 = 0;
    FLAG = -1;
    return;
  elseif (d2 * b2 == 0)
    #H = [1, 0; 0,1];
    B1 = b1;
    D1 = d1;
    D2 = d2;
    FLAG = -2;
    return
  elseif (abs(d1 * b1 * b1) > abs(d2 * b2 * b2))
    FLAG = 0;
    H(1,1) = 1;
    H(1,2) = (d2*b2)/(d1*b1);
    H(2,1) = -(b2/b1);
    H(2,2) = 1;
    u = 1-H(2,1)*H(1,2);
    if (u <= 0)
      H = [0,0;0,0];
      D1 = 0;
      D2 = 0;
      B1 = 0;
      B2 = 0;
      FLAG = -1;
      return;
    endif
    d1 = d1 / u;
    d2 = d2 / u;
    b1 = b1 * u;
    b2 = b2;
  elseif (abs(d1 * b1 * b1) <= abs(d2 * b2 * b2))

    if (d2 * b2 * b2 < 0)
      H = [0,0;0,0];
      D1 = 0;
      D2 = 0;
      B1 = 0;
      B2 = 0;
      FLAG = -1;
      return;
    endif
    FLAG = 1;
    H(1,1) = (d1*b1)/(d2*b2);
    H(1,2) = 1;
    H(2,1) = -1;
    H(2,2) = b1/b2;
    u = 1 + H(1,1) * H(2,2);
    tmp = d1;
    d1 = d2 / u;
    d2 = tmp / u;
    b1 = b2 * u;
    b2 = b1;
  endif

  g = 4096;

  while (d1 != 0 && abs(d1) < 1/g**2)
    FLAG = -1;
    d1 = d1 * g**2;
    H(1,:) = H(1,:) / g;
    b1 = b1 / g;
  endwhile

  while (d1 != 0 && abs(d1) > g**2)
    FLAG = -1;
    d1 = d1 / g**2;
    H(1,:) = H(1,:) * g;
    b1 = b1 * g;
  endwhile

  while (d2 != 0 && abs(d2) < 1/g**2)
    FLAG = -1;
    d2 = d2 * g**2;
    H(2,:) = H(2,:) / g;
  endwhile

  while (d2 != 0 && abs(d2) > g**2)
    FLAG = -1;
    d2 = d2 / g**2;
    H(2,:) = H(2,:) * g;
  endwhile

  if (FLAG == 1)
    H(1,2) = NaN;
    H(2,1) = NaN;
  elseif (FLAG == 0)
    H(1,1) = NaN;
    H(2,2) = NaN;
  endif

  B1 = b1;
  D1 = d1;
  D2 = d2;

endfunction

function [XX, YY] = blas_rotm (N, X, incX, Y, incY, flag, h)
  if (N <= 0)
    XX = X;
    YY = Y;
    return;
  endif
  x = vector(X, incX, N);
  y = vector(Y, incY, N);
  if (flag == -1)
    H = [h(1,1), h(1,2); h(2,1), h(2,2)];
  elseif (flag == 0)
    H = [1, h(1,2); h(2,1), 1];
  elseif (flag == 1)
    H = [h(1,1), 1; -1, h(2,2)];
  elseif (flag == -2)
    H = [1, 0; 0, 1];
  endif  
  M = [x,y] * H';
  x = M(:,1);
  y = M(:,2);
  XX = vout(X, incX, N, x);
  YY = vout(Y, incY, N, y);
endfunction


function XX = blas_scal(N, alpha, X, incX);
  if (N <= 0 || incX <= 0)
    XX = X;
    return;
  endif
  x = vector(X, incX, N);
  x = alpha * x;
  XX = vout(X, incX, N, x);
endfunction


## Level 2 BLAS

function YY = blas_gemv (order, trans, M, N, alpha, A, lda, X, incX, \
                         beta, Y, incY)
  a = matrix (order, A, lda, M, N);
  a = op(a, trans);
  x = vector (X, incX, columns(a));
  y = vector (Y, incY, rows(a));

  y = alpha * a * x + beta * y;
  
  YY = vout(Y, incY, rows(a), y);
endfunction

function YY = blas_gbmv (order, trans, M, N, KL, KU, alpha, A, lda, X, incX, \
                         beta, Y, incY)
  a = bandmatrix (order, A, lda, M, N, KL, KU);
  a = op(a, trans);
  x = vector (X, incX, columns(a));
  y = vector (Y, incY, rows(a));

  y = alpha * a * x + beta * y;
  
  YY = vout(Y, incY, rows(a), y);
endfunction

function XX = blas_trmv (order, uplo, trans, diag, N, A, lda, X, incX)
  a = trmatrix (order, uplo, diag, A, lda, N);
  a = op(a, trans);
  x = vector (X, incX, N);

  y =  a * x ;
  
  XX = vout(X, incX, N, y);
endfunction
 
function XX = blas_tbmv (order, uplo, trans, diag, N, K, A, lda, X, incX)
  a = tbmatrix (order, uplo, diag, A, lda, N, K);
  a = op(a, trans);
  x = vector (X, incX, N);

  y =  a * x ;
  
  XX = vout(X, incX, N, y);
endfunction

function XX = blas_tpmv (order, uplo, trans, diag, N, A, X, incX)
  a = tpmatrix (order, uplo, diag, A, N);
  a = op(a, trans);
  x = vector (X, incX, N);

  y =  a * x ;
  
  XX = vout(X, incX, N, y);
endfunction

function YY = blas_symv (order, uplo, N, alpha, A, lda, X, incX, beta, Y, \
                         incY)
  a = trmatrix (order, uplo, 131, A, lda, N); # nounit
  a = (a + a.' - diag(diag(a)));  # symmetrise
  x = vector (X, incX, N);
  y = vector (Y, incY, N);
  
  y = alpha * a * x + beta * y;
  
  YY = vout(Y, incY, N, y);
endfunction

function YY = blas_hemv (order, uplo, N, alpha, A, lda, X, incX, beta, Y, \
                         incY)
  a = trmatrix (order, uplo, 131, A, lda, N); # nounit
  t = triu(a,1) + tril(a,-1);
  a = diag(real(diag(a))) + t + t';  # make hermitian
  x = vector (X, incX, N);
  y = vector (Y, incY, N);
  
  y = alpha * a * x + beta * y;
  
  YY = vout(Y, incY, N, y);
endfunction

function YY = blas_sbmv (order, uplo, N, K, alpha, A, lda, X, incX, beta, Y, \
                         incY)
  a = tbmatrix (order, uplo, 131, A, lda, N, K); # nounit
  t = triu(a,1) + tril(a,-1);
  a = diag(real(diag(a))) + t + t.';  # make symmetric
  x = vector (X, incX, N);
  y = vector (Y, incY, N);
  
  y = alpha * a * x + beta * y;
  
  YY = vout(Y, incY, N, y);
endfunction


function YY = blas_hbmv (order, uplo, N, K, alpha, A, lda, X, incX, beta, Y, \
                         incY)
  a = tbmatrix (order, uplo, 131, A, lda, N, K); # nounit
  t = triu(a,1) + tril(a,-1);
  a = diag(real(diag(a))) + t + t';  # make hermitian
  x = vector (X, incX, N);
  y = vector (Y, incY, N);
  
  y = alpha * a * x + beta * y;
  
  YY = vout(Y, incY, N, y);
endfunction

function YY = blas_hpmv (order, uplo, N, alpha, A, X, incX, beta, Y, incY)
  a = tpmatrix (order, uplo, 131, A, N); # nounit
  t = triu(a,1) + tril(a,-1);
  a = diag(real(diag(a))) + t + t';  # make hermitian
  x = vector (X, incX, N);
  y = vector (Y, incY, N);
  
  y = alpha * a * x + beta * y;
  
  YY = vout(Y, incY, N, y);
endfunction

function YY = blas_spmv (order, uplo, N, alpha, A, X, incX, beta, Y, incY)
  a = tpmatrix (order, uplo, 131, A, N); # nounit
  t = triu(a,1) + tril(a,-1);
  a = diag(real(diag(a))) + t + t';  # make symmetric
  x = vector (X, incX, N);
  y = vector (Y, incY, N);
  
  y = alpha * a * x + beta * y;
  
  YY = vout(Y, incY, N, y);
endfunction

function XX = blas_trsv (order, uplo, trans, diag, N, A, lda, X, incX)
  a = trmatrix (order, uplo, diag, A, lda, N);
  a = op(a, trans);
  x = vector (X, incX, N);

  y =  a \ x ;
  
  XX = vout(X, incX, N, y);
endfunction

function XX = blas_tbsv (order, uplo, trans, diag, N, K, A, lda, X, incX)
  a = tbmatrix (order, uplo, diag, A, lda, N, K);
  a = op(a, trans);
  x = vector (X, incX, N);

  y =  a \ x ;
  
  XX = vout(X, incX, N, y);
endfunction

function XX = blas_tpsv (order, uplo, trans, diag, N, A, X, incX)
  a = tpmatrix (order, uplo, diag, A, N);
  a = op(a, trans);
  x = vector (X, incX, N);

  y =  a \ x ;
  
  XX = vout(X, incX, N, y);
endfunction

function AA = blas_ger (order, M, N, alpha, X, incX, Y, incY, A, lda)
  a = matrix (order, A, lda, M, N);

  x = vector (X, incX, columns(a));
  y = vector (Y, incY, rows(a));

  if (alpha == 0)
    AA = A;
    return;
  endif

  a = alpha * x * y' + a;
  
  AA = mout(order, A, lda, M, N, a);
endfunction

function AA = blas_geru (order, M, N, alpha, X, incX, Y, incY, A, lda)
  a = matrix (order, A, lda, M, N);

  x = vector (X, incX, columns(a));
  y = vector (Y, incY, rows(a));

  if (alpha == 0)
    AA = A;
    return;
  endif

  a = alpha * x * (y.') + a;
  
  AA = mout(order, A, lda, M, N, a);
endfunction

function AA = blas_gerc (order, M, N, alpha, X, incX, Y, incY, A, lda)
  a = matrix (order, A, lda, M, N);

  x = vector (X, incX, columns(a));
  y = vector (Y, incY, rows(a));

  if (alpha == 0)
    AA = A;
    return;
  endif

  a = alpha * x * y' + a;
  
  AA = mout(order, A, lda, M, N, a);
endfunction

function AA = blas_syr (order, uplo, N, alpha, X, incX, A, lda)
  a = trmatrix (order, uplo, 131, A, lda, N); #nounit
  x = vector (X, incX, N);
  t = triu(a,1) + tril(a,-1);
  a = diag(real(diag(a))) + t + t';  # make symmetric

  if (alpha == 0)
    AA = A;
    return;
  endif

  a = alpha * x * x' + a;
  
  AA = trmatout(order, uplo, 131, A, lda, N, a);
endfunction

function AA = blas_spr (order, uplo, N, alpha, X, incX, A)
  a = tpmatrix (order, uplo, 131, A, N); #nounit
  x = vector (X, incX, N);
  t = triu(a,1) + tril(a,-1);
  a = diag(real(diag(a))) + t + t';  # make symmetric

  if (alpha == 0)
    AA = A;
    return;
  endif

  a = alpha * x * x' + a;
  
  AA = tpmatout(order, uplo, 131, A, N, a);
endfunction

function AA = blas_syr2 (order, uplo, N, alpha, X, incX, Y, incY, A, lda)
  a = trmatrix (order, uplo, 131, A, lda, N); #nounit
  x = vector (X, incX, N);
  y = vector (Y, incY, N);
  t = triu(a,1) + tril(a,-1);
  a = diag(real(diag(a))) + t + t';  # make symmetric

  if (alpha == 0)
    AA = A;
    return;
  endif

  a = alpha * x * y' + alpha * y * x' + a;
  
  AA = trmatout(order, uplo, 131, A, lda, N, a);
endfunction

function AA = blas_spr2 (order, uplo, N, alpha, X, incX, Y, incY, A)
  a = tpmatrix (order, uplo, 131, A, N); #nounit
  x = vector (X, incX, N);
  y = vector (Y, incY, N);
  t = triu(a,1) + tril(a,-1);
  a = diag(real(diag(a))) + t + t';  # make symmetric

  if (alpha == 0)
    AA = A;
    return;
  endif

  a = alpha * x * y' + alpha * y * x' + a;
  
  AA = tpmatout(order, uplo, 131, A, N, a);
endfunction

function AA = blas_her (order, uplo, N, alpha, X, incX, A, lda)
  a = trmatrix (order, uplo, 131, A, lda, N); #nounit
  x = vector (X, incX, N);
  t = triu(a,1) + tril(a,-1);
  a = diag(real(diag(a))) + t + t';  # make symmetric

  if (alpha == 0)
    AA = A;
    return;
  endif

  a = alpha * x * x' + a;
  for i = 1:N
    a(i,i) = real(a(i,i));
  endfor

  AA = trmatout(order, uplo, 131, A, lda, N, a);
endfunction

function AA = blas_hpr (order, uplo, N, alpha, X, incX, A)
  a = tpmatrix (order, uplo, 131, A, N); #nounit
  x = vector (X, incX, N);
  t = triu(a,1) + tril(a,-1);
  a = diag(real(diag(a))) + t + t';  # make symmetric

  if (alpha == 0)
    AA = A;
    return;
  endif

  a = alpha * x * x' + a;
  for i = 1:N
    a(i,i) = real(a(i,i));
  endfor
  
  AA = tpmatout(order, uplo, 131, A, N, a);
endfunction


function AA = blas_her2 (order, uplo, N, alpha, X, incX, Y, incY, A, lda)
  a = trmatrix (order, uplo, 131, A, lda, N); #nounit
  x = vector (X, incX, N);
  y = vector (Y, incY, N);
  t = triu(a,1) + tril(a,-1);
  a = diag(real(diag(a))) + t + t';  # make symmetric

  if (alpha == 0)
    AA = A;
    return;
  endif

  a = alpha * x * y' + conj(alpha) * y * x' + a;
  for i = 1:N
    a(i,i) = real(a(i,i));
  endfor
  
  AA = trmatout(order, uplo, 131, A, lda, N, a);
endfunction

function AA = blas_hpr2 (order, uplo, N, alpha, X, incX, Y, incY, A)
  a = tpmatrix (order, uplo, 131, A, N); #nounit
  x = vector (X, incX, N);
  y = vector (Y, incY, N);
  t = triu(a,1) + tril(a,-1);
  a = diag(real(diag(a))) + t + t';  # make symmetric

  if (alpha == 0)
    AA = A;
    return;
  endif

  a = alpha * x * y' + conj(alpha) * y * x' + a;

  for i = 1:N
    a(i,i) = real(a(i,i));
  endfor
  
  AA = tpmatout(order, uplo, 131, A, N, a);
endfunction


## Level 3 Blas

function CC = blas_gemm (order, transA, transB,  M, N, K, alpha, \
                         A, lda, B, ldb, beta, C, ldc)
  c = matrix (order, C, ldc, M, N);

  if (transA == 111)
    a = matrix (order, A, lda, M, K);
  else
    a = matrix (order, A, lda, K, M);
  endif

  a = op(a, transA);

  if (transB == 111)
    b = matrix (order, B, ldb, K, N);
  else
    b = matrix (order, B, ldb, N, K);
  endif

  b = op(b, transB);


  c = alpha * a * b + beta * c;
  
  CC = mout(order, C, ldc, M, N, c);
endfunction


function CC = blas_symm (order, side, uplo,  M, N, alpha, \
                         A, lda, B, ldb, beta, C, ldc)
  c = matrix (order, C, ldc, M, N);

  if (side == 141)
    a = trmatrix (order, uplo, 131, A, lda, M);
  else
    a = trmatrix (order, uplo, 131, A, lda, N);
  endif

  t = triu(a,1) + tril(a,-1);
  a = t + t.' + diag(diag(a));  # symmetrise

  b = matrix (order, B, ldb, M, N);

  if (side == 141)
    c = alpha * a * b + beta * c;
  else
    c = alpha * b * a + beta * c;
  endif

  CC = mout(order, C, ldc, M, N, c);
endfunction


function CC = blas_hemm (order, side, uplo,  M, N, alpha, \
                         A, lda, B, ldb, beta, C, ldc)
  c = matrix (order, C, ldc, M, N);

  if (side == 141)
    a = trmatrix (order, uplo, 131, A, lda, M);
  else
    a = trmatrix (order, uplo, 131, A, lda, N);
  endif

  t = triu(a,1) + tril(a,-1);
  a = t + t' + diag(real(diag(a)));  # make hermitian

  b = matrix (order, B, ldb, M, N);

  if (side == 141)
    c = alpha * a * b + beta * c;
  else
    c = alpha * b * a + beta * c;
  endif

  CC = mout(order, C, ldc, M, N, c);
endfunction



function CC = blas_syrk (order, uplo, trans, N, K,  alpha, \
                         A, lda, beta, C, ldc)
  c = trmatrix (order, uplo, 131, C, ldc, N);

  if (trans == 111)
    a = matrix (order, A, lda, N, K);
  else
    a = matrix (order, A, lda, K, N);
  endif

  t = triu(c,1) + tril(c,-1);
  c = t + t' + diag(diag(c));  # make symmetric

  if (trans == 111)
    c = alpha * a * a.' + beta * c;
  else
    c = alpha * a.' * a + beta * c;
  endif

  CC = trmatout(order, uplo, 131, C, ldc, N, c);
endfunction


function CC = blas_herk (order, uplo, trans, N, K,  alpha, \
                         A, lda, beta, C, ldc)
  c = trmatrix (order, uplo, 131, C, ldc, N);

  if (beta == 1 && (alpha == 0 || K == 0))
    CC=C;
    return;
  endif

  if (trans == 111)
    a = matrix (order, A, lda, N, K);
  else
    a = matrix (order, A, lda, K, N);
  endif

  t = triu(c,1) + tril(c,-1);
  c = t + t' + diag(real(diag(c)));  # make hermitian

  if (trans == 111)
    c = alpha * a * a' + beta * c;
  else
    c = alpha * a' * a + beta * c;
  endif

  for i = 1:N
    c(i,i) = real(c(i,i));
  endfor

  CC = trmatout(order, uplo, 131, C, ldc, N, c);
endfunction

function CC = blas_syr2k (order, uplo, trans, N, K, alpha, \
                         A, lda, B, ldb, beta, C, ldc)
  c = trmatrix (order, uplo, 131, C, ldc, N);

  t = triu(c,1) + tril(c,-1);
  c = t + t' + diag(diag(c));  # make symmetric

  if (trans == 111)
    a = matrix (order, A, lda, N, K);
    b = matrix (order, B, ldb, N, K);
    c = alpha * (a * b.' + b * a.') + beta * c;
  else
    a = matrix (order, A, lda, K, N);
    b = matrix (order, B, ldb, K, N);
    c = alpha * (a.' * b + b.' * a) + beta * c;
  endif

  CC = trmatout(order, uplo, 131, C, ldc, N, c);
endfunction


function CC = blas_her2k (order, uplo, trans, N, K, alpha, \
                         A, lda, B, ldb, beta, C, ldc)
  c = trmatrix (order, uplo, 131, C, ldc, N);

  if (beta == 1 && (alpha == 0 || K == 0))
    CC=C;
    return;
  endif

  t = triu(c,1) + tril(c,-1);
  c = t + t' + diag(real(diag(c)));  # make hermitian

  if (trans == 111)
    a = matrix (order, A, lda, N, K);
    b = matrix (order, B, ldb, N, K);
    c = alpha * a * b' + alpha' * b * a' + beta * c;
  else
    a = matrix (order, A, lda, K, N);
    b = matrix (order, B, ldb, K, N);
    c = alpha * a' * b + alpha' * b' * a + beta * c;
  endif

  for i = 1:N
    c(i,i) = real(c(i,i));
  endfor

  CC = trmatout(order, uplo, 131, C, ldc, N, c);
endfunction


function BB = blas_trmm (order, side, uplo, trans, diag, M, N, alpha, \
                         A, lda, B, ldb)
  b = matrix (order, B, ldb, M, N);

  if (side == 141)
    a = trmatrix (order, uplo, diag, A, lda, M);
    a = op(a, trans);
    b = alpha * a * b;
  else
    a = trmatrix (order, uplo, diag, A, lda, N);
    a = op(a, trans);
    b = alpha * b * a;
  endif
  
  BB = mout(order, B, ldb, M, N, b);
endfunction

function BB = blas_trsm (order, side, uplo, trans, diag, M, N, alpha, \
                         A, lda, B, ldb)
  b = matrix (order, B, ldb, M, N);

  if (side == 141)
    a = trmatrix (order, uplo, diag, A, lda, M);
    a = op(inv(a), trans);
    b = alpha * a * b;
  else
    a = trmatrix (order, uplo, diag, A, lda, N);
    a = op(inv(a), trans);
    b = alpha * b * a;
  endif
  
  BB = mout(order, B, ldb, M, N, b);
endfunction

######################################################################

                                # testing functions

function begin_file(name,decls)
  global FILE;
  global LIST;
  global HEADER;
  printf("opening %s\n", name) ;
  fprintf(LIST, "  test_%s ();\n", name);
  fprintf(HEADER, "void test_%s (void);\n", name);
  filename = strcat("test_", name, ".c");
  [FILE] = fopen(filename, "w");
  fprintf(FILE,"#include <gsl/gsl_test.h>\n");
  fprintf(FILE,"#include <gsl/gsl_ieee_utils.h>\n");
  fprintf(FILE,"#include <gsl/gsl_math.h>\n");
  fprintf(FILE,"#include <gsl/gsl_cblas.h>\n");
  fprintf(FILE,"\n");
  fprintf(FILE,"#include \"tests.h\"\n");
  fprintf(FILE,"\n");
  fprintf(FILE,"void\n");
  fprintf(FILE,"test_%s (void) {\n", name);
  if (nargin == 1)
    fprintf(FILE,"const double flteps = 1e-4, dbleps = 1e-6;\n");
  else
    # do nothing for now
  endif
endfunction

function end_file()
  global FILE;
  fprintf(FILE, "}\n");
  fclose(FILE);
endfunction

function begin_block(name)
  global FILE;
  fprintf(FILE, "  {\n");
endfunction

function end_block()
  global FILE;
  fprintf(FILE, "  };\n\n\n");
endfunction

function define(S, type, name,x)
  global FILE;
  if (strcmp(type,"scalar"))
    if (S.complex == 0)
      if (nargin == 3)
        fprintf(FILE, "   %s %s;\n", S.precision, name);
      else
        if (strcmp(S.precision,"float"))
          fprintf(FILE, "   %s %s = %#.12gf;\n", S.precision, name, x);
        else 
          fprintf(FILE, "   %s %s = %.12g;\n", S.precision, name, x);
        endif
      endif
    else
      if (nargin == 3)
        fprintf(FILE, "   %s %s[2];\n", S.precision, name);
      else
        if (strcmp(S.precision,"float"))
          fprintf(FILE, "   %s %s[2] = {%#.12gf, %#.12gf};\n", 
                  S.precision, name, real(x), imag(x));
        else
          fprintf(FILE, "   %s %s[2] = {%.12g, %.12g};\n", 
                  S.precision, name, real(x), imag(x));
        endif
      endif
    endif
  elseif (strcmp(type,"int"))
      if (nargin == 3)
        fprintf(FILE, "   %s %s;\n", type, name);
      else
        fprintf(FILE, "   %s %s = %d;\n", type, name, x);
      endif
  elseif (strcmp(type,"vector") || strcmp(type,"matrix"))
    fprintf(FILE, "   %s %s[] = { ", S.precision, name);
    for i = 1:length(x)
      if (i > 1)
        fprintf(FILE, ", ");
      endif
      if ((abs(x(i)) > 1e3 || abs(x(i)) < 1e-3) && abs(x(i)) != 0.0)
        if (strcmp(S.precision,"float"))
          format = "%.6ef";
        else
          format = "%.12e";
        endif
      else
        if (strcmp(S.precision,"float"))
          format = "%#.6gf";
        else
          format = "%#.12g";
        endif
      endif
      if (S.complex == 0)
        fprintf(FILE, format, x(i));
      else
        fprintf(FILE, strcat(format, ", ", format), real(x(i)), imag(x(i)));
      endif
    endfor
    fprintf(FILE, " };\n");
  endif
endfunction

function call (...)
  global FILE;
  fprintf(FILE, "   %s;\n", strcat(all_va_args, ""));
endfunction

function test(S,type,a,b,desc,var)
  global FILE;
  global TEST;
  TEST++;
  desc = strcat(desc, "(case ", int2str(TEST), ")");
  if (strcmp(S.precision, "float"))
    rel = "flteps";
  else
    rel = "dbleps";
  endif

  if (strcmp(type,"scalar"))
    if (S.complex == 0)
      fprintf(FILE, "   gsl_test_rel(%s, %s, %s, \"%s\");\n", a, b, rel, desc);
    else
      fprintf(FILE, "   gsl_test_rel(%s[0], %s[0], %s, \"%s real\");\n", a, b, rel, desc);
      fprintf(FILE, "   gsl_test_rel(%s[1], %s[1], %s, \"%s imag\");\n", a, b, rel, desc);
    endif
  elseif (strcmp(type,"int"))
    fprintf(FILE, "   gsl_test_int(%s, %s, \"%s\");\n", a, b, desc);
  elseif (strcmp(type,"vector") || strcmp(type,"matrix"))
    N = length(var);
    if (S.complex == 0)
#      fprintf(FILE, "   test_rel (%s, %s, %s, \"%s\");\n", N, a, b, rel, desc);
       fprintf(FILE, "   {\n");
       fprintf(FILE, "     int i;\n");
       fprintf(FILE, "     for (i = 0; i < %d; i++) {\n", N);
       fprintf(FILE, "       gsl_test_rel(%s[i], %s[i], %s, \"%s\");\n", a, b, \
              rel, desc);
       fprintf(FILE, "     }\n");
       fprintf(FILE, "   };\n");
    else
#      fprintf(FILE, "   test_zrel (%s, %s, %s, \"%s\");\n", N, a, b, rel, desc);
       fprintf(FILE, "   {\n");
       fprintf(FILE, "     int i;\n");
       fprintf(FILE, "     for (i = 0; i < %d; i++) {\n", N);
       fprintf(FILE, "       gsl_test_rel(%s[2*i], %s[2*i], %s, \"%s real\");\n", a, b, rel, desc);
       fprintf(FILE, "       gsl_test_rel(%s[2*i+1], %s[2*i+1], %s, \"%s imag\");\n", a, b, \
              rel, desc);
       fprintf(FILE, "     };\n");
       fprintf(FILE, "   };\n");
    endif
  endif
  
endfunction


######################################################################

function test_sdsdot (S, fn, N, alpha, X, incX, Y, incY)
  begin_block(fn);
  define(S, "int", "N", N);
  define(S, "scalar", "alpha", alpha);
  define(S, "vector", "X", X);
  define(S, "vector", "Y", Y);
  define(S, "int", "incX", incX);
  define(S, "int", "incY", incY);
  define(S, "scalar", "expected", feval(strcat("blas_", fn), 
                                        N, alpha, X, incX, Y, incY));

  define(S, "scalar", "f");
  if (S.complex == 0)
    call("f = cblas_sdsdot (N, alpha, X, incX, Y, incY)");
  else
    call("cblas_sdsdot (N, alpha, X, incX, Y, incY, &f)");
  endif
  test(S, "scalar", "f", "expected", "sdsdot");
  end_block();
endfunction

function test_dot (S, fn, N, X, incX, Y, incY)
  begin_block(fn);
  define(S, "int", "N", N);
  define(S, "vector", "X", X);
  define(S, "vector", "Y", Y);
  define(S, "int", "incX", incX);
  define(S, "int", "incY", incY);
  define(S, "scalar", "expected", feval(strcat("blas_", fn), 
                                        N, X, incX, Y, incY));

  define(S, "scalar", "f");
  if (S.complex == 0)
    call("f = cblas_", S.prefix, fn, "(N, X, incX, Y, incY)");
  else
    call("cblas_", S.prefix, fn, "_sub(N, X, incX, Y, incY, &f)");
  endif
  test(S, "scalar", "f", "expected", strcat(S.prefix, fn));
  end_block();
endfunction

function test_nrm2 (S, fn, N, X, incX)
  begin_block(fn);
  define(S, "int", "N", N);
  define(S, "vector", "X", X);
  define(S, "int", "incX", incX);
  T = S; T.complex = 0;
  define(T, "scalar", "expected", feval(strcat("blas_", fn), N, X, incX));
  define(T, "scalar", "f");
  if (strcmp(S.prefix,"c"))
    S.prefix = "sc";
  elseif (strcmp(S.prefix,"z"))
    S.prefix = "dz";
  endif
  call("f = cblas_", S.prefix, fn, "(N, X, incX)");
  test(T, "scalar", "f", "expected", strcat(S.prefix, fn));
  end_block();
endfunction

function test_asum (S, fn, N, X, incX)
  begin_block(fn);
  define(S, "int", "N", N);
  define(S, "vector", "X", X);
  define(S, "int", "incX", incX);
  T = S; T.complex = 0;
  define(T, "scalar", "expected", feval(strcat("blas_", fn), N, X, incX));
  define(T, "scalar", "f");
  if (strcmp(S.prefix,"c"))
    S.prefix = "sc";
  elseif (strcmp(S.prefix,"z"))
    S.prefix = "dz";
  endif
  call("f = cblas_", S.prefix, fn, "(N, X, incX)");
  test(T, "scalar", "f", "expected", strcat(S.prefix, fn));
  end_block();
endfunction

function test_amax (S, fn, N, X, incX)
  begin_block(fn);
  define(S, "int", "N", N);
  define(S, "vector", "X", X);
  define(S, "int", "incX", incX);
  T = S; T.complex = 0;
  define(T, "int", "expected", feval(strcat("blas_", fn), N, X, incX));
  define(T, "int", "k");
  call("k = cblas_i", S.prefix, fn, "(N, X, incX)");
  test(T, "int", "k", "expected", strcat(S.prefix, fn));
  end_block();
endfunction

function test_axpy (S, fn, N, alpha, X, incX, Y, incY)
  begin_block(fn);
  define(S, "int", "N", N);
  define(S, "scalar", "alpha", alpha);
  define(S, "vector", "X", X);
  define(S, "int", "incX", incX);
  define(S, "vector", "Y", Y);
  define(S, "int", "incY", incY);
  define(S, "vector", "expected", feval(strcat("blas_", fn), N, alpha, \
                                        X, incX, Y, incY));
  call("cblas_", S.prefix, fn, "(N, alpha, X, incX, Y, incY)");
  test(S, "vector", "Y", "expected", strcat(S.prefix, fn), Y);
  end_block();
endfunction

function test_copy (S, fn, N, X, incX, Y, incY)
  begin_block(fn);
  define(S, "int", "N", N);
  define(S, "vector", "X", X);
  define(S, "int", "incX", incX);
  define(S, "vector", "Y", Y);
  define(S, "int", "incY", incY);
  define(S, "vector", "expected", feval(strcat("blas_", fn), N, \
                                        X, incX, Y, incY));
  call("cblas_", S.prefix, fn, "(N, X, incX, Y, incY)");
  test(S, "vector", "Y", "expected", strcat(S.prefix, fn), Y);
  end_block();
endfunction

function test_swap (S, fn, N, X, incX, Y, incY)
  begin_block(fn);
  define(S, "int", "N", N);
  define(S, "vector", "X", X);
  define(S, "int", "incX", incX);
  define(S, "vector", "Y", Y);
  define(S, "int", "incY", incY);
  [XX, YY] = feval(strcat("blas_", fn), N, X, incX, Y, incY);
  define(S, "vector", "expected1", XX);
  define(S, "vector", "expected2", YY);
  call("cblas_", S.prefix, fn, "(N, X, incX, Y, incY)");
  test(S, "vector", "X", "expected1", strcat(S.prefix, fn), X);
  test(S, "vector", "Y", "expected2", strcat(S.prefix, fn), Y);
  end_block();
endfunction

function test_scal (S, fn, N, alpha, X, incX)
  begin_block(fn);
  define(S, "int", "N", N);
  define(S, "scalar", "alpha", alpha);
  define(S, "vector", "X", X);
  define(S, "int", "incX", incX);
  define(S, "vector", "expected", feval(strcat("blas_", fn), N, alpha, X, incX));
  call("cblas_", S.prefix, fn, "(N, alpha, X, incX)");
  test(S, "vector", "X", "expected", strcat(S.prefix, fn), X);
  end_block();
endfunction

function test_rotg (S, fn, a, b)
  begin_block(fn);
  define(S, "scalar", "a", a);
  define(S, "scalar", "b", b);
  define(S, "scalar", "c");
  define(S, "scalar", "s");
  [r,z,c,s] = feval(strcat("blas_", fn), a, b);
  define(S, "scalar", "r_expected", r);
  define(S, "scalar", "z_expected", z);
  define(S, "scalar", "c_expected", c);
  define(S, "scalar", "s_expected", s);
  call("cblas_", S.prefix, fn, "(&a, &b, &c, &s)");
  test(S, "scalar", "a", "r_expected", strcat(S.prefix, fn));
  test(S, "scalar", "b", "z_expected", strcat(S.prefix, fn));
  test(S, "scalar", "c", "c_expected", strcat(S.prefix, fn));
  test(S, "scalar", "s", "s_expected", strcat(S.prefix, fn));
  end_block();
endfunction

function test_rot (S, fn, N, X, incX, Y, incY, c, s)
  begin_block(fn);
  define(S, "int", "N", N);
  define(S, "scalar", "c", c);
  define(S, "scalar", "s", s);
  define(S, "vector", "X", X);
  define(S, "int", "incX", incX);
  define(S, "vector", "Y", Y);
  define(S, "int", "incY", incY);
  [XX, YY] = feval(strcat("blas_", fn), N, X, incX, Y, incY, c, s);
  define(S, "vector", "x_expected", XX);
  define(S, "vector", "y_expected", YY);
  call("cblas_", S.prefix, fn, "(N, X, incX, Y, incY, c, s)");
  test(S, "vector", "X", "x_expected", strcat(S.prefix, fn), X);
  test(S, "vector", "Y", "y_expected", strcat(S.prefix, fn), Y);
  end_block();
endfunction


function test_rotmg (S, fn, d1, d2, b1, b2)
  begin_block(fn);
  v0 = [-999.0, -999.1, -999.2, -999.3, -999.4];
  define(S, "scalar", "d1", d1);
  define(S, "scalar", "d2", d2);
  define(S, "scalar", "b1", b1);
  define(S, "scalar", "b2", b2);
  define(S, "vector", "h", v0);
  [D1, D2, B1, FLAG, H] = feval(strcat("blas_", fn), d1, d2, b1, b2);
  define(S, "scalar", "d1_expected", D1);
  define(S, "scalar", "d2_expected", D2);
  define(S, "scalar", "b1_expected", B1);
  define(S, "scalar", "h0_expected", FLAG);
  if (!isnan(H(1,1)))
    define(S, "scalar", "h11_expected", H(1,1));
  else
    define(S, "scalar", "h11_expected", v0(2));
  endif
  if (!isnan(H(2,1)))
    define(S, "scalar", "h21_expected", H(2,1));
  else
    define(S, "scalar", "h21_expected", v0(3));
  endif
  if (!isnan(H(1,2)))
    define(S, "scalar", "h12_expected", H(1,2));
  else
    define(S, "scalar", "h12_expected", v0(4));
  endif
  if (!isnan(H(2,2)))
    define(S, "scalar", "h22_expected", H(2,2));
  else
    define(S, "scalar", "h22_expected", v0(5));
  endif
  call("cblas_", S.prefix, fn, "(&d1, &d2, &b1, b2, h)");
  test(S, "scalar", "d1", "d1_expected", strcat(S.prefix, fn));
  test(S, "scalar", "d2", "d2_expected", strcat(S.prefix, fn));
  test(S, "scalar", "b1", "b1_expected", strcat(S.prefix, fn));
  test(S, "scalar", "h[0]", "h0_expected", strcat(S.prefix, fn));
  test(S, "scalar", "h[1]", "h11_expected", strcat(S.prefix, fn));
  test(S, "scalar", "h[2]", "h21_expected", strcat(S.prefix, fn));
  test(S, "scalar", "h[3]", "h12_expected", strcat(S.prefix, fn));
  test(S, "scalar", "h[4]", "h22_expected", strcat(S.prefix, fn));

  end_block();
endfunction

function test_rotm (S, fn, N, X, incX, Y, incY, flag, h)
  begin_block(fn);
  define(S, "int", "N", N);
  define(S, "vector", "h", [flag, h(1,1), h(2,1), h(1,2), h(2,2)]);
  define(S, "vector", "X", X);
  define(S, "int", "incX", incX);
  define(S, "vector", "Y", Y);
  define(S, "int", "incY", incY);
  [XX, YY] = feval(strcat("blas_", fn), N, X, incX, Y, incY, flag, h);
  define(S, "vector", "x_expected", XX);
  define(S, "vector", "y_expected", YY);
  call("cblas_", S.prefix, fn, "(N, X, incX, Y, incY, h)");
  test(S, "vector", "X", "x_expected", strcat(S.prefix, fn), X);
  test(S, "vector", "Y", "y_expected", strcat(S.prefix, fn), Y);
  end_block();
endfunction


function test_gemv (S, fn, order, trans, M, N, alpha, A, lda, X, incX, \
                    beta, Y, incY)
  begin_block(fn);
  define(S, "int", "order", order);
  define(S, "int", "trans", trans);
  define(S, "int", "M", M);
  define(S, "int", "N", N);
  define(S, "int", "lda", lda);
  define(S, "scalar", "alpha", alpha);
  define(S, "scalar", "beta", beta);
  define(S, "matrix", "A", A);
  define(S, "vector", "X", X);
  define(S, "int", "incX", incX);
  define(S, "vector", "Y", Y);
  define(S, "int", "incY", incY);
  YY = feval(strcat("blas_", fn), order, trans, M, N, alpha, A, lda, \
             X, incX, beta, Y, incY);
  define(S, "vector", "y_expected", YY);
  call("cblas_", S.prefix, fn, "(order, trans, M, N, alpha, A, lda, X, \
                                 incX, beta, Y, incY)");
  test(S, "vector", "Y", "y_expected", strcat(S.prefix, fn), Y);
  end_block();
endfunction

function test_gbmv (S, fn, order, trans, M, N, KL, KU, alpha, A, lda, \
                    X, incX, beta, Y, incY)
  begin_block(fn);
  define(S, "int", "order", order);
  define(S, "int", "trans", trans);
  define(S, "int", "M", M);
  define(S, "int", "N", N);
  define(S, "int", "KL", KL);
  define(S, "int", "KU", KU);
  define(S, "int", "lda", lda);
  define(S, "scalar", "alpha", alpha);
  define(S, "scalar", "beta", beta);
  define(S, "matrix", "A", A);
  define(S, "vector", "X", X);
  define(S, "int", "incX", incX);
  define(S, "vector", "Y", Y);
  define(S, "int", "incY", incY);
  YY = feval(strcat("blas_", fn), order, trans, M, N, KU, KL, alpha, A, lda, \
             X, incX, beta, Y, incY);
  define(S, "vector", "y_expected", YY);
  call("cblas_", S.prefix, fn, "(order, trans, M, N, KU, KL, alpha, A, \
                                 lda, X, incX, beta, Y, incY)");
  test(S, "vector", "Y", "y_expected", strcat(S.prefix, fn), Y);
  end_block();
endfunction


function test_trmv (S, fn, order, uplo, trans, diag, N, A, lda, X, incX)
  begin_block(fn);
  define(S, "int", "order", order);
  define(S, "int", "trans", trans);
  define(S, "int", "uplo", uplo);
  define(S, "int", "diag", diag);
  define(S, "int", "N", N);
  define(S, "int", "lda", lda);
  define(S, "matrix", "A", A);
  define(S, "vector", "X", X);
  define(S, "int", "incX", incX);

  XX = feval(strcat("blas_", fn), order, uplo, trans, diag, N, A, lda, X, incX);
  define(S, "vector", "x_expected", XX);
  call("cblas_", S.prefix, fn, "(order, uplo, trans, diag, N, A, lda, X, incX)");
  test(S, "vector", "X", "x_expected", strcat(S.prefix, fn), X);
  end_block();
endfunction

function test_tbmv (S, fn, order, uplo, trans, diag, N, K, A, lda, X, incX)
  begin_block(fn);
  define(S, "int", "order", order);
  define(S, "int", "trans", trans);
  define(S, "int", "uplo", uplo);
  define(S, "int", "diag", diag);
  define(S, "int", "N", N);
  define(S, "int", "K", K);
  define(S, "int", "lda", lda);
  define(S, "matrix", "A", A);
  define(S, "vector", "X", X);
  define(S, "int", "incX", incX);

  XX = feval(strcat("blas_", fn), order, uplo, trans, diag, N, K, A, lda, X, incX);
  define(S, "vector", "x_expected", XX);
  call("cblas_", S.prefix, fn, "(order, uplo, trans, diag, N, K, A, lda, X, incX)");
  test(S, "vector", "X", "x_expected", strcat(S.prefix, fn), X);
  end_block();
endfunction


function test_tpmv (S, fn, order, uplo, trans, diag, N, A, X, incX)
  begin_block(fn);
  define(S, "int", "order", order);
  define(S, "int", "trans", trans);
  define(S, "int", "uplo", uplo);
  define(S, "int", "diag", diag);
  define(S, "int", "N", N);
  define(S, "matrix", "A", A);
  define(S, "vector", "X", X);
  define(S, "int", "incX", incX);

  XX = feval(strcat("blas_", fn), order, uplo, trans, diag, N, A, X, incX);
  define(S, "vector", "x_expected", XX);
  call("cblas_", S.prefix, fn, "(order, uplo, trans, diag, N, A, X, incX)");
  test(S, "vector", "X", "x_expected", strcat(S.prefix, fn), X);
  end_block();
endfunction

function test_hesymv (S, fn, order, uplo, N, alpha, A, lda,  X, incX, \
                    beta, Y, incY)
  begin_block(fn);
  define(S, "int", "order", order);
  define(S, "int", "uplo", uplo);
  define(S, "scalar", "alpha", alpha);
  define(S, "scalar", "beta", beta);
  define(S, "int", "N", N);
  define(S, "int", "lda", lda);
  define(S, "matrix", "A", A);
  define(S, "vector", "X", X);
  define(S, "int", "incX", incX);
  define(S, "vector", "Y", Y);
  define(S, "int", "incY", incY);

  YY = feval(strcat("blas_", fn), order, uplo, N, alpha, A, lda, X, incX, \
             beta, Y, incY);
  define(S, "vector", "y_expected", YY);
  call("cblas_", S.prefix, fn, "(order, uplo, N, alpha, A, lda, X, incX, beta, Y, incY)");
  test(S, "vector", "Y", "y_expected", strcat(S.prefix, fn), Y);
  end_block();
endfunction

function test_hbsbmv (S, fn, order, uplo, N, k, alpha, A, lda,  X, incX, \
                    beta, Y, incY)
  begin_block(fn);
  define(S, "int", "order", order);
  define(S, "int", "uplo", uplo);
  define(S, "scalar", "alpha", alpha);
  define(S, "scalar", "beta", beta);
  define(S, "int", "N", N);
  define(S, "int", "k", k);
  define(S, "int", "lda", lda);
  define(S, "matrix", "A", A);
  define(S, "vector", "X", X);
  define(S, "int", "incX", incX);
  define(S, "vector", "Y", Y);
  define(S, "int", "incY", incY);

  YY = feval(strcat("blas_", fn), order, uplo, N, k, alpha, A, lda, X, incX, \
             beta, Y, incY);
  define(S, "vector", "y_expected", YY);
  call("cblas_", S.prefix, fn, "(order, uplo, N, k, alpha, A, lda, X, incX, beta, Y, incY)");
  test(S, "vector", "Y", "y_expected", strcat(S.prefix, fn), Y);
  end_block();
endfunction

function test_hpspmv (S, fn, order, uplo, N, alpha, A,  X, incX, \
                    beta, Y, incY)
  begin_block(fn);
  define(S, "int", "order", order);
  define(S, "int", "uplo", uplo);
  define(S, "scalar", "alpha", alpha);
  define(S, "scalar", "beta", beta);
  define(S, "int", "N", N);
  define(S, "matrix", "A", A);
  define(S, "vector", "X", X);
  define(S, "int", "incX", incX);
  define(S, "vector", "Y", Y);
  define(S, "int", "incY", incY);

  YY = feval(strcat("blas_", fn), order, uplo, N, alpha, A, X, incX, \
             beta, Y, incY);
  define(S, "vector", "y_expected", YY);
  call("cblas_", S.prefix, fn, "(order, uplo, N, alpha, A, X, incX, beta, Y, incY)");
  test(S, "vector", "Y", "y_expected", strcat(S.prefix, fn), Y);
  end_block();
endfunction


function test_ger (S, fn, order, M, N, alpha, X, incX, \
                   Y, incY, A, lda)
  begin_block(fn);
  define(S, "int", "order", order);
  define(S, "int", "M", M);
  define(S, "int", "N", N);
  define(S, "int", "lda", lda);
  define(S, "scalar", "alpha", alpha);
  define(S, "matrix", "A", A);
  define(S, "vector", "X", X);
  define(S, "int", "incX", incX);
  define(S, "vector", "Y", Y);
  define(S, "int", "incY", incY);
  AA = feval(strcat("blas_", fn), order, M, N, alpha, X, incX, Y, incY, A, lda);
  define(S, "matrix", "A_expected", AA);
  call("cblas_", S.prefix, fn, "(order, M, N, alpha, X, incX, Y, incY, A, lda)");
  test(S, "vector", "A", "A_expected", strcat(S.prefix, fn), A);
  end_block();
endfunction

function test_syr (S, fn, order, uplo, N, alpha, X, incX, A, lda)
  begin_block(fn);
  define(S, "int", "order", order);
  define(S, "int", "uplo", uplo);
  define(S, "int", "N", N);
  define(S, "int", "lda", lda);
  define(S, "scalar", "alpha", alpha);
  define(S, "matrix", "A", A);
  define(S, "vector", "X", X);
  define(S, "int", "incX", incX);
  AA = feval(strcat("blas_", fn), order, uplo, N, alpha, X, incX, A, lda);
  define(S, "matrix", "A_expected", AA);
  call("cblas_", S.prefix, fn, "(order, uplo, N, alpha, X, incX, A, lda)");
  test(S, "vector", "A", "A_expected", strcat(S.prefix, fn), A);
  end_block();
endfunction

function test_spr (S, fn, order, uplo, N, alpha, X, incX, Ap)
  begin_block(fn);
  define(S, "int", "order", order);
  define(S, "int", "uplo", uplo);
  define(S, "int", "N", N);
  define(S, "scalar", "alpha", alpha);
  define(S, "matrix", "Ap", Ap);
  define(S, "vector", "X", X);
  define(S, "int", "incX", incX);
  AA = feval(strcat("blas_", fn), order, uplo, N, alpha, X, incX, Ap);
  define(S, "matrix", "Ap_expected", AA);
  call("cblas_", S.prefix, fn, "(order, uplo, N, alpha, X, incX, Ap)");
  test(S, "vector", "Ap", "Ap_expected", strcat(S.prefix, fn), Ap);
  end_block();
endfunction

function test_syr2 (S, fn, order, uplo, N, alpha, X, incX, Y, incY, A, lda)
  begin_block(fn);
  define(S, "int", "order", order);
  define(S, "int", "uplo", uplo);
  define(S, "int", "N", N);
  define(S, "int", "lda", lda);
  define(S, "scalar", "alpha", alpha);
  define(S, "matrix", "A", A);
  define(S, "vector", "X", X);
  define(S, "int", "incX", incX);
  define(S, "vector", "Y", Y);
  define(S, "int", "incY", incY);
  AA = feval(strcat("blas_", fn), order, uplo, N, alpha, X, incX, Y, incY, A, lda);
  define(S, "matrix", "A_expected", AA);
  call("cblas_", S.prefix, fn, "(order, uplo, N, alpha, X, incX, Y, incY, A, lda)");
  test(S, "vector", "A", "A_expected", strcat(S.prefix, fn), A);
  end_block();
endfunction

function test_spr2 (S, fn, order, uplo, N, alpha, X, incX, Y, incY, Ap)
  begin_block(fn);
  define(S, "int", "order", order);
  define(S, "int", "uplo", uplo);
  define(S, "int", "N", N);
  define(S, "scalar", "alpha", alpha);
  define(S, "matrix", "Ap", Ap);
  define(S, "vector", "X", X);
  define(S, "int", "incX", incX);
  define(S, "vector", "Y", Y);
  define(S, "int", "incY", incY);
  AA = feval(strcat("blas_", fn), order, uplo, N, alpha, X, incX, Y, \
             incY, Ap);
  define(S, "matrix", "Ap_expected", AA);
  call("cblas_", S.prefix, fn, "(order, uplo, N, alpha, X, incX, Y, incY, Ap)");
  test(S, "vector", "Ap", "Ap_expected", strcat(S.prefix, fn), Ap);
  end_block();
endfunction


function test_her (S, fn, order, uplo, N, alpha, X, incX, A, lda)
  begin_block(fn);
  define(S, "int", "order", order);
  define(S, "int", "uplo", uplo);
  define(S, "int", "N", N);
  define(S, "int", "lda", lda);
  T = S ; T.complex = 0;
  define(T, "scalar", "alpha", alpha);
  define(S, "matrix", "A", A);
  define(S, "vector", "X", X);
  define(S, "int", "incX", incX);
  AA = feval(strcat("blas_", fn), order, uplo, N, alpha, X, incX, A, lda);
  define(S, "matrix", "A_expected", AA);
  call("cblas_", S.prefix, fn, "(order, uplo, N, alpha, X, incX, A, lda)");
  test(S, "vector", "A", "A_expected", strcat(S.prefix, fn), A);
  end_block();
endfunction

function test_hpr (S, fn, order, uplo, N, alpha, X, incX, Ap)
  begin_block(fn);
  define(S, "int", "order", order);
  define(S, "int", "uplo", uplo);
  define(S, "int", "N", N);
  T = S ; T.complex = 0;
  define(T, "scalar", "alpha", alpha);
  define(S, "matrix", "Ap", Ap);
  define(S, "vector", "X", X);
  define(S, "int", "incX", incX);
  AA = feval(strcat("blas_", fn), order, uplo, N, alpha, X, incX, Ap);
  define(S, "matrix", "Ap_expected", AA);
  call("cblas_", S.prefix, fn, "(order, uplo, N, alpha, X, incX, Ap)");
  test(S, "vector", "Ap", "Ap_expected", strcat(S.prefix, fn), Ap);
  end_block();
endfunction

function test_gemm (S, fn, order, transA, transB, M, N, K, alpha, A, \
                    lda, B, ldb, beta, C, ldc)
  begin_block(fn);
  define(S, "int", "order", order);
  define(S, "int", "transA", transA);
  define(S, "int", "transB", transB);
  define(S, "int", "M", M);
  define(S, "int", "N", N);
  define(S, "int", "K", K);
  define(S, "scalar", "alpha", alpha);
  define(S, "scalar", "beta", beta);
  define(S, "matrix", "A", A);
  define(S, "int", "lda", lda);
  define(S, "matrix", "B", B);
  define(S, "int", "ldb", ldb);
  define(S, "matrix", "C", C);
  define(S, "int", "ldc", ldc);

  CC = feval(strcat("blas_", fn), order, transA, transB, M, N, K, \
             alpha, A, lda, B, ldb, beta, C, ldc);
  define(S, "matrix", "C_expected", CC);
  call("cblas_", S.prefix, fn, "(order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)");
  test(S, "matrix", "C", "C_expected", strcat(S.prefix, fn), C);
  end_block();
endfunction


function test_hesymm (S, fn, order, side, uplo, M, N, alpha, A, \
                    lda, B, ldb, beta, C, ldc)
  begin_block(fn);
  define(S, "int", "order", order);
  define(S, "int", "side", side);
  define(S, "int", "uplo", uplo);
  define(S, "int", "M", M);
  define(S, "int", "N", N);
  define(S, "scalar", "alpha", alpha);
  define(S, "scalar", "beta", beta);
  define(S, "matrix", "A", A);
  define(S, "int", "lda", lda);
  define(S, "matrix", "B", B);
  define(S, "int", "ldb", ldb);
  define(S, "matrix", "C", C);
  define(S, "int", "ldc", ldc);

  CC = feval(strcat("blas_", fn), order, side, uplo, M, N, \
             alpha, A, lda, B, ldb, beta, C, ldc);
  define(S, "matrix", "C_expected", CC);
  call("cblas_", S.prefix, fn, "(order, side, uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc)");
  test(S, "matrix", "C", "C_expected", strcat(S.prefix, fn), C);
  end_block();
endfunction


function test_syrk (S, fn, order, uplo, trans, N, K, alpha, A, \
                    lda, beta, C, ldc)
  begin_block(fn);
  define(S, "int", "order", order);
  define(S, "int", "uplo", uplo);
  define(S, "int", "trans", trans);
  define(S, "int", "N", N);
  define(S, "int", "K", K);
  define(S, "scalar", "alpha", alpha);
  define(S, "scalar", "beta", beta);
  define(S, "matrix", "A", A);
  define(S, "int", "lda", lda);
  define(S, "matrix", "C", C);
  define(S, "int", "ldc", ldc);

  CC = feval(strcat("blas_", fn), order, uplo, trans, N, K, \
             alpha, A, lda, beta, C, ldc);
  define(S, "matrix", "C_expected", CC);
  call("cblas_", S.prefix, fn, "(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc)");
  test(S, "matrix", "C", "C_expected", strcat(S.prefix, fn), C);
  end_block();
endfunction

function test_herk (S, fn, order, uplo, trans, N, K, alpha, A, \
                    lda, beta, C, ldc)
  begin_block(fn);
  define(S, "int", "order", order);
  define(S, "int", "uplo", uplo);
  define(S, "int", "trans", trans);
  define(S, "int", "N", N);
  define(S, "int", "K", K);
  T = S ; T.complex = 0;
  define(T, "scalar", "alpha", alpha);
  define(T, "scalar", "beta", beta);
  define(S, "matrix", "A", A);
  define(S, "int", "lda", lda);
  define(S, "matrix", "C", C);
  define(S, "int", "ldc", ldc);

  CC = feval(strcat("blas_", fn), order, uplo, trans, N, K, \
             alpha, A, lda, beta, C, ldc);
  define(S, "matrix", "C_expected", CC);
  call("cblas_", S.prefix, fn, "(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc)");
  test(S, "matrix", "C", "C_expected", strcat(S.prefix, fn), C);
  end_block();
endfunction


function test_syr2k (S, fn, order, uplo, trans,  N, K, alpha, A, \
                    lda, B, ldb, beta, C, ldc)
  begin_block(fn);
  define(S, "int", "order", order);
  define(S, "int", "uplo", uplo);
  define(S, "int", "trans", trans);
  define(S, "int", "N", N);
  define(S, "int", "K", K);
  define(S, "scalar", "alpha", alpha);
  define(S, "scalar", "beta", beta);
  define(S, "matrix", "A", A);
  define(S, "int", "lda", lda);
  define(S, "matrix", "B", B);
  define(S, "int", "ldb", ldb);
  define(S, "matrix", "C", C);
  define(S, "int", "ldc", ldc);

  CC = feval(strcat("blas_", fn), order, uplo, trans, N, K, \
             alpha, A, lda, B, ldb, beta, C, ldc);
  define(S, "matrix", "C_expected", CC);
  call("cblas_", S.prefix, fn, "(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc)");
  test(S, "matrix", "C", "C_expected", strcat(S.prefix, fn), C);
  end_block();
endfunction

function test_her2k (S, fn, order, uplo, trans,  N, K, alpha, A, \
                    lda, B, ldb, beta, C, ldc)
  begin_block(fn);
  define(S, "int", "order", order);
  define(S, "int", "uplo", uplo);
  define(S, "int", "trans", trans);
  define(S, "int", "N", N);
  define(S, "int", "K", K);
  define(S, "scalar", "alpha", alpha);
  T=S; T.complex = 0;
  define(T, "scalar", "beta", beta);
  define(S, "matrix", "A", A);
  define(S, "int", "lda", lda);
  define(S, "matrix", "B", B);
  define(S, "int", "ldb", ldb);
  define(S, "matrix", "C", C);
  define(S, "int", "ldc", ldc);

  CC = feval(strcat("blas_", fn), order, uplo, trans, N, K, \
             alpha, A, lda, B, ldb, beta, C, ldc);
  define(S, "matrix", "C_expected", CC);
  call("cblas_", S.prefix, fn, "(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc)");
  test(S, "matrix", "C", "C_expected", strcat(S.prefix, fn), C);
  end_block();
endfunction

function test_trmm (S, fn, order, side, uplo, trans, diag, M, N, alpha, A, \
                    lda, B, ldb)
  begin_block(fn);
  define(S, "int", "order", order);
  define(S, "int", "side", side);
  define(S, "int", "uplo", uplo);
  define(S, "int", "trans", trans);
  define(S, "int", "diag", diag);
  define(S, "int", "M", M);
  define(S, "int", "N", N);
  define(S, "scalar", "alpha", alpha);
  define(S, "matrix", "A", A);
  define(S, "int", "lda", lda);
  define(S, "matrix", "B", B);
  define(S, "int", "ldb", ldb);

  BB = feval(strcat("blas_", fn), order, side, uplo, trans, diag, M, N, \
             alpha, A, lda, B, ldb);
  define(S, "matrix", "B_expected", BB);
  call("cblas_", S.prefix, fn, "(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb)");
  test(S, "matrix", "B", "B_expected", strcat(S.prefix, fn), B);
  end_block();
endfunction

######################################################################

s=1;d=2;c=3;z=4;

testcases=[1:3];
testcases2=[3];
testcases3=[3];

begin_file("dot");
for j = testcases
  for i = [s]
    S = context(i);
    T = test_vectors(S, j);
    for alpha = [0.0, 0.1, 1.0]
      test_sdsdot (S, "dsdot", T.n, alpha, T.v1, T.s1, T.v2, T.s2);
    endfor
  endfor
endfor

for j = testcases
  for i = [s,d,c,z]
    S = context(i);
    T = test_vectors(S, j);
    if (S.complex == 0)
      test_dot (S, "dot", T.n, T.v1, T.s1, T.v2, T.s2);
    else
      test_dot (S, "dotu", T.n, T.v1, T.s1, T.v2, T.s2);
      test_dot (S, "dotc", T.n, T.v1, T.s1, T.v2, T.s2);
    endif
  endfor
endfor
end_file();

begin_file("nrm2");
for j = testcases
  for i = [s,d,c,z]
    S = context(i);
    T = test_vector(S, j);
    test_nrm2 (S, "nrm2", T.n, T.v, T.s);
  endfor
endfor
end_file();

begin_file("asum");
for j = testcases
  for i = [s,d,c,z]
    S = context(i);
    T = test_vector(S, j);
    test_nrm2 (S, "asum", T.n, T.v, T.s);
  endfor
endfor
end_file();

begin_file("amax", "int");
for j = testcases
  for i = [s,d,c,z]
    S = context(i);
    T = test_vector(S, j);
    test_amax (S, "amax", T.n, T.v, T.s);
  endfor
endfor
end_file();

begin_file("axpy");
for j = testcases
  for i = [s,d,c,z]
    S = context(i);
    T = test_vectors(S, j);
    for alpha = coeff(S)
      test_axpy (S, "axpy", T.n, alpha, T.v1, T.s1, T.v2, T.s2);
    endfor
  endfor
endfor
end_file();

begin_file("copy");
for j = testcases
  for i = [s,d,c,z]
    S = context(i);
    T = test_vectors(S, j);
    test_copy (S, "copy", T.n, T.v1, T.s1, T.v2, T.s2);
  endfor
endfor
end_file();

begin_file("swap");
for j = testcases
  for i = [s,d,c,z]
    S = context(i);
    T = test_vectors(S, j);
    test_swap (S, "swap", T.n, T.v1, T.s1, T.v2, T.s2);
  endfor
endfor
end_file();

begin_file("scal");
for j = testcases
  for i = [s,d,c,z]
    S = context(i);
    T = test_vector(S, j);
    for alpha = [0, 0.1, 1.0]
      test_scal (S, "scal", T.n, alpha, T.v, T.s);
    endfor
    if (S.complex)
      for alpha = [0.1*I, 0.1+0.2*I, 1.0 + 0.3 * I]
        test_scal (S, "scal", T.n, alpha, T.v, T.s);
      endfor
    endif
  endfor
endfor
end_file();

begin_file("rotg");
for i = [s,d]
  S = context(i);
  for a = [-1.5, -1, -0.1, 0, 0.1, 1, 1.5]
    for b = [-1.5, -1, -0.1, 0, 0.1, 1, 1.5]
      test_rotg (S, "rotg", a, b);
    endfor
  endfor
endfor
end_file();

begin_file("rot");
for j = testcases
  for i = [s,d]
    S = context(i);
    T = test_vectors(S, j);
    for zz = [0, exp(I*pi/6), 0-I, -1]
      test_rot (S, "rot", T.n, T.v1, T.s1, T.v2, T.s2, \
                 real(zz), imag(zz));
    endfor
  endfor
endfor
end_file();

# for i = [s,d]
#   S = context(i);
#   for d1 = [0.1, 1, 1e8, 1e-32]
#     for d2 = [0.1, 1, 1e-8, 1e32]
#       for b1 = [-1.5,  0.1, 1, 1.5, 1e8, 1e-32]
#         for b2 = [-1.5, 0, 0.1, 1, 1.5, 3e9, 1e-27]
#           test_rotmg (S, "rotmg", d1, d2, b1, b2);
#         endfor
#       endfor
#     endfor
#   endfor
# endfor

begin_file("rotmg");
for j = testcases
  for i = [s,d]
    S = context(i);
    test_rotmg (S, "rotmg", rnd(), rnd(), rnd(), rnd());
  endfor
endfor
end_file();

begin_file("rotm");
for j = testcases
  for i = [s,d]
    S = context(i);
    T = test_vectors(S, j);
    v = [-2, -1, 0, 1];
    for k = 1:10
      h = [rnd(), rnd(); rnd(), rnd()];
      flag = v(rem(k,4)+1);
      test_rotm (S, "rotm", T.n, T.v1, T.s1, T.v2, T.s2, flag, h);
    endfor
  endfor
endfor
end_file();

## Level 2

begin_file("gemv");
for j = testcases2
  for i = [s,d,c,z]
    S = context(i);
    T = test_matvectors(S, j);
    for trans = Trans(S)
      for alpha = coeff(S)
        for beta = coeff(S)
          for order = [101, 102]
            test_gemv (S, "gemv", order, trans, T.m, T.n, alpha, T.A, \
                       T.lda, T.v1, T.s1, beta, T.v2, T.s2);
          endfor
        endfor
      endfor
    endfor
  endfor
endfor
end_file();

begin_file("gbmv");
for j = testcases2
  for i = [s,d,c,z]
    S = context(i);
    for trans = Trans(S)
      T = test_bmatvectors(S, j, trans);
      for alpha = coeff(S)
        for beta = coeff(S)
          for order = [101, 102]
            test_gbmv (S, "gbmv", order, trans, T.m, T.n, T.kl, T.ku, \
                       alpha, T.A, T.lda, T.v1, T.s1, beta, T.v2, T.s2);
          endfor
        endfor
      endfor
    endfor
  endfor
endfor
end_file();

begin_file("trmv");
for j = testcases2
  for i = [s,d,c,z]
    S = context(i);
    for trans = Trans(S)
      T = test_trmatvector(S, j);
      for order = [101, 102]
        for uplo = [121, 122]
          for diag = [131, 132]
            test_trmv (S, "trmv", order, uplo, trans, diag, T.n,
                       T.A, T.lda, T.v, T.s);
          endfor
        endfor
      endfor
    endfor
  endfor
endfor
end_file();


begin_file("tbmv");
for j = testcases2
  for i = [s,d,c,z]
    S = context(i);
    for trans = Trans(S)
      T = test_tbmatvector(S, j);
      for order = [101, 102]
        for uplo = [121, 122]
          for diag = [131, 132]
            test_tbmv (S, "tbmv", order, uplo, trans, diag, T.n, T.k,
                       T.A, T.lda, T.v, T.s);
          endfor
        endfor
      endfor
    endfor
  endfor
endfor
end_file();

begin_file("tpmv");
for j = testcases2
  for i = [s,d,c,z]
    S = context(i);
    for trans = Trans(S)
      T = test_tpmatvector(S, j);
      for order = [101, 102]
        for uplo = [121, 122]
          for diag = [131, 132]
            test_tpmv (S, "tpmv", order, uplo, trans, diag, T.n, 
                       T.A, T.v, T.s);
          endfor
        endfor
      endfor
    endfor
  endfor
endfor
end_file();

begin_file("symv");
for j = testcases2
  for i = [s,d]
    S = context(i);
    T = test_symatvectors(S, j);
    for alpha = coeff(S)
      for beta = coeff(S)
        for order = [101, 102]
          for uplo = [121, 122]
            for diag = [131, 132]
              test_hesymv (S, "symv", order, uplo, T.n, alpha, 
                           T.A, T.lda, 
                           T.v1, T.s1, beta, T.v2, T.s2);
            endfor
          endfor
        endfor
      endfor
    endfor
  endfor
endfor
end_file();

begin_file("hemv");
for j = testcases2
  for i = [c,z]
    S = context(i);
    T = test_symatvectors(S, j);
    for alpha = coeff(S)
      for beta = coeff(S)
        for order = [101, 102]
          for uplo = [121, 122]
            for diag = [131, 132]
              test_hesymv (S, "hemv", order, uplo, T.n, alpha, 
                           T.A, T.lda, T.v1, T.s1, beta, T.v2, T.s2);
            endfor
          endfor
        endfor
      endfor
    endfor
  endfor
endfor
end_file();

begin_file("hbmv");
for j = testcases2
  for i = [c,z]
    S = context(i);
    T = test_sbmatvectors(S, j);
    for alpha = coeff(S)
      for beta = coeff(S)
        for order = [101, 102]
          for uplo = [121, 122]
            for diag = [131, 132]
              test_hbsbmv (S, "hbmv", order, uplo, T.n, T.k, alpha, 
                           T.A, T.lda, T.v1, T.s1, beta, T.v2, T.s2);
            endfor
          endfor
        endfor
      endfor
    endfor
  endfor
endfor
end_file();

begin_file("sbmv");
for j = testcases2
  for i = [s,d]
    S = context(i);
    T = test_sbmatvectors(S, j);
    for alpha = coeff(S)
      for beta = coeff(S)
        for order = [101, 102]
          for uplo = [121, 122]
            for diag = [131, 132]
              test_hbsbmv (S, "sbmv", order, uplo, T.n, T.k, alpha, 
                           T.A, T.lda, T.v1, T.s1, beta, T.v2, T.s2);
            endfor
          endfor
        endfor
      endfor
    endfor
  endfor
endfor
end_file();

begin_file("hpmv");
for j = testcases2
  for i = [c,z]
    S = context(i);
    T = test_spmatvectors(S, j);
    for alpha = coeff(S)
      for beta = coeff(S)
        for order = [101, 102]
          for uplo = [121, 122]
            for diag = [131, 132]
              test_hpspmv (S, "hpmv", order, uplo, T.n, alpha, T.A, T.v1, \
                           T.s1, beta, T.v2, T.s2);
            endfor
          endfor
        endfor
      endfor
    endfor
  endfor
endfor
end_file();

begin_file("spmv");
for j = testcases2
  for i = [s,d]
    S = context(i);
    T = test_spmatvectors(S, j);
    for alpha = coeff(S)
      for beta = coeff(S)
        for order = [101, 102]
          for uplo = [121, 122]
            for diag = [131, 132]
              test_hpspmv (S, "spmv", order, uplo, T.n, alpha, T.A, T.v1, \
                           T.s1, beta, T.v2, T.s2);
            endfor
          endfor
        endfor
      endfor
    endfor
  endfor
endfor
end_file();

begin_file("trsv");
for j = testcases2
  for i = [s,d,c,z]
    S = context(i);
    for trans = Trans(S);
      T = test_trmatvector(S, j);
      for order = [101, 102]
        for uplo = [121, 122]
          for diag = [131, 132]
            test_trmv (S, "trsv", order, uplo, trans, diag, T.n,
                       T.A, T.lda, T.v, T.s);
          endfor
        endfor
      endfor
    endfor
  endfor
endfor
end_file();

begin_file("tbsv");
for j = testcases2
  for i = [s,d,c,z]
    S = context(i);
    for trans = Trans(S);
      T = test_tbmatvector(S, j);
      for order = [101, 102]
        for uplo = [121, 122]
          for diag = [131, 132]
            test_tbmv (S, "tbsv", order, uplo, trans, diag, T.n, T.k,
                       T.A, T.lda, T.v, T.s);
          endfor
        endfor
      endfor
    endfor
  endfor
endfor
end_file();

begin_file("tpsv");
for j = testcases2
  for i = [s,d,c,z]
    S = context(i);
    T = test_tpmatvector(S, j);
    for trans = Trans(S)
      for order = [101, 102]
        for uplo = [121, 122]
          for diag = [131, 132]
            test_tpmv (S, "tpsv", order, uplo, trans, diag, T.n, 
                       T.A, T.v, T.s);
          endfor
        endfor
      endfor
    endfor
  endfor
endfor
end_file();

begin_file("ger");
for j = testcases2
  for i = [s,d,c,z];
    S = context(i);
    T = test_matvectors(S, j);
    for alpha = coeff(S)
      for order = [101, 102]
        if (S.complex == 0)
          test_ger (S, "ger", order, T.m, T.n, alpha, T.v1, T.s1, \
                    T.v2, T.s2, T.A, T.lda);
        else
          test_ger (S, "geru", order, T.m, T.n, alpha, T.v1, T.s1, \
                    T.v2, T.s2, T.A, T.lda);
          test_ger (S, "gerc", order, T.m, T.n, alpha, T.v1, T.s1, \
                    T.v2, T.s2, T.A, T.lda);
        endif
      endfor
    endfor
  endfor
endfor
end_file();

begin_file("syr");
for j = testcases2
  for i = [s,d];
    S = context(i);
    T = test_trmatvector(S, j);
    for alpha = coeff(S)
      for order = [101, 102]
        for uplo = [121, 122]
          test_syr (S, "syr", order, uplo, T.n, alpha, T.v, T.s, \
                    T.A, T.lda);
        endfor
      endfor
    endfor
  endfor
endfor
end_file();

begin_file("her");
for j = testcases2
  for i = [c,z];
    S = context(i);
    T = test_trmatvector(S, j);
    R = S ; R.complex = 0;
    for alpha = coeff(R)
      for order = [101, 102]
        for uplo = [121, 122]
          test_her (S, "her", order, uplo, T.n, alpha, T.v, T.s, \
                      T.A, T.lda);
        endfor
      endfor
    endfor
  endfor
endfor
end_file();

begin_file("hpr");
for j = testcases2
  for i = [c,z];
    S = context(i);
    T = test_tpmatvector(S, j);
    R = S ; R.complex = 0;
    for alpha = coeff(R)
      for order = [101, 102]
        for uplo = [121, 122]
          test_hpr (S, "hpr", order, uplo, T.n, alpha, T.v, T.s, T.A);
        endfor
      endfor
    endfor
  endfor
endfor
end_file();


begin_file("spr");
for j = testcases2
  for i = [s,d];
    S = context(i);
    T = test_tpmatvector(S, j);
    for alpha = coeff(S)
      for order = [101, 102]
        for uplo = [121, 122]
          test_spr (S, "spr", order, uplo, T.n, alpha, T.v, T.s, T.A);
        endfor
      endfor
    endfor
  endfor
endfor
end_file();


begin_file("syr2");
for j = testcases2
  for i = [s,d];
    S = context(i);
    T = test_trmatvectors(S, j);
    for alpha = coeff(S)
      for order = [101, 102]
        for uplo = [121, 122]
          test_syr2 (S, "syr2", order, uplo, T.n, alpha, T.v1, T.s1, \
                     T.v2, T.s2, T.A, T.lda);
        endfor
      endfor
    endfor
  endfor
endfor
end_file();

begin_file("spr2");
for j = testcases2
  for i = [s,d];
    S = context(i);
    T = test_tpmatvectors(S, j);
    for alpha = coeff(S)
      for order = [101, 102]
        for uplo = [121, 122]
          test_spr2 (S, "spr2", order, uplo, T.n, alpha, T.v1, T.s1, \
                     T.v2, T.s2, T.A);
        endfor
      endfor
    endfor
  endfor
endfor
end_file();

begin_file("her2");
for j = testcases2
  for i = [c,z];
    S = context(i);
    T = test_trmatvectors(S, j);
    for alpha = coeff(S)
      for order = [101, 102]
        for uplo = [121, 122]
          test_syr2 (S, "her2", order, uplo, T.n, alpha, T.v1, T.s1, \
                     T.v2, T.s2, T.A, T.lda);
        endfor
      endfor
    endfor
  endfor
endfor
end_file();

begin_file("hpr2");
for j = testcases2
  for i = [c,z];
    S = context(i);
    T = test_trmatvectors(S, j);
    for alpha = coeff(S)
      for order = [101, 102]
        for uplo = [121, 122]
          test_spr2 (S, "hpr2", order, uplo, T.n, alpha, T.v1, T.s1, \
                     T.v2, T.s2, T.A, T.lda);
        endfor
      endfor
    endfor
  endfor
endfor
end_file();

## Level 3

begin_file("gemm");
for j = testcases3
  for i = [s,d,c,z]
    S = context(i);
    for transA = Trans(S)
      for transB = Trans(S)
        for alpha = coeff(S)
          for beta = coeff(S)
            for order = [101, 102]
              T = test_matmat(S, j, order, transA, transB);
              test_gemm (S, "gemm", order, transA, transB, T.m, T.n,
                         T.k, alpha, 
                         T.A, T.lda, T.B, T.ldb, beta, T.C, T.ldc);
            endfor
          endfor
        endfor
      endfor
    endfor
  endfor
endfor
end_file();

begin_file("symm")  ;
for j = testcases3
  for i = [s,d,c,z]
    S = context(i);
    for uplo = [121, 122]
      for side = [141, 142]
        for alpha = coeff(S)
          for beta = coeff(S)
            for order = [101, 102]
              T = test_symatmat(S, j, order, side);
              test_hesymm (S, "symm", order, side, uplo, T.m, T.n,
                         alpha, T.A, T.lda, T.B, T.ldb, beta, T.C, T.ldc);
            endfor
          endfor
        endfor
      endfor
    endfor
  endfor
endfor
end_file();

begin_file("hemm")  ;
for j = testcases3
  for i = [c,z]
    S = context(i);
    for uplo = [121, 122]
      for side = [141, 142]
        for alpha = coeff(S)
          for beta = coeff(S)
            for order = [101, 102]
              T = test_symatmat(S, j, order, side);
              test_hesymm (S, "hemm", order, side, uplo, T.m, T.n,
                         alpha, T.A, T.lda, T.B, T.ldb, beta, T.C, T.ldc);
            endfor
          endfor
        endfor
      endfor
    endfor
  endfor
endfor
end_file();

begin_file("syrk");
for j = testcases3
  for i = [s,d,c,z]
    S = context(i);
    for uplo = [121, 122]
      for trans = Trans(S)
        if (S.complex && trans == 113)
          continue; # ConjTrans not allowed for complex case, 
        endif
        for alpha = coeff(S)
          for beta = coeff(S)
            for order = [101, 102]
              T = test_syrkmatmat(S, j, order, trans);
              test_syrk (S, "syrk", order, uplo, trans, T.n, T.k,
                         alpha, T.A, T.lda, beta, T.C, T.ldc);
            endfor
          endfor
        endfor
      endfor
    endfor
  endfor
endfor
end_file();

begin_file("herk");
for j = testcases3
  for i = [c,z] #[s,d] #c,z]
    S = context(i);
    for uplo = [121, 122]
      for trans = Trans(S)
        if (S.complex && trans == 112)
          continue; # Trans not allowed for complex case, 
        endif
        T = S ; T.complex = 0;
        for alpha = coeff(T)
          for beta = coeff(T)
            for order = [101, 102]
              T = test_syrkmatmat(S, j, order, trans);
              test_herk (S, "herk", order, uplo, trans, T.n, T.k,
                         alpha, T.A, T.lda, beta, T.C, T.ldc);
            endfor
          endfor
        endfor
      endfor
    endfor
  endfor
endfor
end_file();


begin_file("syr2k");
for j = testcases3
  for i = [s,d,c,z]
    S = context(i);
    for trans = Trans(S)
      if (S.complex && trans == 113)
        continue; # ConjTrans not allowed for complex case, 
      endif
      for alpha = coeff(S)
        for beta = coeff(S)
          for order = [101, 102]
            for uplo = [121, 122]
              T = test_syr2kmatmat(S, j, order, trans);
              test_syr2k (S, "syr2k", order, uplo, trans, T.n, T.k, alpha, 
                          T.A, T.lda, T.B, T.ldb, beta, T.C, T.ldc);
            endfor
          endfor
        endfor
      endfor
    endfor
  endfor
endfor
end_file();

begin_file("her2k");
for j = testcases3
  for i = [c,z]
    S = context(i);
    for trans = Trans(S)
      if (S.complex && trans == 112)
        continue; # Trans not allowed for complex case, 
      endif
      for alpha = coeff(S)
        R = S; R.complex = 0;
        for beta = coeff(R)
          for order = [101, 102]
            for uplo = [121, 122]
              T = test_syr2kmatmat(S, j, order, trans);
              test_her2k (S, "her2k", order, uplo, trans, T.n, T.k, alpha, 
                          T.A, T.lda, T.B, T.ldb, beta, T.C, T.ldc);
            endfor
          endfor
        endfor
      endfor
    endfor
  endfor
endfor
end_file();


begin_file("trmm");
for j = testcases3
  for i = [s,d,c,z]
    S = context(i);
    for trans = Trans(S)
      for alpha = coeff(S)
        for side = [141, 142]
          for order = [101, 102]
            for uplo = [121, 122]
              for diag = [131, 132]
                T = test_trmatmat(S, j, order, side, trans);
                test_trmm (S, "trmm", order, side, uplo, trans, diag, T.m, T.n,
                           alpha, T.A, T.lda, T.B, T.ldb);
              endfor
            endfor
          endfor
        endfor
      endfor
    endfor
  endfor
endfor
end_file();

begin_file("trsm");
for j = testcases3
  for i = [s,d,c,z]
    S = context(i);
    for trans = Trans(S)
      for alpha = coeff(S)
        for side = [141, 142]
          for order = [101, 102]
            for uplo = [121, 122]
              for diag = [131, 132]
                T = test_trmatmat(S, j, order, side, trans);
                test_trmm (S, "trsm", order, side, uplo, trans, diag, T.m, T.n,
                           alpha, T.A, T.lda, T.B, T.ldb);
              endfor
            endfor
          endfor
        endfor
      endfor
    endfor
  endfor
endfor
end_file();
