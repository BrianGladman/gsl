function j = svdjacobi (A)
  [m,n] = size(A);
  A
  sweep = 0;
  tol = 10 * eps;
  count = n*(n-1)/2
  for j=1:(n-1)
    for k = (j+1):n
      sweep++;
      disp("=========================================");
      sweep
      count
      j
      k
      x = A(:,j);
      y = A(:,k);
      Z = sum(abs(x.*y))
      P = 2*dot(x,y)
      Q = dot(x,x) - dot(y,y)
      v = sqrt(P*P+Q*Q)

      if (norm(y) <= tol * norm(x))
        disp("norm y is small")
        count--;
        continue
      elseif (abs(P) <= tol * norm(x) * norm(y))
        disp("P is small")
        count--;
        continue
      endif

      if (v == 0)
        sine = 1
        cosine = 0
      elseif Q < 0
        sine = 1
        cosine = 0
      else
        cosine = sqrt((v+Q)/(2*v))
        sine = P/(2*v*cosine)
      endif

      disp ("rotating...")

      X = cosine * x + sine * y;
      Y = -sine * x + cosine * y;

      A(:,j) = X;
      A(:,k) = Y;
      A
      count
    endfor
  endfor
  j = A;
endfunction


