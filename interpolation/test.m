
function s = cspline (x, y)
  n = length(x);
  
  a = y;
  h = diff(x);

  A = zeros(n-2,n-2);
  A = A + 2*(diag(h(1:n-2)) + diag(h(2:n-1)));
  A = A + diag(h(2:n-2),1);
  A = A + diag(h(2:n-2),-1);

  da = diff(a);
  
  for i = 2:(n-1)
    g(i) = 3 * (da(i+1)/ h(i+1) -  di(i)/h(i));
  endfor

  c = A \ g;
  
  b = diff(a)./h - (1/3) * h ./ c;
  d = (1/3)*diff(c)./h ;


endfunction
  
  
   