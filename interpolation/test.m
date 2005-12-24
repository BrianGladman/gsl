
function [xx,yy,dy,iy] = cspline (x, y, factor)
  n = length(x);
  
  a = y;
  h = diff(x);

##### natural cspline
#   A = zeros(n-2,n-2);
#   A = A + (diag(h(1:n-2)) + diag(h(2:n-1)))/3.0;
#   A = A + diag(h(2:n-2),1)/6;
#   A = A + diag(h(2:n-2),-1)/6;

#   da = diff(a);
  
#   g = zeros(n-2,1);
#   for i = 1:(n-2)
#     g(i) = da(i+1)/h(i+1) - da(i)/h(i);
#   endfor

#   M = [0; A \ g; 0];

##### periodic cspline
#    A = zeros(n-1,n-1);
#    A = A + (diag(h(1:n-1)) + diag(shift(h(1:n-1),1)))/3.0;
#    A = A + diag(h(1:n-2),1)/6;
#    A = A + diag(h(1:n-2),-1)/6;
#    A(1,n-1) = h(n-1)/6
#    A(n-1,1) = h(n-1)/6

#    da = diff(a);
  
#    g = zeros(n-1,1);
#    g(1) = da(1)/h(1) - da(n-1)/h(n-1)
#    for i = 2:(n-1)
#      g(i) = da(i)/h(i) - da(i-1)/h(i-1);
#    endfor

#    M = A \ g;
#    M(n) = M(1);

# periodic 3pt cubic
  
  
  da = diff(a);
  M = 6 * [1;-1] * da(1) / (h(1)*h(2))
  M(3) = M(1);

  area=zeros(n-1,1);
  for k = 1:n-1
    z=x(k+1)
    area(k) = M(k+1)*(((z-x(k))^2*((z-x(k))^2-2*h(k)^2))/(24*h(k))) - M(k)*((x(k+1)-z)^2-h(k)^2)^2/(24*h(k))+ y(k+1)*(z-x(k))^2/(2*h(k)) + y(k)*(h(k)^2-(x(k+1)-z)^2)/(2*h(k));
  end

  integral=[0;cumsum(area)];

  xx = zeros(factor*(n-1),1);
  yy = zeros(factor*(n-1),1);
  dy = zeros(factor*(n-1),1);
  iy = zeros(factor*(n-1),1);

  for k = 1:n-1
    for j = 0:factor
      i = (k-1)*factor + j + 1;
      z = x(k) + h(k)*j/factor;
      xx(i) = z;
      yy(i) = M(k+1)*(((z-x(k))*((z-x(k))^2-h(k)^2))/(6*h(k))) + M(k)*(((x(k+1)-z)*((x(k+1)-z)^2-h(k)^2))/(6*h(k)))+ y(k+1)*(z-x(k))/h(k) + y(k)*(x(k+1)-z)/h(k);
      dy(i) = M(k+1)*(3*(z-x(k))^2-h(k)^2)/(6*h(k)) - M(k)*(3*(x(k+1)-z)^2-h(k)^2)/(6*h(k)) +(y(k+1)-y(k)) / h(k);
      iy(i) = integral(k) + M(k+1)*(((z-x(k))^2*((z-x(k))^2-2*h(k)^2))/(24*h(k))) - M(k)*((x(k+1)-z)^2-h(k)^2)^2/(24*h(k))+ y(k+1)*(z-x(k))^2/(2*h(k)) + y(k)*(h(k)^2-(x(k+1)-z)^2)/(2*h(k));
    end
  end

  #b = diff(a)./h - (1/3) * h ./ c;
  #d = (1/3)*diff(c)./h ;
endfunction
  
  
   
