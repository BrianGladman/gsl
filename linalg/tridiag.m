d = [4, 2, 1, 2, 4];
f = [1, 1, 1, 1, 1];


#b = [1, 1, 1, 1, 1];
b = [30, -24, 3, 21, -30];

x     = [0,0,0,0,0];
alpha = [0,0,0,0,0];
gamma = [0,0,0,0,0];
delta = [0,0,0,0,0];
z     = [0,0,0,0,0];

n = length(d);

A = diag(d) + diag(f(1:n-1),1) + diag(f(1:n-1),-1) + diag(f(n),n-1) + diag(f(n),-(n-1))

alpha(1) = d(1)
gamma(1) = f(1)/alpha(1)
delta(1) = f(n)/alpha(1)
for i = 2:n-2 
  alpha(i) = d(i) - f(i-1)*gamma(i-1)
  gamma(i) = f(i) / alpha(i)
  delta(i) = -delta(i-1)*f(i-1)/alpha(i)
endfor

s = sum(alpha(1:n-2).*(delta(1:n-2).^2))

alpha(n-1) = d(n-1) - f(n-2)*gamma(n-2)
gamma(n-1) = (f(n-1) - f(n-2)*delta(n-2))/alpha(n-1)
alpha(n) = d(n) - s - alpha(n-1) * gamma(n-1) * gamma(n-1)

z(1) = b(1)
for i = 2:n-1
  z(i) = b(i) - z(i-1)*gamma(i-1)
endfor
s = sum(delta(1:n-2).*z(1:n-2))
z(n) = b(n) - s - gamma(n-1)*z(n-1)
c = z ./ alpha

x(n) = c(n)
x(n-1) = c(n-1) - gamma(n-1)*x(n)
for i = n-2:-1:1
  x(i) = c(i) - gamma(i)*x(i+1) - delta(i)*x(n)
endfor


