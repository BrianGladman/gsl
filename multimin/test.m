function y = wood (x)
  t1 = x(1)^2 - x(2);
  t2 = x(3)^2 - x(4);
  
  y =  100 * t1 * t1 + (1 - x(1))^2 + 90 * t2 * t2 + (1 - x(3))^2 \
      + 10.1 * ( (1 - x(2))^2 + (1 - x(4))^2) \
      + 19.8 * (1 - x(2)) * (1 - x(4));
endfunction

function dy = dwood(x)
  t1 = x(1)^2 - x(2);
  t2 = x(3)^2 - x(4);
  dy = ones(1,4);
  dy(1) = 400 * x(1) * t1 - 2 * (1 - x(1));
  dy(2) = -200 * t1 - 20.2 * (1 - x(2)) - 19.8 * (1 - x(4));
  dy(3) = 360 * x(3) * t2 - 2 * (1 - x(3));
  dy(4) = -180 * t2 - 20.2 * (1 - x(4)) - 19.8 * (1 - x(2));
endfunction

format long 

x0 = [-3,-1,-3,-1];
f0 = wood(x0)
g0 = dwood(x0)
g0norm = norm(g0)
p0 = dwood(x0);

global x; global p;
x = x0; p = p0;

function f = fmin(lambda); 
  global x; global p; 
  g=dwood(x-lambda*p); 
  f=dot(g,p)/(norm(g)*norm(p)); 
endfunction

function beta = beta(g1, g0)
  beta = - dot(g1,g1)/dot(g0,g0);
endfunction

k0 = fsolve("fmin", 0.0025)

x = x0 - k0 * p0;

x1 = x
f1 = wood(x1)
g1 = dwood(x1)
g1norm = norm(g1)

beta0 = beta(g1, g0)

p = g1 - beta0 * p0
p1 = p

k1 = fsolve("fmin", 0.0041)

x = x1 - k1 * p1;

x2 = x
f2 = wood(x2)
g2 = dwood(x2)
g2norm = norm(g2)

beta1 = beta(g2, g1)

p = g2 - beta1 * p1
p2 = p
