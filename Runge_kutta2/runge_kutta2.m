function y_final = runge_kutta2(x_0, x_n, n)
%constants
h = (x_n-x_0)/n;
x = zeros(1,n);
y = zeros(1,n);
c2=1/2;
a21=1/2;
b2=1;
%initial conditions
x(1) = x_0;
y(1) = 1;
%algorithm
for i = 1:n
    x(i+1) = x(1) + i*h;
    k1 = dydx(x(i),y(i));
    k2 = dydx(x(i)+c2*h, y(i)+h*a21*k1);
    y(i+1) = y(i) + b2*h*k2;
end
y_final = y(n+1);
end

