function y_final = sophisticated_version_of_euler(x_0, x_n, n)
%constants
h = single((x_n-x_0)/(n));
x=zeros(1,n);
y=zeros(1,n);
%initial conditions
x(1) = x_0;
y(1) = 1;
z = single(0);
%algorithm
for i = 1:n
    term = single(h * dydx(x(i), y(i)) + z);
    newy = single(y(i) + term);
    z = single(term - (newy - y(i)));
    y(i+1) = single(newy);
    x(i+1) = single(x(1) + i*h);
end

plot(x,y);
y_final = single(y(n+1));
end