function y_final = simple_version_of_euler(x_0, x_n, n)

%Numerical constants
h = single((x_n-x_0)/n);
x=zeros(1,n);
y=zeros(1,n);
x(1) = x_0;
y(1) = 1;

for i = 1:n
    term = single(h * dydx(x(i), y(i)));
    y(i+1) = single(y(i) + term);
    x(i+1) = single(x(1) + i*h);
end

y_final = single(y(n+1));
end