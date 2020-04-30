function y_final = sophisticated_version_of_euler(x_0, x_n, n)

%Numerical constants
%n = int16((x_n-x_0) / h)+1;% numbers of steps
h = single((x_n-x_0)/(n-1));
%n = 10;
x=zeros(1,n);
y=zeros(1,n);
x(1) = x_0;
y(1) = 1;
z = zeros(1,n);

for i = 1:n
    term = single(h * dydx(y(i), x(i)) + z);
    newy = single(y(i) + term);
    z = single(term - (newy - y(i)));
    y(i+1) = single(newy(i));
    x(i+1) = single(x(i) + h);
end

%plot(x,y);
%disp(vpa(x,9));
%disp(vpa(y,9));
y_final = single(y(n));
%disp(vpa(y_final,7));
end