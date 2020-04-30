function y_final = simple_version_of_euler(x_0, x_n, n)

%Numerical constants
%n = int16((x_n-x_0) / h)+1;% numbers of steps
%n = 10;
h = single((x_n-x_0)/(n-1));
x=zeros(1,n);
y=zeros(1,n);
x(1) = x_0;
y(1) = 1;

for i = 1:(n)
    term = single(h * dydx(y(i), x(i)));
    y(i+1) = single(y(i) + term);
    x(i+1) = single(x(i) + h);
end

%plot(x,y);
%disp(vpa(x,9));
%disp(vpa(y,9));
y_final = single(y(n));
%disp(vpa(y_final,7));
end