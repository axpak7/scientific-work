function [] = sophisticated_version_of_euler(x_0, x_n, h)

%Numerical constants
n = int16((x_n-x_0) / h)+1;% numbers of steps
x=zeros(1,n);
y=zeros(1,n);
%y_prev = 1;

x(1) = x_0;
y(1) = 1;
z = zeros(size(y));
for i = 2:n
    term = h * dydx(y(i-1), x(i-1)) + z;
    %disp(term);
    newy = y(i-1) + term;
    %disp(y);
    z = term - (newy - y(i-1));
    y(i) = newy(i);
    x(i) = x(i-1) + h;
end
plot(x,y);
disp(vpa(x,9));
disp(vpa(y,9));
end