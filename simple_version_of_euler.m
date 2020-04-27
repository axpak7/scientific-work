function [] = simple_version_of_euler(x_0, x_n, h)

%Numerical constants
n = int16((x_n-x_0) / h)+1;% numbers of steps
x=zeros(1,n);
y=zeros(1,n);
%y_prev = 1;
x(1) = x_0;
y(1) = 1;
for i = 2:n
    %x(i)=x_0;
    %y(i)=y_prev;
    %term = h * (y_prev + x_0)/(y_prev - x_0);
    %y_prev = y_prev + term;
    %x_0 = x_0 + h;
    
    term = h * dydx(y(i-1), x(i-1));
    %disp(term);
    y(i) = y(i-1) + term;
    x(i) = x(i-1) + h;
    %disp(y);
end
plot(x,y);
disp(vpa(x,9));
disp(vpa(y,9));
end