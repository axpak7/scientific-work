function y_final = sophisticated_runge_kutta3(x_0, x_n, n)
%constants
h = (x_n-x_0)/n;
x = zeros(1,n);
y = zeros(1,n);
c2=1/3;
c3=2/3;
a21=1/3;
a31=0;
a32=2/3;
b1=1/4;
b2=0;
b3=3/4;
x(1) = x_0;
y(1) = 1;
z = 0;

for i = 1:n
    x(i+1) = x(1) + i*h;
    k1 = dydx(x(i),y(i));
    k2 = dydx(x(i)+c2*h, y(i)+h*a21*k1);
    k3 = dydx(x(i)+c3*h, y(i)+h*(a31*k1+a32*k2));
    term = h*(b1*k1+b2*k2+b3*k3) + z;
    newy = y(i) + term;
    z = term - (newy - y(i));
    y(i+1) = newy;
end
y_final = y(n+1);
end