function y_final = sophisticated_runge_kutta3(x_0, x_n, n)
%constants
h = (x_n-x_0)/n;
x = zeros(1,n);
y1 = zeros(1,n);
y2 = zeros(1,n);
y3 = zeros(1,n);
y4 = zeros(1,n);
c2=1/3;
c3=2/3;
a21=1/3;
a31=0;
a32=2/3;
b1=1/4;
b2=0;
b3=3/4;
%initial conditions
x(1) = x_0;
y1(1) = 1;
y2(1) = 1;
y3(1) = 1;
y4(1) = 1;
z = [0 0 0 0];

for i = 1:n
    x(i+1) = x(1) + i*h;
    k11 = dy1dx(x(i),y1(i),y2(i),y3(i),y4(i));
    k12 = dy2dx(x(i),y1(i),y2(i),y3(i),y4(i));
    k13 = dy3dx(x(i),y1(i),y2(i),y3(i),y4(i));
    k14 = dy4dx(x(i),y1(i),y2(i),y3(i),y4(i));
    k21 = dy1dx(x(i)+c2*h, y1(i)+h*a21*k11, y2(i)+h*a21*k12, y3(i)+h*a21*k13, y4(i)+h*a21*k14);
    k22 = dy2dx(x(i)+c2*h, y1(i)+h*a21*k11, y2(i)+h*a21*k12, y3(i)+h*a21*k13, y4(i)+h*a21*k14);
    k23 = dy3dx(x(i)+c2*h, y1(i)+h*a21*k11, y2(i)+h*a21*k12, y3(i)+h*a21*k13, y4(i)+h*a21*k14);
    k24 = dy4dx(x(i)+c2*h, y1(i)+h*a21*k11, y2(i)+h*a21*k12, y3(i)+h*a21*k13, y4(i)+h*a21*k14);
    k31 = dy1dx(x(i)+c3*h, y1(i)+h*(a31*k11+a32*k21), y2(i)+h*(a31*k12+a32*k22), y3(i)+h*(a31*k13+a32*k23), y4(i)+h*(a31*k14+a32*k24));
    k32 = dy2dx(x(i)+c3*h, y1(i)+h*(a31*k11+a32*k21), y2(i)+h*(a31*k12+a32*k22), y3(i)+h*(a31*k13+a32*k23), y4(i)+h*(a31*k14+a32*k24));
    k33 = dy3dx(x(i)+c3*h, y1(i)+h*(a31*k11+a32*k21), y2(i)+h*(a31*k12+a32*k22), y3(i)+h*(a31*k13+a32*k23), y4(i)+h*(a31*k14+a32*k24));
    k34 = dy4dx(x(i)+c3*h, y1(i)+h*(a31*k11+a32*k21), y2(i)+h*(a31*k12+a32*k22), y3(i)+h*(a31*k13+a32*k23), y4(i)+h*(a31*k14+a32*k24));
    term1 = h*(b1*k11+b2*k21+b3*k31) + z(1);
    term2 = h*(b1*k12+b2*k22+b3*k32) + z(2);
    term3 = h*(b1*k13+b2*k23+b3*k33) + z(3);
    term4 = h*(b1*k14+b2*k24+b3*k34) + z(4);
    newy1 = y1(i) + term1;
    newy2 = y2(i) + term2;
    newy3 = y3(i) + term3;
    newy4 = y4(i) + term4;
    z(1) = term1 - (newy1 - y1(i));
    z(2) = term2 - (newy2 - y2(i));
    z(3) = term3 - (newy3 - y3(i));
    z(4) = term4 - (newy4 - y4(i));
    y1(i+1) = newy1;
    y2(i+1) = newy2;
    y3(i+1) = newy3;
    y4(i+1) = newy4;
end
y_final = [y1(n+1) y2(n+1) y3(n+1) y4(n+1)];
end