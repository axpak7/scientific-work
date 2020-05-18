function y_final = sophisticated_structural_method3(x_0, x_n, n)
%constants
h = (x_n-x_0)/n;
x = zeros(1,n);
y1 = zeros(1,n);
y2 = zeros(1,n);
y3 = zeros(1,n);
y4 = zeros(1,n);
%initial conditions
x(1) = x_0;
y1(1) = 1;
y2(1) = 1;
y3(1) = 1;
y4(1) = 1;
z = [0 0 0 0];
%algorithm
for i = 1:n
    x(i+1) = x(1) + i*h;
    k11 = h * dy1dx(x(i), 0, y2(i), 0, y4(i));
    k31 = h * dy3dx(x(i), 0, 0, 0, y4(i));
    k41 = h * dy4dx(x(i) + h/3, y1(i)+1/3*k11, 0, 0, 0);
    k21 = h * dy2dx(x(i) + h/3, 0, 0, y3(i)+1/3*k31, y4(i)+1/3*k41);
    k12 = h * dy1dx(x(i) + 2*h/3, 0, y2(i)+2/3*k21, 0, y4(i)+2/3*k41);
    k32 = h * dy3dx(x(i) + 2*h/3, 0, 0, 0, y4(i)+2/3*k41);
    k42 = h * dy4dx(x(i) + h, y1(i)+k12, 0, 0, 0);
    k22 = h * dy2dx(x(i) + h, 0, 0, y3(i)+k32, y4(i)+k41);
    term1 = 1/4*(k11+3*k12) + z(1);
    term3 = 1/4*(k31+3*k32) + z(3);
    term4 = 1/4*(3*k41+k42) + z(4);
    term2 = 1/4*(3*k21+k22) + z(2);
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