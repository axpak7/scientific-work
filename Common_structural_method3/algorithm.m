%h = [1e-23 1e-22 1e-21 1e-20 1e-19 1e-18 1e-17 1e-16 1e-15 1e-14 1e-13 1e-12 1e-11 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1];
%n = [16777216 8388608 4194304 2097152 1048576 524288 262144 131072 65536 32768 16384 8192 4096 2048 1024 512 256 128 64 32 16 8 4 2];
n = [131072 65536 32768 16384 8192 4096 2048 1024 512 256 128 64 32 16 8 4 2];
x_start = 0;
x_end = 1;
h = (x_end - x_start) * ones(size(n))./n;
E1 = zeros(size(n));
E2 = zeros(size(n));
f = {@f1, @f2, @f3, @f4};
y0 = [1; 1; 1; 1];

for i = 1:17
    E1(i) = norm(structural_method3(f, x_start, x_end, y0, n(i), 0, 2, 2)' - solution(x_end),Inf);
    E2(i) = norm(sophisticated_structural_method3(f, x_start, x_end, y0, n(i), 0, 2, 2)' - solution(x_end),Inf);
    disp(i);
end
%E2(E2 == 0) = min(E2(E2 ~= 0));

loglog(h,E1,'-s');
hold on;
loglog(h,E2,'-o');
title('3-order structural method for system');
xlabel('h');
ylabel('|E|');
hold off;

function y = f1(x,y)
   y=2*x*y(4)^(1/5)*y(3);
end
function y = f2(x, y)
   y=2*x*y(3);  
end
function y = f3(x, y)
   y=-2*x*log(y(1));  
end
function y = f4(x,y)
   y=10*x*exp(5*(y(2)-1))*y(3);
end
function y = solution(x)
    y = [exp(sin(x^2)) sin(x^2)+1 cos(x^2) exp(5*sin(x^2))];
end