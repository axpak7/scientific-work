%h = [1e-23 1e-22 1e-21 1e-20 1e-19 1e-18 1e-17 1e-16 1e-15 1e-14 1e-13 1e-12 1e-11 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1];
n = [8388608 4194304 2097152 1048576 524288 262144 131072 65536 32768 16384 8192 4096 2048 1024 512 256 128 64 32 16 8 4 2];
h = (x_end - x_start) * ones(size(n))./n;
E1 = zeros(size(n));
E2 = zeros(size(n));
x_start = 0;
x_end = 0.5;


for i = 1:23
    E1(i) = abs(runge_kutta3(x_start,x_end,n(i)) - solution(x_end));
    E2(i) = abs(sophisticated_runge_kutta3(x_start,x_end,n(i)) - solution(x_end));
end
E2(E2 == 0) = min(E2(E2 ~= 0));


loglog(h,E1,'-s');
hold on;
loglog(h,E2,'-o');
title('runge kutta 3');
xlabel('h');
ylabel('|E|');
hold off;