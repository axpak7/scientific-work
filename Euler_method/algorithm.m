h = [1e-23 1e-22 1e-21 1e-20 1e-19 1e-18 1e-17 1e-16 1e-15 1e-14 1e-13 1e-12 1e-11 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1];
n = [8388608 4194304 2097152 1048576 524288 262144 131072 65536 32768 16384 8192 4096 2048 1024 512 256 128 64 32 16 8 4 2];
E1 = zeros(size(n));
E2 = zeros(size(n));
x_start = 0;
x_end = 0.5;
%n = [10000000 1000000 100000 10000 1000 100 10 1];
%n = 2.^[1:23];

for i = 1:23
    E1(i) = single(abs(simple_version_of_euler(x_start,x_end,n(i)) - solution(x_end)));
    E2(i) = single(abs(sophisticated_version_of_euler(x_start,x_end,n(i)) - solution(x_end)));
end

loglog(h,E1,'-s');
hold on;
loglog(h,E2,'-o');
title('euler method');
xlabel('h');
ylabel('|E|');
hold off;