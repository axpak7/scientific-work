h = [1e-5 1e-4 1e-3 1e-2 1e-1];
n = [1000000 100000 10000 1000 100 10];
x_start = 0;
x_end = 0.5;
diff_x = x_end - x_start;
E1 = zeros(1,5);
E2 = zeros(1,5);
for i = 1:5
    E1(i) = single(abs(simple_version_of_euler(x_start,x_end,n(i)) - solution(x_end)));
    E2(i) = single(abs(sophisticated_version_of_euler(x_start,x_end,n(i)) - solution(x_end)));
    %disp('1');
    %disp(E1);
    %disp('2');
    %disp(E2);
end
loglog(h,E1);
hold on;
loglog(h,E2);
xlabel('h');
ylabel('|E|');
hold off;