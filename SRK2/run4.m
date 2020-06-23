function run4()
format long
K = 20
n = 2.^(1:K)
E1 = zeros(size(n));
E2 = zeros(size(n));
E3 = zeros(size(n));

x_start = 0;
x_end = 1;
h = (x_end - x_start) * ones(size(n))./n;

N = 4;
neq = 2 * N + 2;
n0 = 2;
n11 = n0 + 1;
n12 = n0 + N;
n21 = n12 + 1;
n22 = n12 + N;

y0 = initial;

for k = 1:K
    solOL = OL2(@odefun, x_start, x_end, y0, n(k), 2, N, N);
    solRK = RK2(@odefun, x_start, x_end, y0, n(k));
    solSOL = SOL2(@odefun, x_start, x_end, y0, n(k), 2, N, N);
    E1(k) = norm(solOL - solution(x_end),Inf);
    E2(k) = norm(solRK - solution(x_end),Inf);
    E3(k) = norm(solSOL - solution(x_end),Inf);
    disp(k);
end

%E2(E2 == 0) = min(E2(E2 ~= 0));
figure
    loglog(h,E1,'-s');
    hold on;
    loglog(h,E2,'-o');
    title('Structural method2 and RK2');
    xlabel('h');
    ylabel('|E|');
    hold off;
figure
    loglog(h,E1,'-s');
    hold on;
    loglog(h,E3,'-o');
    title('Structural method2 and Sophisticated Structural method2');
    xlabel('h');
    ylabel('|E|');
    hold off;

function y0 = initial()
   y0 = ones(neq, 1);
end

function y = solution(t)
   y = zeros(neq,1);
   tt = t;
   for i = 1:neq
       y(i) = exp(tt);
       tt = -tt;
   end
end

function f = odefun(t, y, i)
    if i ==1
        term = 0;
        for j=n11:n22
            term = term + (-1)^(rem(ceil(j/2),2)) * y(j)^2;
        end
        f = y(1)^2*y(2) + term;
    elseif i == 2
        term = 0;
        for j=n11:n22
            term = term + (-1)^(rem(ceil((j+1)/2),2)) * y(j)^2;
        end
        f = -y(2)^2*y(1) + term;
    elseif i < n21
        term = 0;
        for j=n21:n22
            term = term + (-1)^(rem(ceil((j+i)/2),2)) * y(j)^2;
        end
        f = (-1)^(rem(i+1,2)) * y(i-2)^2*y(i-1) + term;        
    else
        term = 0;
        for j=n11:n12
            term = term + (-1)^(rem(ceil((j+i)/2),2)) * y(j)^2;
        end
        f = (-1)^(rem(i+1,2)) * y(i-2)^2*y(i-1) + term;        
    end
end

end