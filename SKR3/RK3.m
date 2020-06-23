function y_final = RK3(f, x0, xend, y0, N)
%f - ������ ���������� �� ������� ������ ������
% f{i}(x,y) � i-� ��������� f
%x0 - ��������� �������� �������
%xend - �������� �������� ������
%y0 - ������ ��������� ����������
%N - ����� �����
%[n0 n1 n2] - ���������� ��������� � �����, ������ � ������ �����������
%������� ��������������

%constants
h = (xend-x0)/N;
neq = size(y0,1);

y1 = zeros(neq, 1);


% ���� ���, ��������, �� ��� �����:
% ��������� K � X, Y
for j = 1:N
 % ������ ����
    K1 = zeros(neq,1);
    X = x0;
    Y = y0;
    for i = 1:neq
        K1(i) = h * f(X,Y,i);
    end

    K2 = zeros(neq,1);
    X = x0 + 1/2 * h;
    Y = y0 + 1/2 * K1;
    for i = 1:neq
        K2(i) = h * f(X,Y,i);
    end

    K3 = zeros(neq,1);
    X = x0 + h;
    Y = y0 - K1 + 2 * K2;
    for i = 1:neq
        K3(i) = h * f(X,Y,i);
    end
    
    y1 = y0 + 1/6 * ( K1 + 4 * K2 + K3);
    
    x0 = x0 + h;
    y0 = y1;
end

y_final = y1;

end

