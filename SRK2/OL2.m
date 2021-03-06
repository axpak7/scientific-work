function y_final = OL2(f, x0, xend, y0, N, n0, n1, n2)
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
neq = n0 + n1 + n2;

%1 - ������ ������ � ������� ������
% n0 � ��������� ������ �� ������� ������
n11 = n0 + 1; % ������ ������ � ������ ������
n12 = n0 + n1; % ��������� ������ � ������ ������
n21 = n12 + 1; % ������ ������ �� ������ ������
n22 = n12 + n2; % n22 = neq � ����� ����� ���������, ��������� ������ �� ������ ������
y1 = zeros(neq, 1);


% ���� ���, ��������, �� ��� �����:
% ��������� K � X, Y
for j = 1:N
    % ������ ����
    K1 = zeros(neq,1);
    % ������ 0
    X = x0; % ��������� �����
    Y = y0; % ��������� �����
    for i = 1:n0
        K1(i) = h * f(X,Y,i);
    end

    % ������ 1
    X = x0; % ��������� �����
    Y = y0; % ��������� �����
    for i = n11:n12
        K1(i) = h * f(X,Y,i);
    end

    % ������ 2
    X = x0 + 1/2 * h;
    Y = y0 + 1/2 * K1;
    for i = n21:n22
        K1(i) = h * f(X,Y,i);
        Y(i) = y0(i) + 1/2 * K1(i); % ��� ����� ��� ���� ������ ��� ��������� ����������� ������ ������
    end

    % ������ ����
    K2 = zeros(neq,1);
    % ������ 0
    X = x0 + 1/2 * h;
    Y = y0 + 1/2 * K1;
    for i = 1:n0
        K2(i) = h * f(X,Y,i);
    end

    % ������ 1
    X = x0 + 1/2 * h;
    Y(1:n0) = y0(1:n0) + 1/2 * K2(1:n0); % �� ������ 0
    Y(n21:n22) = y0(n21:n22) + 1/2 * K1(n21:n22); % �� ������ 2

    for i = n11:n12
        K2(i) = h * f(X,Y,i);
        Y(i) = y0(i) + 1/2 * K2(i); % ��� ����� ��� ���� ������ ��� ��������� ����������� ������ ������
    end

    % ����:
    y1(1:n12) = y0(1:n12) + K2(1:n12);
    y1(n21:n22) = y0(n21:n22) + K1(n21:n22);
    
    x0 = x0 + h;
    y0 = y1; % ��� ���������� ����
end

y_final = y1;

end

