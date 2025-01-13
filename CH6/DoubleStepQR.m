%07222113邓俊烽
%算法6.4.2双重步位移的QR
function [P, Hess] = DoubleStepQR(A)
%用于将Hessberg矩阵进一步迭代为Hessberg矩阵
[n, ~] = size(A);
[~, Hess] = HessenbergDec(A);

s = Hess(n - 1, n - 1) + Hess(n, n);
t = Hess(n - 1, n - 1)*Hess(n, n) - Hess(n - 1, n)*Hess(n, n - 1);

% M的第一列只有这三个元素
x = Hess(1, 1)*Hess(1, 1) + Hess(1, 2)*Hess(2, 1) - s*Hess(1, 1) + t;
y = Hess(2, 1)*(Hess(1, 1) + Hess(2, 2) - s);
z = Hess(2, 1)*Hess(3, 2);

% 确定Q第一列
P = HouseTrans([x, y, z]');
P = [P, zeros(3, n - 3); zeros(n - 3, 3), eye(n - 3)];
B = P*Hess*P;

for k = 1 : n - 3
    PTEMP = HouseTrans(B(k + 1 : k + 3, k));
    P(:, k + 1 : k + 3) = P(:, k + 1 : k + 3)*PTEMP;
    % P0 = [eye(k), zeros(k, 3), zeros(k, n - k - 3);...
    %     zeros(3, k), P, zeros(3, n - k - 3);...
    %     zeros(n - k - 3, k), zeros(n - k - 3, 3), eye(n - k - 3)]*P0;
    B(k + 1 : k + 3, :) = PTEMP*B(k + 1 : k + 3, :);
    B(:, k + 1 : k + 3) = B(:, k + 1 : k + 3)*PTEMP;
end

if n == 3
    k = 1;
else
    k = k + 1;
end
% 最后一步只需要消一个二维向量
PTEMP = HouseTrans(B(k + 1 : k + 2, k));
P(:, end - 1 : end) = P(:, end - 1 : end)*PTEMP;
% P0(end - 1 : end, :) = [eye(k), zeros(k, 2); zeros(2, k), P]*P0(end - 1 : end, :);
B(k + 1 : k + 2, :) = PTEMP*B(k + 1 : k + 2, :);
B(:, k + 1 : k + 2) = B(:, k + 1 : k + 2)*PTEMP;

Hess = B;
for i = 1 : n
    Hess(i + 2 : end, i) =  zeros(n - i - 1, 1);
end
end