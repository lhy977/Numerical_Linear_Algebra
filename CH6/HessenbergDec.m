% 07222113邓俊烽
%算法6.4.1.Hessenberg分解
function [Q, H] = HessenbergDec(H)

[n, ~] = size(H);
Q = eye(n);
for k = 1 : n - 2
    % 求第k列k + 1 : n部分的Householder变换Hk
    Hk = HouseTrans(H(k + 1 : n, k));
    % 左边作用在A的第k + 1 : n行上, 消第k列元素
    H(k + 1 : n, k : n) = Hk*H(k + 1 : n, k : n);
    % 将A部分赋成严格的0
    H(k + 2 : n, k) = zeros;
    H(1 : n, k + 1 : n) = H(1 : n, k + 1 : n)*Hk;
    % 记录Hk
    Q = Q*[eye(k), zeros(k, n - k); zeros(n - k, k), Hk];
end