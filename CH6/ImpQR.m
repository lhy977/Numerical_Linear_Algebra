%07222113邓俊烽
%算法6.4.3隐式QR
function [H,Q,lambda] = ImpQR(A)
% 用于求解矩阵 A 的所有特征值
% 输入: A - 输入的方阵
% 输出: H - 上 Hessenberg 矩阵
%       Q - 正交变换矩阵
%       lambda - 特征值向量
[n, ~] = size(A);
if ~min(min(tril(A, -1) == 0))    % 判断A是否为Hessenberg阵
    [Q, H] = HessenbergDec(A);
else
    Q = eye(n);
    H = A;
end
tic
k = 0;
while 1
    k = k + 1;
    % 收敛性判定
    for i = 2 : n
        if abs(H(i, i - 1)) < (abs(H(i, i)) + abs(H(i - 1, i - 1)))*1e-15
            H(i, i - 1) = 0;
        end
    end
    % 寻找不可约Hessenberg阵: 次对角线均不为0
    d = diag(H(2 : end, 1 : end - 1));
    breakpoints = find(d == 0);
    breakpoints = [0; breakpoints; n];
    for j = length(breakpoints) : -1 : 2
            l = breakpoints(j - 1);
            m = breakpoints(j);
        if breakpoints(j) - breakpoints(j - 1) < 3
            % 在右下角分出一个或一对特征值
        else
            break;
        end
    end

    if m - l <= 2
        % 已经是Schur阵
        break;
    else
        [P, H(l + 1 : m, l + 1 : m)] = DoubleStepQR(H(l + 1 : m, l + 1 : m));
        Q(:, l + 1 : m) = Q(:, l + 1 : m)*P;
        H(1 : l, l + 1 : m) = H(1 : l, l + 1 : m)*P;
        H(l + 1 : m, m + 1 : end) = P'*H(l + 1 : m, m + 1 : end);
    end
end
T = H;

lambda = diag(T);
d = diag(T(2 : end, 1 : end - 1));
breakpoints = find(d == 0);
breakpoints = [0; breakpoints; n];
for i = 1 : length(breakpoints) - 1
    if breakpoints(i + 1) > breakpoints(i) + 1
        r = breakpoints(i);
        [lambda(r + 1), lambda(r + 2)] = eigDim2(T(r + 1 : r + 2, r + 1 : r + 2));
    end
end
lambda = sort(lambda);
disp(['迭代次数为: ', num2str(k)]);
disp(['迭代时间为',num2str(toc),'s']);
end