%07222113邓俊烽
%Page169反幂法
function [Eig, Value] = InversePowerMethod(A)
%求解最小的特征值
n = size(A, 1);
maxtime = 100*n;
x0=zeros(n,1);
eps=1e-12;
x = x0;
x(1) = x(1) + 1; 
t = 1;
while norm(x0 - x) > eps
    if t > maxtime
        disp('达到最大迭代次数.');
        break;
    else
        x0 = x;
        y = A\x;
        index = find(abs(y) == max(abs(y)));
        mu = y(index(1));
        x = y/mu; % 归一化处理
        t = t + 1;
    end
end
Value = x; % 最终得到的近似特征向量
Eig = 1/mu;
end