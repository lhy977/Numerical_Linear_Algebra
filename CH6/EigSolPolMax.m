%07222113邓俊烽
%调用系统eig求解多项式的最大根
function x = EigSolPolMax(a)
n = length(a);
A = diag(ones(1,n-1),-1);% 构造下对角线为 1 的矩阵
A(:,end) = -a'; % 将最后一列设为 -a
% 计算矩阵 A 的所有特征值
x = eig(A);
% 找到绝对值最大的特征值
[~, idx] = max(abs(x)); % 找到绝对值最大的索引
x = x(idx); % 获取对应的原值
end
