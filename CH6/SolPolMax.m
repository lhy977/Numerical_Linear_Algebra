%07222113邓俊烽
%求解多项式的最大根
function x = SolPolMax(a)
n = length(a);
A = diag(ones(1,n-1),-1);% 构造下对角线为 1 的矩阵
A(:,end) = -a'; % 将最后一列设为 -a
x = PowerMethod(A);
end
