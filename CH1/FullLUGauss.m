% 07222113邓俊烽
% 算法1.2.1全主元Gauss消去
function [L,U,P,Q] = FullLUGauss(A)
%L是一个n*n下三角矩阵,U是一个n*n的上三角矩阵,P和Q是置换矩阵,A是一个n*n的矩阵
%用于解决PAQ=LU;
n = length(A);
L = zeros(n);
U = A;
p = zeros(n,1);
q = zeros(n,1);
for k = 1:n-1
	%找索引
	MaxValue = max(max(abs(U(k:n,k:n))));%首先寻找最大元
	[pk,qk] = find(MaxValue == abs(U(k:n,k:n)));%然后寻找最大元的位置
    %可能有多个最大值只取第一个
    pk = pk(1);
	qk = qk(1);
    %将子矩阵中的索引映射回原始矩阵的索引。
	pk = pk + k - 1;
	qk = qk + k - 1;
    %记录交换的行和列
    p(k) = pk;
	q(k) = qk;
	%交换行和列
	[U(k,:),U(pk,:)] = deal(U(pk,:),U(k,:));
	[U(:,k),U(:,qk)] = deal(U(:,qk),U(:,k));
	[L(k,:),L(pk,:)] = deal(L(pk,:),L(k,:));
	[L(:,k),L(:,qk)] = deal(L(:,qk),L(:,k));
    % 更新 L 和 U
    L(k+1:n,k) = U(k+1:n,k)/U(k,k);
    U(k+1:n,k) = 0;
    L(k,k) = 1;
    U(k+1:n,k+1:n) = U(k+1:n,k+1:n) - L(k+1:n,k)*U(k,k+1:n);
end
L(n,n) = 1;
%构造置换矩阵P,Q
P = eye(n);
Q = eye(n);
for k = 1:n-1
	[P(k,:),P(p(k),:)] = deal(P(p(k),:),P(k,:));
	[Q(:,k),Q(:,q(k))] = deal(Q(:,q(k)),Q(:,k));
end
end