% 07222113邓俊烽
% 算法1.2.2列主元Gauss消去
function [L,U,P] = ColLUGauss(A)
%L是一个n*n上三角矩阵,U是一个n*n的下三角矩阵,P是置换矩阵,A是一个n*n的矩阵
%用于解决PA=LU;
n = length(A);
L = zeros(n);
U = A;
p = zeros(n,1);
for k = 1:n-1
	%找索引
	MaxValue = max(abs(U(k:n,k)));
	pk = find(MaxValue == abs(U(k:n,k)));
	pk = pk(1);%可能有多个最大值只取第一个
	pk = pk + k - 1;%将子矩阵中的索引映射回原始矩阵的索引。
	p(k) = pk; 
	%交换
	[U(k,:),U(pk,:)] = deal(U(pk,:),U(k,:));
	[L(k,:),L(pk,:)] = deal(L(pk,:),L(k,:));
    L(k+1:n,k) = U(k+1:n,k)/U(k,k);
    U(k+1:n,k) = 0;
    L(k,k) = 1;
    U(k+1:n,k+1:n) = U(k+1:n,k+1:n) - L(k+1:n,k)*U(k,k+1:n);

end
L(n,n) = 1;
P = eye(n);
for k = 1:n-1
	[P(k,:),P(p(k),:)] = deal(P(p(k),:),P(k,:));
end
end