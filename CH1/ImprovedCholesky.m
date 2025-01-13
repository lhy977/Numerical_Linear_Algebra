% 07222113邓俊烽
% 算法1.3.2改进的Cholesky分解
function [L,D] = ImprovedCholesky(A)
%L是一个n*n单位下三角矩阵，A是一个对称正定n*n的矩阵，D是一个对角矩阵
%用于解决A=L*D*L^T;
if ~all(all(A==A'))
    disp("A ~= A',不能进行Cholsky分解");
	D = 0;
	L = 0;
    return;
end
n = length(A);
L = eye(n) + tril(A,-1);%生成单位下三角
D = diag(A);
%D转换为对角矩阵。
for k = 1:n-1
    if k == 1
        L(k+1:n,k) = L(k+1:n,k)/D(k);
    else
        L(k+1:n,k) = (L(k+1:n,k) - L(k+1:n,1:k-1)*v)/D(k);
    end
    v = L(k+1,1:k)'.*D(1:k);
    D(k+1) = D(k+1) - L(k+1,1:k)*v;
end
D = diag(D);
end