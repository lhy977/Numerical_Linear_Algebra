% 07222113邓俊烽
% 算法1.3.1Cholesky分解
function L = Cholesky(A)
%L是一个n*n下三角矩阵A是一个对称正定n*n的矩阵
%用于解决A=L*L^T;
if ~all(all(A==A'))
    disp("A ~= A',不能进行Cholsky分解");
	L = 0;
    return;
end
n = length(A);
L = tril(A);
for k = 1:n
    L(k,k) = sqrt(L(k,k));
    L(k+1:n,k) = L(k+1:n,k)/L(k,k);
    for j = k+1:n
        L(j:n,j) = L(j:n,j) - L(j:n,k)*L(j,k);
    end
end
end