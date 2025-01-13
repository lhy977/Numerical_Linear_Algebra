% 07222113邓俊烽
% 算法1.1.1前代法
function y = ForwardMethod(L,b)
%L是一个n*n下三角矩阵,b是一个n*1的列向量
%用于解决Ly=b
n = length(L);
y = zeros(n,1);
for j = 1:n
    y(j) = b(j)/L(j,j);
    b(j+1:n) = b(j+1:n) - y(j)*L(j+1:n,j);
end
end