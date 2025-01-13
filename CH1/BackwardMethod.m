% 07222113邓俊烽
% 算法1.1.2回代法
function x = BackwardMethod(U,y)
%U是一个n*n上三角矩阵,y是一个n*1的列向量
%用于解决Ux=y;
n = length(U);
x = zeros(n,1);
for j = n:-1:1
    x(j) = y(j)/U(j,j);
    y(1:j-1) = y(1:j-1) - x(j)*U(1:j-1,j);
end
end
