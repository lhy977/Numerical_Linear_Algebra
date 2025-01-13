% 07222113邓俊烽
% 算法1.1.3Gauss消去计算三角分解
function [L,U] = LUGauss(A)
%L是一个n*n下三角矩阵,U是一个n*n上三角矩阵,A是一个n*n矩阵
%解决LU=A
n = length(A);
L = eye(n); %L=I+[l_1,l_2,...,l_(n-1),0]
U = A;
for k = 1:n-1
    L(k+1:n,k) = U(k+1:n,k)/U(k,k);%这里计算的是lk
    U(k+1:n,k) = 0;%让U是一个上三角矩阵
    U(k+1:n,k+1:n) = U(k+1:n,k+1:n) - L(k+1:n,k)*U(k,k+1:n);%再继续更新U
end
end