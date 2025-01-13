% 07222113邓俊烽
% 算法3.3.1.基于House变换的QR分解
function [Q,R] = QRDec(A)
%A是一个m>=n的m*n的矩阵,Q是一个m*m的正交矩阵,R是一个n*n的上三角矩阵
% 解决A = Q*[R;0];
[m,n] = size(A);
Q = eye(m);
R = A;
for j = 1:n
    %逐列进行 Householder 变换
	if j<m
		H = HouseTrans(R(j:m,j)); % 对 R 的第 j 列从第 j行到第 m 行进行 Householder 变换，生成 Householder 矩阵H
		R(j:m,j:n) = H*R(j:m,j:n);%将 H 应用到R 的子矩阵上，使得 R(j+1:m,j) 变为 0。
		H2 = eye(m); % H2 = diag(I(j-1))+Hj, such that size H2 is (m,m)
		H2(j:m,j:m) = H;%将 H 嵌入到 H2 的右下角，使得 H2 是一个 m×m 的矩阵。
		Q = Q*H2;
		R(j+1:m,j) = 0;
	end
end
R = R(1:n,:);
end