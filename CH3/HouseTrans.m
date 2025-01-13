% 07222113邓俊烽
% 算法3.2.1.Householder变换
function H = HouseTrans(x)
%H是一个n阶方阵
%用于进行Householder变换
n = length(x);
x = x/max(abs(x));%归一化处理
sigma = x(2:n)'*x(2:n);
v = x;
%如果 σ=0，说明向量x的第二个到最后一个分量都为 0，此时不需要进行变换。直接返回单位矩阵
if sigma == 0
    beta = 0;
    H = eye(n) - beta*(v*v');
	return;
end
a=sqrt(x(1)^2 + sigma);
if x(1)>0
	v(1) = -sigma/(x(1) + a); % x(1)^2 + sigma = norm2(x)^2
else
	v(1) = x(1) - a;
end
beta = 2*v(1)^2/(sigma + v(1)^2);
v = v/v(1);
H = eye(n) - beta*(v*v');
end