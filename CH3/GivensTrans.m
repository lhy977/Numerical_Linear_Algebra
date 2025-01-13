% 07222113邓俊烽
% 算法3.2.2.Givens变换
function G = GivensTrans(x,i,k)
%用于进行Givens变换
n = length(x);
%如果如果 x(k)=0，则不需要进行旋转，直接返回单位矩阵 G=I
if x(k) == 0
	c = 1;
	s = 0;
    G = eye(n);
    G(i,i) = c;
    G(i,k) = s;
    G(k,i) = -s;
    G(k,k) = c;
    return;
end
%计算 Givens 旋转参数
if abs(x(k))>abs(x(i))
    tgn = x(i)/x(k);
    s = 1/sqrt(1 + tgn^2);
    c = s*tgn;
else
    tgn = x(k)/x(i);
    c = 1/sqrt(1 + tgn^2);
    s = c*tgn;
end
%生成 Givens 矩阵 G
G = eye(n);
G(i,i) = c;
G(i,k) = s;
G(k,i) = -s;
G(k,k) = c;
end