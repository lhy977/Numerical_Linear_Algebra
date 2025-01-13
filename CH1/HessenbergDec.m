%07222113邓俊烽
%算法6.4.1Hessenberg分解
function [Hess] = HessenbergDec(A)
n = length(A);
Hess = A;
Q = eye(n);
for k = 1:n-2
    [v,beta]= HouseTrans2(Hess(k+1:n,k));
    u=beta*Hess(k+1:n,k+1:n)*v;
    w=u-(beta*u'*v/2)*v;
    Hess(k+1:n,k:n) = H*Hess(k+1:n,k:n);
    Hess(1:n,k+1:n) = Hess(1:n,k+1:n)*H;
    Q(1:n,k+1:n) = Q(1:n,k+1:n)*H;
    Hess(k+2:n, k) = 0; % 将下三角部分的子对角线以下元素设置为 0
end
end